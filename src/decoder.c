/* ATRAC3 frame decoder.  Mirrors spec/reference/atrac3_decoder.py (double precision
 * throughout); output matches the reference decoder apart from summation order. */
#include <math.h>
#include <stdlib.h>
#include "internal.h"

static const int QBITS[8] = {0, 2, 3, 3, 4, 4, 5, 6};
static const double MAX_QUANT[8] = {0.0, 1.5, 2.5, 3.5, 4.5, 7.5, 15.5, 31.5};
/* joint stereo: matrix selector -> coefficients of the channel-0 and channel-1 signals */
static const double MATRIX[4][2] = {{0.0, 1.0}, {1.0, 1.0}, {0.0, 0.0}, {0.5, 0.5}};

typedef struct { int n; int lev[8], loc[8]; } dec_gain;

typedef struct {
    dec_gain prev_gain[4];
    double overlap[4][256];
    double qmf_delay[3][46];
} dec_unit;

typedef struct { int w; int mat[4]; int last_band; } js_params;   /* w < 0: none */

typedef struct { uint8_t len, idx; } lut_entry;

struct at3_decoder {
    int ch, fb, js;
    dec_unit units[2];
    js_params js_prev, js_cur, js_next;
    lut_entry lut[2][8][256];               /* next 8 bits -> code length and value index */
    double imdct[256][512];                 /* transposed: imdct[k][n] */
    double spec[2][1024];
    double bands[2][4][256];
};

/* ---------------------------------------------------------------- bit reader */
typedef struct { const uint8_t *data; int nbytes; int pos; } bitreader;

static int br_bit(const bitreader *br, int pos)
{
    return (pos >> 3) < br->nbytes ? (br->data[pos >> 3] >> (7 - (pos & 7))) & 1 : 0;
}

static unsigned br_read(bitreader *br, int n)
{
    unsigned v = 0;
    for (int i = 0; i < n; i++)
        v = (v << 1) | br_bit(br, br->pos++);
    return v;
}

static int sext(int v, int bits) { return v >= 1 << (bits - 1) ? v - (1 << bits) : v; }

/* Read `count` mantissas for a quantiser selector (1..7) and coding mode. */
static int read_mantissas(at3_decoder *d, bitreader *br, int sel, int mode, int count, int *out)
{
    int n = 0;
    while (n < count) {
        unsigned peek = 0;
        for (int i = 0; i < 8; i++)
            peek = (peek << 1) | br_bit(br, br->pos + i);
        const lut_entry *e = &d->lut[mode][sel][peek];
        if (!e->len)
            return -1;
        br->pos += e->len;
        if (sel == 1) {
            out[n++] = sext(e->idx >> 2, 2);
            if (n < count)
                out[n++] = sext(e->idx & 3, 2);
        } else {
            out[n++] = sext(e->idx, QBITS[sel]);
        }
    }
    return 0;
}

/* Parse one sound unit into its gain data and spectrum.  Returns nbands, or -1. */
static int decode_unit(at3_decoder *d, bitreader *br, int js_second, dec_gain *gain, double *spec)
{
    int m[1024];
    if (js_second) {
        if (br_read(br, 2) != 3)
            return -1;
    } else if (br_read(br, 6) != 0x28) {
        return -1;
    }
    int nbands = br_read(br, 2) + 1;
    for (int b = 0; b < 4; b++) {
        gain[b].n = 0;
        if (b >= nbands)
            continue;
        gain[b].n = br_read(br, 3);
        for (int k = 0; k < gain[b].n; k++) {
            gain[b].lev[k] = br_read(br, 4);
            gain[b].loc[k] = br_read(br, 5);
            if (k > 0 && gain[b].loc[k] <= gain[b].loc[k - 1])
                return -1;
        }
    }
    memset(spec, 0, 1024 * sizeof(double));
    int ngroups = br_read(br, 5);
    if (ngroups) {                           /* tonal components */
        int sel_mode = br_read(br, 2);
        if (sel_mode == 2)
            return -1;
        for (int gi = 0; gi < ngroups; gi++) {
            int flags[4];
            for (int b = 0; b < nbands; b++)
                flags[b] = br_read(br, 1);
            int nvalues = br_read(br, 3) + 1;
            int qstep = br_read(br, 3);
            if (qstep <= 1)
                return -1;
            int mode = sel_mode == 3 ? (int)br_read(br, 1) : sel_mode;
            for (int b = 0; b < nbands; b++) {
                if (!flags[b])
                    continue;
                for (int g = 0; g < 4; g++) {
                    int ncomp = br_read(br, 3);
                    for (int c = 0; c < ncomp; c++) {
                        int sfi = br_read(br, 6);
                        int pos = b * 256 + g * 64 + br_read(br, 6);
                        int cnt = imin(nvalues, 1024 - pos);
                        if (read_mantissas(d, br, qstep, mode, cnt, m))
                            return -1;
                        double scale = dec_sf[sfi] / MAX_QUANT[qstep];
                        for (int k = 0; k < cnt; k++)
                            spec[pos + k] += m[k] * scale;
                    }
                }
            }
        }
    }
    int nsub = br_read(br, 5) + 1;           /* spectrum */
    int mode = br_read(br, 1);
    int sels[32], sfs[32];
    for (int i = 0; i < nsub; i++)
        sels[i] = br_read(br, 3);
    for (int i = 0; i < nsub; i++)
        sfs[i] = sels[i] ? (int)br_read(br, 6) : 0;
    for (int i = 0; i < nsub; i++) {
        if (!sels[i])
            continue;
        int a = at3_sb_start[i], n = at3_sb_start[i + 1] - a;
        if (read_mantissas(d, br, sels[i], mode, n, m))
            return -1;
        double scale = dec_sf[sfs[i]] / MAX_QUANT[sels[i]];
        for (int k = 0; k < n; k++)
            spec[a + k] += m[k] * scale;
    }
    return nbands;
}

/* ---------------------------------------------------------------- synthesis */
static void gain_curve(const dec_gain *gp, double *g)
{
    for (int i = 0; i < 256; i++)
        g[i] = 1.0;
    int pos = 0;
    for (int i = 0; i < gp->n; i++) {
        int lev = gp->lev[i], start = gp->loc[i] * 8;
        int next = i + 1 < gp->n ? gp->lev[i + 1] : 4;
        double cur = ldexp(1.0, 4 - lev);
        for (int k = pos; k < start; k++)
            g[k] = cur;
        for (int k = 0; k < 8; k++)
            if (start + k < 256)
                g[start + k] = cur * dec_inc_pow[lev - next + 15][k];
        pos = start + 8;
    }
}

static void synth_band(at3_decoder *d, dec_unit *u, int band, const double *coefs,
                       const dec_gain *gain_now, const dec_gain *gain_next, double *out)
{
    double y[512] = {0}, g[256];
    for (int k = 0; k < 256; k++) {
        double c = coefs[band & 1 ? 255 - k : k];
        if (c == 0.0)
            continue;
        const double *row = d->imdct[k];
        for (int n = 0; n < 512; n++)
            y[n] += row[n] * c;
    }
    for (int n = 0; n < 512; n++)
        y[n] *= dec_window[n];
    double scale_next = gain_next->n ? ldexp(1.0, 4 - gain_next->lev[0]) : 1.0;
    gain_curve(gain_now, g);
    for (int i = 0; i < 256; i++)
        out[i] = (y[i] * scale_next + u->overlap[band][i]) * g[i];
    memcpy(u->overlap[band], y + 256, 256 * sizeof(double));
}

static void iqmf(const double *lo, const double *hi, int n, double *delay, double *out)
{
    double tmp[46 + 1024];
    memcpy(tmp, delay, 46 * sizeof(double));
    for (int j = 0; j < n; j++) {
        tmp[46 + 2 * j] = lo[j] + hi[j];
        tmp[47 + 2 * j] = lo[j] - hi[j];
    }
    for (int j = 0; j < n; j++) {
        const double *p = tmp + 2 * j;
        double a = 0.0, b = 0.0;
        for (int i = 0; i < 24; i++) {
            a += p[2 * i + 1] * dec_qmf_win[2 * i + 1];
            b += p[2 * i] * dec_qmf_win[2 * i];
        }
        out[2 * j] = a;
        out[2 * j + 1] = b;
    }
    memcpy(delay, tmp + 2 * n, 46 * sizeof(double));
}

/* Inverse matrixing + channel weighting on the QMF-band time signals (in place). */
static void apply_js(at3_decoder *d)
{
    const js_params *P = &d->js_prev, *C = &d->js_cur;
    for (int b = 0; b < 4; b++) {
        double *L = d->bands[0][b], *R = d->bands[1][b];
        double p0 = MATRIX[P->mat[b]][0], p1 = MATRIX[P->mat[b]][1];
        double c0 = MATRIX[C->mat[b]][0], c1 = MATRIX[C->mat[b]][1];
        double pw0 = 2.0, pw1 = 2.0, cw0 = 2.0, cw1 = 2.0;
        if (b > P->last_band) {
            pw0 = dec_weights[P->w][0];
            pw1 = dec_weights[P->w][1];
        }
        if (b > C->last_band) {
            cw0 = dec_weights[C->w][0];
            cw1 = dec_weights[C->w][1];
        }
        int start = 0;
        if (pw0 != cw0 || P->mat[b] != C->mat[b]) {
            for (int n = 0; n < 8; n++) {
                double t = (n - 8) / 8.0;
                double s = (c0 + t * (c0 - p0)) * L[n] + (c1 + t * (c1 - p1)) * R[n];
                R[n] = (L[n] - s) * (cw1 + t * (cw1 - pw1));
                L[n] = (cw0 + t * (cw0 - pw0)) * s;
            }
            start = 8;
        }
        for (int n = start; n < 256; n++) {
            double s = c0 * L[n] + c1 * R[n];
            R[n] = (L[n] - s) * cw1;
            L[n] = s * cw0;
        }
    }
    d->js_prev = d->js_cur;
    d->js_cur = d->js_next;
}

/* ---------------------------------------------------------------- public API */
at3_decoder *at3_decoder_create(int channels, int frame_bytes, int joint_stereo)
{
    if (channels < 1 || channels > 2 || frame_bytes < 8 || frame_bytes > 4096 ||
        (joint_stereo && channels != 2))
        return NULL;
    at3_decoder *d = calloc(1, sizeof *d);
    if (!d)
        return NULL;
    d->ch = channels;
    d->fb = frame_bytes;
    d->js = joint_stereo;
    js_params init = {-1, {3, 3, 3, 3}, 4};
    d->js_prev = d->js_cur = d->js_next = init;
    for (int mode = 0; mode < 2; mode++)
        for (int sel = 1; sel < 8; sel++)
            for (int idx = 0; idx < 64; idx++) {
                const at3_code *c = &at3_codes[mode][sel][idx];
                if (!c->len)
                    continue;
                int sh = 8 - c->len;
                for (int r = 0; r < 1 << sh; r++) {
                    lut_entry *e = &d->lut[mode][sel][(c->code << sh) | r];
                    e->len = c->len;
                    e->idx = (uint8_t)idx;
                }
            }
    const double pi = 3.14159265358979323846;
    for (int k = 0; k < 256; k++)
        for (int n = 0; n < 512; n++)
            d->imdct[k][n] = cos(pi / 256.0 * (n + 0.5 + 128) * (k + 0.5));
    return d;
}

void at3_decoder_destroy(at3_decoder *d)
{
    free(d);
}

int at3_decoder_decode_frame(at3_decoder *d, const uint8_t *data, int16_t *pcm)
{
    dec_gain gain[2][4];
    int per = d->fb / d->ch;
    if (!d->js) {
        for (int c = 0; c < d->ch; c++) {
            bitreader br = {data + c * per, per, 0};
            if (decode_unit(d, &br, 0, gain[c], d->spec[c]) < 0)
                return -1;
        }
    } else {
        uint8_t rev[4096];
        bitreader br = {data, d->fb, 0};
        if (decode_unit(d, &br, 0, gain[0], d->spec[0]) < 0)
            return -1;
        for (int i = 0; i < d->fb; i++)
            rev[i] = data[d->fb - 1 - i];
        int k = 0;
        while (k < d->fb && rev[k] == 0xF8)
            k++;
        bitreader br1 = {rev + k, d->fb - k, 0};
        int wflag = br_read(&br1, 1), widx = br_read(&br1, 3);
        int mat[4];
        for (int b = 0; b < 4; b++)
            mat[b] = br_read(&br1, 2);
        int nb = decode_unit(d, &br1, 1, gain[1], d->spec[1]);
        if (nb < 0)
            return -1;
        for (int b = 0; b < 4; b++)
            if (mat[b] == 2)
                return -1;
        d->js_next.w = wflag * 8 + widx;
        memcpy(d->js_next.mat, mat, sizeof mat);
        d->js_next.last_band = nb - 1;
    }
    for (int c = 0; c < d->ch; c++) {
        dec_unit *u = &d->units[c];
        for (int b = 0; b < 4; b++)
            synth_band(d, u, b, d->spec[c] + 256 * b, &u->prev_gain[b], &gain[c][b], d->bands[c][b]);
        memcpy(u->prev_gain, gain[c], sizeof u->prev_gain);
    }
    if (d->js)
        apply_js(d);
    for (int c = 0; c < d->ch; c++) {
        double low[512], high[512], out[1024];
        double (*dl)[46] = d->units[c].qmf_delay;
        iqmf(d->bands[c][0], d->bands[c][1], 256, dl[0], low);
        iqmf(d->bands[c][3], d->bands[c][2], 256, dl[1], high);
        iqmf(low, high, 512, dl[2], out);
        for (int i = 0; i < 1024; i++) {
            double v = out[i];
            if (v > 32767.0)
                v = 32767.0;
            if (v < -32767.0)
                v = -32767.0;
            pcm[i * d->ch + c] = (int16_t)nearbyint(v);
        }
    }
    return 0;
}
