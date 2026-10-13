/* LP2 encoder: bit counting, serialisation (pack_unit) and the library call
 * (LP2Encoder.encode_call).  Mirrors spec/reference/lp2_bitstream.py and the
 * LP2Encoder class of atrac3_encoder_lp2.py. */
#include "internal.h"

static const int QBITS[8] = {0, 2, 3, 3, 4, 4, 5, 6};
static const int TONE_NVALS[8] = {1, 2, 3, 4, 5, 6, 7, 8};

int lp2_mantissa_bits(int mode, int sel, const int *q, int n)
{
    int b = QBITS[sel], m = (1 << b) - 1, tot = 0;
    const at3_code *L = at3_codes[mode][sel];
    if (sel != 1) {
        for (int k = 0; k < n; k++)
            tot += L[q[k] & m].len;
    } else {
        for (int k = 0; k + 1 < n; k += 2)
            tot += L[((q[k] & m) << b) | (q[k + 1] & m)].len;
    }
    return tot;
}

static int component_bits(const lp2_unit *u)
{
    int n = u->ngroups ? 7 : 5;
    for (int gi = 0; gi < u->ngroups; gi++) {
        const lp2_group *G = &u->groups[gi];
        n += u->nbands + 6;
        for (int c = 0; c < 4 * u->nbands; c++) {
            if (!G->flags[c >> 2])
                continue;
            n += 3;
            for (int k = 0; k < G->ncell[c]; k++) {
                const lp2_comp *cp = &G->cells[c][k];
                int cnt = imin(TONE_NVALS[cp->len_code], 1024 - cp->pos);
                n += 12 + lp2_mantissa_bits(cp->mode, cp->sel, cp->mant, cnt);
            }
        }
    }
    return n;
}

static int spectrum_bits(const lp2_unit *u)
{
    int n = 3 * u->nsub + 6;
    for (int i = 0; i < u->nsub; i++)
        if (u->sel[i] > 0)
            n += 6 + lp2_mantissa_bits(u->mode, u->sel[i], u->q + at3_sb_start[i],
                                  at3_sb_start[i + 1] - at3_sb_start[i]);
    return n;
}

int lp2_gain_bits(const lp2_gain *gain, int nbands)
{
    int pts = 0;
    for (int b = 0; b < nbands; b++)
        pts += gain[b].n;
    return 3 * (nbands + 3 * pts);
}

int lp2_packed_bits(const lp2_unit *u)
{
    return component_bits(u) + spectrum_bits(u) + lp2_gain_bits(u->gain, u->nbands) + 8;
}

static void pack_mantissas(bitwriter *bw, int mode, int sel, const int *q, int n)
{
    int b = QBITS[sel], m = (1 << b) - 1;
    const at3_code *tab = at3_codes[mode][sel];
    if (sel != 1) {
        for (int k = 0; k < n; k++)
            bw_put(bw, tab[q[k] & m].code, tab[q[k] & m].len);
    } else {
        for (int k = 0; k + 1 < n; k += 2) {
            const at3_code *c = &tab[((q[k] & m) << b) | (q[k + 1] & m)];
            bw_put(bw, c->code, c->len);
        }
    }
}

void lp2_pack_unit(const lp2_unit *u, uint8_t *out)
{
    bitwriter bw = {out, LP2_UNIT_BITS, 0};
    memset(out, 0, LP2_UNIT_BYTES);
    bw_put(&bw, 0x28, 6);
    bw_put(&bw, u->nbands - 1, 2);
    for (int b = 0; b < u->nbands; b++) {
        bw_put(&bw, u->gain[b].n, 3);
        for (int k = 0; k < u->gain[b].n; k++) {
            bw_put(&bw, u->gain[b].lev[k], 4);
            bw_put(&bw, u->gain[b].loc[k], 5);
        }
    }
    bw_put(&bw, u->ngroups, 5);
    if (u->ngroups)
        bw_put(&bw, u->tone_mode, 2);
    for (int gi = 0; gi < u->ngroups; gi++) {
        const lp2_group *G = &u->groups[gi];
        for (int b = 0; b < u->nbands; b++)
            bw_put(&bw, G->flags[b], 1);
        bw_put(&bw, G->len_code, 3);
        bw_put(&bw, G->sel, 3);
        for (int c = 0; c < 4 * u->nbands; c++) {
            if (!G->flags[c >> 2])
                continue;
            bw_put(&bw, G->ncell[c], 3);
            for (int k = 0; k < G->ncell[c]; k++) {
                const lp2_comp *cp = &G->cells[c][k];
                bw_put(&bw, cp->sfi, 6);
                bw_put(&bw, cp->pos & 63, 6);
                pack_mantissas(&bw, G->mode, G->sel, cp->mant, imin(TONE_NVALS[G->len_code], 0x400 - cp->pos));
            }
        }
    }
    bw_put(&bw, u->nsub - 1, 5);
    bw_put(&bw, u->mode, 1);
    for (int i = 0; i < u->nsub; i++)
        bw_put(&bw, u->sel[i], 3);
    for (int i = 0; i < u->nsub; i++)
        if (u->sel[i] > 0)
            bw_put(&bw, u->sfi[i], 6);
    for (int i = 0; i < u->nsub; i++)
        if (u->sel[i] > 0)
            pack_mantissas(&bw, u->mode, u->sel[i], u->q + at3_sb_start[i],
                       at3_sb_start[i + 1] - at3_sb_start[i]);
}

/* ================================================================ library call */
#define END_THRESHOLD 0x376      /* a call whose block had fewer samples ends the stream two calls later */

void lp2_encoder_init(lp2_encoder *e, int channels)
{
    memset(e, 0, sizeof *e);
    e->channels = channels;
    for (int c = 0; c < channels; c++)
        lp2_analysis_init(&e->ana[c]);
    e->counts[0] = e->counts[1] = e->counts[2] = AT3_FRAME_SAMPLES;
}

/* block: interleaved int16, nsamples <= 1024 (zero-padded internally).  The frame of call k
 * describes the input of call k-1 (one block of internal delay); call 0 yields no frame. */
int lp2_encoder_call(lp2_encoder *e, const int16_t *block, int nsamples, uint8_t *out, int *end)
{
    float pcm[AT3_FRAME_SAMPLES];
    e->counts[2] = e->counts[1];
    e->counts[1] = e->counts[0];
    e->counts[0] = nsamples;
    *end = e->counts[2] < END_THRESHOLD;
    for (int c = 0; c < e->channels; c++) {
        lp2_analysis *a = &e->ana[c];
        for (int i = 0; i < AT3_FRAME_SAMPLES; i++)
            pcm[i] = i < nsamples ? (float)block[i * e->channels + c] : 0.0f;
        lp2_analysis_process(a, pcm, e->spec_gc, e->spec_plain);
        lp2_encode_unit(&e->unit, e->spec_gc, e->spec_plain, a->cur_gain, a->prev_gain);
        lp2_pack_unit(&e->unit, out + c * LP2_UNIT_BYTES);
    }
    return e->calls++ > 0;
}
