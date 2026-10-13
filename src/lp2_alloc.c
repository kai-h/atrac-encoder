/* LP2 encoder core: tonal component extraction, bit allocation and quantisation of one sound
 * unit (encode_unit and its helpers).  Mirrors spec/reference/lp2_alloc.py. */
#include <math.h>
#include "internal.h"

#define NSUB_DEFAULT 29      /* coded_layout(15100 Hz): lines 0..703 */
#define NBANDS_DEFAULT 3

static const int MAXQ[8] = {0, 1, 2, 3, 4, 7, 15, 31};
static const double T2[8] = {1.0, 3.0, 5.0, 7.0, 9.0, 15.0, 31.0, 63.0};
static const int TONE_NVALS[8] = {1, 2, 3, 4, 5, 6, 7, 8};
static const int EST_BITS_CLC[8] = {0, 2, 3, 3, 4, 4, 5, 6};
enum { TONE_SEL = 7, TONE_LEN_CODE = 3, TONE_MODE = 1 };

static double offset_divisor(int i)
{
    if (i < 8) return 0.3333333333333333;
    if (i < 12) return 0.30303030303030304;
    if (i < 16) return 0.29411764705882354;
    if (i < 18) return 0.2857142857142857;
    if (i < 26) return 0.2777777777777778;
    if (i < 28) return 0.2631578947368421;
    if (i < 30) return 0.25;
    return 0.23809523809523808;
}

static double offset_step_weight(int i)      /* < 0: added unweighted */
{
    if (i < 1) return 0.2;
    if (i < 2) return 0.3;
    if (i < 8) return 0.4;
    if (i < 18) return 0.6;
    return -1.0;
}

static int dead_zone_limit(double density)
{
    if (!(density <= 1.15))
        return 0;
    double d = (double)(float)density;
    if (d > 1.1) return 3;
    if (d > 1.05) return 5;
    if (d > 1.0) return 7;
    if (d > 0.95) return 9;
    if (d > 0.9) return 11;
    return 12;
}

/* ================================================================ scale factors */
static int sf_index(float v)
{
    if (!(at3_sf[0] <= v))
        return 0;
    int e = 32, c = 16;
    while (c > 0) {
        e = at3_sf[e] <= v ? e + c : e - c;
        c >>= 1;
    }
    if (e < 63 && at3_sf[e] <= v)
        e++;
    return e;
}

static int sf_index_tone(float v)
{
    v = fabsf(v);
    for (int i = 0; i < 63; i++)
        if (v < at3_sf[i])
            return i;
    return 63;
}

static float abs_max(const float *a, int n)
{
    float m = fabsf(a[0]);
    for (int i = 1; i < n; i++)
        if (fabsf(a[i]) > m)
            m = fabsf(a[i]);
    return m;
}

static void group_sf(const float *spec, int *G)
{
    for (int g = 0; g < 256; g++)
        G[g] = sf_index(abs_max(spec + 4 * g, 4));
}

static void subband_sf(const int *G, int nsub, int *S)
{
    for (int i = 0; i < nsub; i++) {
        int m = G[at3_sb_start[i] / 4];
        for (int g = at3_sb_start[i] / 4 + 1; g < at3_sb_start[i + 1] / 4; g++)
            m = imax(m, G[g]);
        S[i] = m;
    }
}

/* ================================================================ quantisation */
static float dead_zone(int e, int sel)
{
    if (e == 0)
        return 0.0f;
    return (float)((double)lp2_t1[imin(imax(e, 0), 12)] / T2[sel]);
}

/* quantise_subband: returns bits including the 6-bit scale factor; q untouched if sel 0. */
static int quantise_subband(int sel, float thr, const float *x, int n, int *q)
{
    if (sel == 0)
        return 0;
    int mq = MAXQ[sel];
    double s = (double)(float)(mq + 0.5);
    for (int k = 0; k < n; k++) {
        int v = (int)(s * (double)x[k] + 31.5) - 31;
        v = imin(imax(v, -mq), mq);
        if (!(thr <= fabsf(x[k])))
            v = 0;
        q[k] = v;
    }
    return lp2_mantissa_bits(0, sel, q, n) + 6;
}

static int quant_sb(const float *A, int i, int sel, int e, int *q)
{
    int a = at3_sb_start[i];
    return quantise_subband(sel, dead_zone(e, sel), A + a, at3_sb_start[i + 1] - a, q);
}

static int tone_mantissa(float x, float sfval, int mq)
{
    int v = (int)(((double)mq + 0.5) * ((double)x / (double)sfval) + 31.5) - 31;
    return imin(imax(v, -mq), mq);
}

/* ================================================================ tonal components */
static int quantise_component(const float *spec, lp2_comp *c)
{
    int n = TONE_NVALS[c->len_code];
    c->sfi = sf_index_tone(abs_max(spec + c->pos, imin(n, 1024 - c->pos)));
    float sfv = at3_sf[c->sfi];
    int mq = MAXQ[c->sel];
    for (int k = 0; k < n; k++)
        c->mant[k] = c->pos + k <= 1023 ? tone_mantissa(spec[c->pos + k], sfv, mq) : 0;
    return lp2_mantissa_bits(c->mode, c->sel, c->mant, n) + 12;
}

static void subtract_component(float *spec, const lp2_comp *c)
{
    float step = (float)((double)at3_sf[c->sfi] / ((double)MAXQ[c->sel] + 0.5));
    for (int k = 0; k < TONE_NVALS[c->len_code]; k++) {
        int p = c->pos + k;
        if (p <= 0x3ff)
            spec[p] = (float)((double)spec[p] - (double)c->mant[k] * (double)step);
    }
}

static lp2_group *new_group(lp2_unit *u)
{
    lp2_group *G = &u->groups[u->ngroups++];
    memset(G, 0, sizeof *G);
    G->sel = TONE_SEL;
    G->len_code = TONE_LEN_CODE;
    G->mode = TONE_MODE;
    return G;
}

static void new_component(lp2_comp *c, const lp2_group *G, int g, const int *GA)
{
    memset(c, 0, sizeof *c);
    c->pos = 4 * g;
    c->sfi = GA[g];
    c->len_code = G->len_code;
    c->sel = G->sel;
    c->mode = G->mode;
}

/* single_tone_class: 1 or 2 if the strongest 4-line group stands >= 17 scale steps above
 * everything except +-2 groups around it and its QMF alias image. */
static int single_tone_class(const int *GA)
{
    int p = 0;
    for (int g = 0; g < 256; g++)
        if (GA[p] < GA[g])
            p = g;
    int r = p & 63, blk = p >> 6;
    int bands[2] = {blk, r < 32 ? blk - 1 : blk + 1};
    int centre[2] = {r, 63 - r};
    int mx = 0;
    for (int b = 0; b < 4; b++) {
        int e = bands[0] == b ? 0 : bands[1] == b ? 1 : -1;
        int base = 64 * b;
        if (e < 0) {
            for (int i = 0; i < 64; i++)
                mx = imax(mx, GA[base + i]);
            continue;
        }
        int u = centre[e];
        for (int i = 0; i < imax(u - 2, 0); i++)
            mx = imax(mx, GA[base + i]);
        if (u + 2 < 63)
            for (int i = u + 3; i < 64; i++)
                mx = imax(mx, GA[base + i]);
    }
    int v = GA[p];
    if (mx < imax(v - 17, 9)) {
        if (v >= 36) return 2;
        if (v >= 21) return 1;
    }
    return 0;
}

/* single_tone_path: code the lines around a dominant peak (and its alias image in the
 * neighbouring QMF band) as tonal components, in ngroups_tone groups. */
static int single_tone_path(int budget, int ngroups_tone, int width, int peak, int nbands,
                                float *spec, int *GA, lp2_unit *u)
{
    int blk = peak >> 6, r = peak & 63;
    int bands[2] = {blk, r < 32 ? blk - 1 : blk + 1};
    int centre[2] = {r, 63 - r};
    int bits = (nbands + 6) * ngroups_tone;
    int go_on = 1;
    for (int t = 0; t < ngroups_tone; t++) {
        lp2_group *G = new_group(u);
        if (!go_on)
            continue;
        for (int e = 0; e < 2; e++) {
            int b = bands[e];
            if (b < 0 || b >= nbands)
                continue;
            if (!G->flags[b]) {
                G->flags[b] = 1;
                bits += 12;
            }
            if (bits >= budget) {
                G->flags[b] = 0;
                bits -= 12;
                go_on = 1;
                continue;
            }
            int c0 = centre[e] + 64 * b;
            for (int g = imax(c0 - width, 0); g <= imin(c0 + width, 255); g++) {
                int cell = g >> 4;
                if (G->ncell[cell] >= 7)
                    continue;
                lp2_comp c;
                new_component(&c, G, g, GA);
                int cb = quantise_component(spec, &c);
                if (bits + cb > budget) {
                    go_on = 0;
                    continue;
                }
                u->ncomp++;
                G->cells[cell][G->ncell[cell]++] = c;
                subtract_component(spec, &c);
                GA[g] = sf_index(abs_max(spec + 4 * g, 4));
                bits += cb;
            }
        }
    }
    return bits;
}

/* multitone_path: up to two passes; every 4-line group whose plain-spectrum scale factor
 * exceeds mean + 32 becomes a component (max 7 per cell, 64 in total). */
static int multitone_path(int budget, int ngroups, int nbands, float mean, int *GA, int *GB,
                             float *spec, lp2_unit *u)
{
    int bits = 0, stop = 0;
    double thr = (double)(float)((double)mean + 32.0);
    for (int pass = 0; pass < 2; pass++) {
        if (stop)
            continue;
        int cnt[16] = {0}, nc = 0;
        for (int g = 0; g < ngroups; g++)
            if ((double)GB[g] > thr && cnt[g >> 4] < 7 && u->ncomp + nc < 64) {
                cnt[g >> 4]++;
                nc++;
            }
        if (nc <= 0)
            continue;
        if (nbands + bits + 6 > budget) {
            stop = 1;
            continue;
        }
        bits += nbands + 6;
        lp2_group *G = new_group(u);
        for (int g = 0; g < ngroups; g++) {
            int cell = g >> 4, band = g >> 6;
            if (!((double)GB[g] > thr && G->ncell[cell] < 7 && u->ncomp < 64))
                continue;
            if (!G->flags[band]) {
                if (bits + 12 <= budget) {
                    bits += 12;
                    G->flags[band] = 1;
                } else {
                    stop = 1;
                }
            }
            if (G->flags[band]) {
                if (bits + 4 * EST_BITS_CLC[G->sel] + 12 > budget) {
                    stop = 1;
                    continue;
                }
                lp2_comp *c = &G->cells[cell][G->ncell[cell]++];
                new_component(c, G, g, GA);
                u->ncomp++;
                bits += quantise_component(spec, c);
                subtract_component(spec, c);
                int s = sf_index(abs_max(spec + 4 * g, 4));
                GB[g] += s - GA[g];
                GA[g] = s;
            }
        }
    }
    return bits;
}

/* ================================================================ word-length decisions */
static int is_attack(const lp2_gain *g)
{
    int mx = 4, mn = 4;
    for (int k = 0; k < g->n; k++) {
        mx = imax(mx, g->lev[k]);
        mn = imin(mn, g->lev[k]);
    }
    return mx - mn > 2;
}

/* word_lengths: word lengths from offsets, then zero the subbands that are absolutely
 * quiet or masked by one of the 8 subbands below them.  Returns the quiet threshold. */
static int word_lengths(int slope, const float *offs, const int *SB, int n, int *W)
{
    for (int i = 0; i < n; i++)
        W[i] = imin(imax((int)((double)offs[i] + 0.5), 1), 7);
    int mx = 0;
    if (n > 1) {
        mx = SB[1];
        for (int i = 2; i < n; i++)
            mx = imax(mx, SB[i]);
    }
    mx = imax(mx, 0);
    int quiet = (n > 1 && mx >= 30) ? 6 : imax(mx / 6, 1);
    for (int i = 0; i < n; i++) {
        if (SB[i] < quiet)
            W[i] = 0;
        for (int c = 0; c < 8; c++) {
            int j = i - 1 - c;
            if (j >= 0 && SB[i] < SB[j] - slope)
                W[i] = 0;
        }
    }
    return quiet;
}

/* order_descending: Shell sort (gaps ..., 13, 4, 1), descending; returns indices. */
static void order_descending(const int *keys, int n, int *order)
{
    int K[32];
    for (int i = 0; i < n; i++) {
        K[i] = keys[i];
        order[i] = i;
    }
    int h = 1;
    while (n >= 1 && h <= n)
        h = 3 * h + 1;
    for (;;) {
        h /= 3;
        if (h < 1)
            return;
        if (h >= n)
            continue;
        for (int t = 0; t < n - h; t++) {
            int v = K[h + t], j = t;
            while (j >= 0 && K[j] < v) {
                K[j + h] = K[j];
                int tmp = order[j + h];
                order[j + h] = order[j];
                order[j] = tmp;
                j -= h;
            }
            K[j + h] = v;
        }
    }
}

/* ================================================================ encode_unit */
typedef struct {
    const float *A;
    int nsub, bits, total, used;
    int *W, *E, *sbits, *Q, *sel;
} alloc_ctx;

static int count_bits(alloc_ctx *c, const int *flags)
{
    int tot = 0;
    for (int i = 0; i < c->nsub; i++) {
        if (flags[i])
            c->sbits[i] = quant_sb(c->A, i, c->W[i], c->E[i], c->Q + at3_sb_start[i]);
        tot += c->sbits[i];
    }
    return tot;
}

static void commit(alloc_ctx *c, int i, int b_new, const int *q)
{
    c->total += b_new - c->sbits[i];
    c->used += b_new - c->sbits[i];
    c->sbits[i] = b_new;
    memcpy(c->Q + at3_sb_start[i], q, (at3_sb_start[i + 1] - at3_sb_start[i]) * sizeof(int));
}

static int try_raise(alloc_ctx *c, int i)
{
    int q[128];
    c->sel[i]++;
    int b_new = quant_sb(c->A, i, c->sel[i], c->E[i], q);
    if (c->used + b_new - c->sbits[i] > c->bits) {
        c->sel[i]--;
        return 0;
    }
    commit(c, i, b_new, q);
    return 1;
}

void lp2_encode_unit(lp2_unit *u, const float *spec_gc, const float *spec_plain,
                       const lp2_gain *gain, const lp2_gain *prev_gain)
{
    const int bits = LP2_UNIT_BITS;
    int nsub = NSUB_DEFAULT, nbands = NBANDS_DEFAULT;
    float A[1024];               /* modified: tones removed, then normalised */
    int GA[256], GB[256], SA[32], SB[32], W[32], W0[32], W_old[32], E[32] = {0}, sbits[32] = {0};
    int fixed[32] = {0}, flags[32], keys[32], order[32], rank[32], q[128];
    float off[32], offs[32], margin[32];
    int *Q = u->q;

    memcpy(A, spec_gc, sizeof A);
    memset(u, 0, sizeof *u);
    u->gain = gain;
    group_sf(A, GA);
    group_sf(spec_plain, GB);
    int peak = 0;
    for (int g = 0; g < 256; g++)
        if (GB[peak] < GB[g])
            peak = g;

    /* ---- 1. tonal components and the spectrum bit budget */
    int tone_class = single_tone_class(GA);
    int tonal = 0, tone_bits = 0, budget;
    if (tone_class > 0 && nsub > 25) {
        tonal = 1;
        if (peak >= 160) {
            nsub = 32;
            nbands = 4;
        } else if (peak >= 96) {
            nsub = 30;
            nbands = 3;
        }
        budget = bits - 3 * nsub - 8 - lp2_gain_bits(gain, nbands) - 13;
        tone_bits = single_tone_path(budget, tone_class, (bits > 1200) + 1, peak, nbands,
                                         A, GA, u);
        budget -= tone_bits;
    } else {
        budget = bits - 3 * nsub - 8 - lp2_gain_bits(gain, nbands) - 13;
    }
    int ngroups = at3_sb_start[nsub] / 4;
    int gsum = 0;
    for (int g = 0; g < ngroups; g++)
        gsum += GB[g];
    float mean = (float)((double)gsum / (double)ngroups);
    int attack = -1;                         /* 0: attack now, 1: attack in previous frame */
    for (int b = 0; b < nbands; b++)
        if (is_attack(&gain[b]))
            attack = 0;
    if (attack < 0)
        for (int b = 0; b < nbands; b++)
            if (is_attack(&prev_gain[b]))
                attack = 1;
    int slope = tonal ? 60 : 30;
    int multi_bits = 0;
    if (tone_bits == 0)
        multi_bits = multitone_path(budget, ngroups, nbands, mean, GA, GB, A, u);
    budget -= multi_bits;
    if (!u->ngroups)
        budget += 2;                         /* no 2-bit tonal mode selector */
    subband_sf(GA, nsub, SA);
    subband_sf(GB, nsub, SB);

    /* ---- 2. initial word lengths from the scale-factor profile */
    for (int i = 0; i < nsub; i++)
        off[i] = (float)(((double)SB[i] - (double)mean) * offset_divisor(i));
    if (attack != -1) {
        for (int i = 0; i < 8; i++)
            off[i] = (float)((double)off[i] + 0.7);
        for (int i = 8; i < 18; i++)
            off[i] = (float)((double)off[i] + 0.5);
        if (off[0] < 6.0f)
            off[0] = 6.0f;
    } else {
        if (off[0] < 6.0f)
            off[0] = 6.0f;
        for (int i = 1; i < 4; i++)
            if (off[i] < 3.0f)
                off[i] = 3.0f;
    }
    int quiet = word_lengths(slope, off, SB, nsub, W);
    for (int i = 0; i < nsub; i++) {
        if (SA[i] < quiet) {
            fixed[i] = 1;
            W[i] = 0;
        }
        if (multi_bits == 0 && (W[i] == 0 || (double)SA[i] < (double)mean))
            fixed[i] = 1;
    }
    memcpy(W0, W, sizeof W);                 /* upper limit for later increments */

    /* ---- 3. normalise each subband by its scale factor */
    for (int i = 0; i < nsub; i++) {
        float inv = (float)(1.0 / (double)at3_sf[SA[i]]);
        for (int k = at3_sb_start[i]; k < at3_sb_start[i + 1]; k++)
            A[k] = (float)((double)inv * (double)A[k]);
    }
    int elim = dead_zone_limit((double)budget / (double)(at3_sb_start[nsub] + 256));

    alloc_ctx c = {A, nsub, bits, 0, 0, W, E, sbits, Q, u->sel};
    for (int i = 0; i < nsub; i++)
        flags[i] = 1;
    int total = count_bits(&c, flags);

    /* ---- 4. global offset search (bisection-like, 15 iterations max) */
    float offset = 2.0f, step = 4.0f;
    memcpy(offs, off, sizeof off);
    for (int it = 0; it < 15; it++) {
        if (offset > 0) {
            for (int i = 0; i < nsub; i++)
                offs[i] = (float)((double)offset + (double)off[i]);
        } else {
            for (int i = 0; i < nsub; i++) {
                double w = offset_step_weight(i);
                offs[i] = w >= 0 ? (float)((double)off[i] + (double)offset * w)
                                 : (float)((double)offset + (double)off[i]);
            }
        }
        memcpy(W_old, W, sizeof W);
        word_lengths(slope, offs, SB, nsub, W);
        for (int i = 0; i < nsub; i++) {
            if (W0[i] == 0)
                W[i] = 0;
            flags[i] = W[i] != W_old[i];
        }
        total = count_bits(&c, flags);
        if (total > budget) {
            offset = (float)((double)offset - (double)step);
        } else {
            if (it >= 6)
                break;
            offset = total < budget ? (float)((double)step + (double)offset)
                                    : (float)((double)offset - (double)step);
        }
        step = (float)((double)step * (it < 7 ? 0.5 : 1.5));
    }

    /* ---- 5. if still too large: widen dead zones, then lower word lengths */
    for (int i = 0; i < nsub; i++)
        flags[i] = 1;
    for (int rep = 0; rep < 4; rep++)
        for (int i = nsub - 1; i >= 0 && total > budget; i--)
            if (W[i] > 0 && E[i] < elim) {
                E[i]++;
                total = count_bits(&c, flags);
            }
    if (total > budget) {
        for (int i = 0; i < nsub; i++)
            keys[i] = tonal ? (SB[i] + 1) * 32 - i : SB[i] - i / 2;
        order_descending(keys, nsub, order);
        int k = nsub;
        while (total > budget && k > 0) {
            int i = order[--k];              /* lowest-priority subband first */
            while (W[i] > 0 && total > budget) {
                W[i]--;
                for (int j = 0; j < nsub; j++)
                    flags[j] = j == i;
                total = count_bits(&c, flags);
            }
        }
    }

    /* ---- 6. fix the coded subband count; spend the remaining bits */
    int nused = nsub;
    while (W[nused - 1] == 0 && nused > 1)
        nused--;
    int *sel = u->sel, *sfi = u->sfi;
    memcpy(sel, W, nused * sizeof(int));
    memcpy(sfi, SA, nused * sizeof(int));
    c.total = total;
    c.used = bits - (budget + 3 * (nsub - nused) - total);   /* all bits of the unit so far */
    for (int i = 0; i < nused; i++)
        margin[i] = W[i] ? (float)((double)offs[i] - (double)W[i]) : (float)-9.99;
    float lo = margin[0], hi = margin[0];
    for (int i = 1; i < nused; i++) {
        if (!(lo <= margin[i]))
            lo = margin[i];
        if (margin[i] > hi)
            hi = margin[i];
    }
    double span = (double)hi - (double)lo;
    int norder = 0;                          /* subbands by margin, best first (10 classes) */
    for (int i = 0; i < nused; i++)
        rank[i] = -1;
    for (int k = 0; k < 10; k++) {
        float t = (float)((double)hi - (double)(k + 1) * span * 0.1);
        for (int i = nused - 1; i >= 0; i--)
            if (!(margin[i] < t) && rank[i] == -1) {
                rank[i] = norder;
                order[norder++] = i;
            }
    }

    /* 6a. one step per subband in margin order: raise W (up to W0) or narrow the dead zone */
    for (int oi = 0; oi < norder; oi++) {
        int i = order[oi];
        if (fixed[i] || !(0 < sel[i] && sel[i] < 7))
            continue;
        int s0 = sel[i], e0 = E[i];
        if (s0 < W0[i]) {
            sel[i]++;
        } else if (E[i] > 0) {
            E[i]--;
        } else {
            if (c.used > bits - 8)
                break;
            continue;
        }
        int b_new = quant_sb(A, i, sel[i], E[i], q);
        if (c.used + b_new - sbits[i] > bits) {
            sel[i] = s0;
            E[i] = e0;
            continue;
        }
        commit(&c, i, b_new, q);
        if (c.used > bits - 8)
            break;
    }
    /* 6b. remove remaining dead zones where affordable */
    for (int i = 0; i < nused; i++)
        if (E[i] > 0 && sel[i] > 0) {
            int b_new = quant_sb(A, i, sel[i], 0, q);
            if (c.used + b_new - sbits[i] <= bits) {
                commit(&c, i, b_new, q);
                E[i] = 0;
            }
        }
    /* 6c. raise word lengths towards W0 (7 passes) */
    for (int rep = 0; rep < 7 && !(c.used > bits - 2); rep++)
        for (int i = 0; i < nused && !(c.used > bits - 2); i++)
            if (0 < sel[i] && sel[i] < W0[i] && !fixed[i])
                if (!try_raise(&c, i))
                    fixed[i] = 1;
    /* 6d. raise low word lengths in scale-factor order (4 passes, limit 3..6) */
    if (nused > 1) {
        for (int i = 0; i < nused; i++)
            keys[i] = (SB[i] + 1) * 32 - i;
        order_descending(keys, nused, order);
        for (int p = 0; p < 4 && !(c.used > bits - 2); p++)
            for (int k = 0; k < nused && !(c.used > bits - 2); k++) {
                int i = order[k];
                if (0 < sel[i] && sel[i] < p + 3)
                    try_raise(&c, i);
            }
    }
    /* 6e. raise any word length below 7 (7 passes) */
    for (int rep = 0; rep < 7 && !(c.used > bits - 2); rep++)
        for (int i = 0; i < nused && !(c.used > bits - 2); i++)
            if (0 < sel[i] && sel[i] < 7)
                try_raise(&c, i);

    /* ---- 7. scale-factor correction: if the quantised energy exceeds 1.25x the original,
     *         use the next smaller scale factor */
    for (int i = 0; i < nused; i++) {
        float st = (float)(1.0 / ((double)MAXQ[sel[i]] + 0.5));
        double e_orig = 0.0, e_quant = 0.0;
        for (int k = at3_sb_start[i]; k < at3_sb_start[i + 1]; k++) {
            double x = A[k];
            e_orig = e_orig + x * x;
            double y = (double)Q[k] * (double)st;
            e_quant = e_quant + y * y;
        }
        if (e_orig * 1.25 < e_quant && sfi[i] > 0)
            sfi[i]--;
    }

    u->nbands = nbands;
    u->nsub = nused;
    u->mode = 0;
    u->tone_mode = TONE_MODE;
    u->bits = c.used;
}
