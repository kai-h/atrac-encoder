/* Bit allocation, tonal component extraction and quantisation of the LP4 encoder
 * (allocate and helpers).  Mirrors spec/reference/lp4_alloc.py.  Integer logic is
 * literal; the tone quantiser is x87 code (doubles, float stores) and the spectral
 * quantiser/bit counter is SSE code (float multiply, then a float add of 1.5 * 2**23 whose
 * bit pattern holds the rounded value). */
#include <math.h>
#include <stdlib.h>
#include "internal.h"

#define MAGIC 12582912.0                     /* 1.5 * 2**23 */
#define POS(rec) ((rec)->r[4])
#define SF(rec) ((rec)->r[5])

/* ---------------------------------------------------------------- basic helpers */
static int low16(float x)                    /* signed low 16 bits of the float pattern */
{
    int v = (int)(f32_bits(x) & 0xffff);
    return v >= 0x8000 ? v - 0x10000 : v;
}

/* sf_index4: smallest index whose 2**((i-15)/3) is not below the peak of 4 values,
 * computed from the bit pattern of the peak. */
static int sf_index4(const float *v)
{
    uint32_t u = 0;
    for (int k = 0; k < 4; k++) {
        uint32_t b = f32_bits(v[k]) << 1;
        if (b > u)
            u = b;
    }
    int e = (int)(u >> 24), m = (int)(u & 0xffffff);
    int i = 3 * e - 0x16c;
    if (m > 0x965fe9)
        i++;
    if (m < 0x428a30)
        i--;
    return (i >= 0 && i <= 63) ? i : 0;
}

/* qscale: MAX_QUANT[sel] / 2**((sf-15)/3) as float32, built like the binary does. */
static float qscale(int sf, int sel)
{
    uint32_t u = (uint32_t)sf * 0x2b0000u;
    return f32_from_bits(lp4_qscale_bits[((u >> 23) + sel) * 3 - sf] - (u & 0x7f800000u));
}

/* quantise: rounded integer (as bit pattern - 0x4b400000) of x * scale. */
uint32_t lp4_quantise(float x, int sel, int sf)
{
    float v = x * qscale(sf, sel);
    v = v + (float)MAGIC;
    return f32_bits(v) - 0x4b400000u;
}

/* bitcount: bits of one subband, 6 (scale factor) + VLC lengths.  Values outside the
 * selector's range wrap through the index mask exactly as in Sony's code. */
static int bitcount(int sel, int sf, int n, const float *x)
{
    float s = qscale(sf, sel);
    int tot = 6;
    if (sel == 1) {
        for (int g = 0; g < n; g += 4) {
            int idx = 0;
            for (int k = 0; k < 4; k++) {
                float v = x[g + k] * s;
                v = v + (float)MAGIC;
                if (v != (float)MAGIC)
                    idx |= 1 << k;
            }
            tot += lp4_sel1_len[idx];
        }
    } else {
        for (int k = 0; k < n; k++) {
            float v = x[k] * s;
            v = v + (float)MAGIC;
            tot += lp4_vlc[0][sel][(f32_bits(v) - 0x4b400000u) & lp4_qmask[sel]].len;
        }
    }
    return tot;
}

static int clamp_sel(int v, int shift)
{
    v = asr(v, shift);
    return v < 1 ? 1 : (v > 7 ? 7 : v);
}

void lp4_alloc_init(lp4_alloc_state *S, const float *spec)
{
    memset(S, 0, sizeof *S);
    memcpy(S->spec, spec, sizeof S->spec);
}

/* ---------------------------------------------------------------- tones */
/* best_sf_tone: try the 3 scale factors above the peak's index and keep the one with the
 * smallest maximum rounding error, then move up in octaves while all values are even. */
static void best_sf_tone(const float *spec, int pos, lp4_comp *rec, int qstep, int n)
{
    double best = lp4_tone_best_init;
    double x[4];
    for (int k = 0; k < 4; k++)
        x[k] = spec[pos + k];
    int lo = sf_index4(spec + pos);
    int hi = imin(lo + 3, 64) - 1;
    for (int sf = hi; sf >= lo; sf--) {
        double scale = qscale(sf, qstep);
        double err = 0.0;
        int acc = 0;
        for (int k = 0; k < n; k++) {
            double p = scale * x[k];
            int q = low16((float)(p + MAGIC));
            acc |= q;
            double d = fabs(p - q);
            if (!(err >= d))
                err = (float)d;
        }
        if (!(scale * best <= err)) {
            best = err / scale;
            SF(rec) = sf;
            if ((acc & 1) == 0)
                while (SF(rec) < 61) {
                    acc = asr(acc, 1);
                    SF(rec) += 3;
                    if (acc & 1)
                        break;
                }
        }
    }
}

/* quant_tone: quantise a tonal component at its position and subtract it from the
 * spectrum.  Returns its VLC cost in bits (6 sf + 6 position + codes). */
static int quant_tone(float *spec, lp4_comp *rec, int qstep, int n)
{
    best_sf_tone(spec, POS(rec), rec, qstep, n);
    double scale = qscale(SF(rec), qstep);
    double inv = (float)(1.0 / scale);
    int nb = 12;
    for (int k = 0; k < n; k++) {
        int p = POS(rec) + k;
        double x = spec[p];
        int q = low16((float)(scale * x + MAGIC));
        int idx = q & lp4_qmask[qstep];
        spec[p] = (float)(x - q * inv);
        rec->r[k] = idx;
        nb += lp4_vlc[0][qstep][idx].len;
    }
    return nb;
}

/* restore_tone: add a component back into the spectrum (uses the unrounded reciprocal). */
static void restore_tone(float *spec, const lp4_comp *rec, int nvm1, int qstep)
{
    double inv = 1.0 / qscale(SF(rec), qstep);
    int nbits = lp4_vlc[1][qstep][0].len;
    for (int k = nvm1; k >= 0; k--) {
        int v = rec->r[k];
        if (v >= 1 << (nbits - 1))
            v -= 1 << nbits;
        int p = POS(rec) + k;
        spec[p] = (float)(v * inv + (double)spec[p]);
    }
}

/* find_peaks: positions with |x| >= SF value of sfthr, skipping 3 after each hit.
 * Returns the count, or -1 if more than maxn peaks exist. */
static int find_peaks(const float *spec, int sfthr, int start, int end, int maxn, int *out)
{
    uint32_t thr2 = lp4_sf_bits[sfthr] * 2u;
    int n = 0, k = start;
    for (;;) {
        if ((f32_bits(spec[k]) << 1) >= thr2) {
            if (n == maxn)
                return -1;
            out[n++] = k;
            k += 3;
        }
        k++;
        if (k >= end)
            return n;
    }
}

/* shift_comp: drop leading zero values of a component by moving it right. */
static int shift_comp(lp4_comp *rec, int i)
{
    int e = 1;
    while (e < 4 && rec->r[e] == 0)
        e++;
    int np_ = POS(rec) + e;
    if (np_ >= at3_sb_start[i + 1])
        return -1;
    if (np_ > 0x3fc)
        return 0;
    int ncopy = imax(4 - e, 1);
    for (int c = 0; c < ncopy; c++)
        rec->r[c] = rec->r[c + e];
    for (int c = ncopy; c < 4; c++)
        rec->r[c] = 0;
    POS(rec) = np_;
    return e;
}

/* Python list indexing of a record (q[-1] is the sf field, and so on). */
#define QI(rec, k) ((rec)->r[((k) + 6) % 6])

/* shrink_values: reduce the number of values per component while the trailing value is zero
 * for all components.  Returns how many values were removed. */
static int shrink_values(lp4_alloc_state *S, int n)
{
    lp4_group *g = &S->groups[0];
    for (;;) {
        int flag = 0, nv = g->nvm1;
        if (nv < -5)
            return 3 - nv;
        for (int j = n - 1; j >= 0; j--) {
            lp4_comp *rec = &S->comps[j];
            if (QI(rec, nv) != 0) {
                POS(rec)++;
                for (int c = 0; c < nv; c++)
                    rec->r[c] = rec->r[c + 1];
            }
            if (QI(rec, nv - 1) != 0 && (rec->r[0] != 0 || (POS(rec) & 0x3f) == 0x3f))
                flag = 1;
        }
        g->nvm1 = nv - 1;
        if (flag)
            return 3 - g->nvm1;
    }
}

/* split_comps: split every 4-value component whose last value is nonzero into two 2-value
 * components.  Returns the new component count, or 0 if not possible. */
static int split_comps(lp4_alloc_state *S, int n)
{
    lp4_group *g = &S->groups[0];
    int cnt[16];
    for (int b = 0; b < 16; b++)
        cnt[b] = g->nlst[b];
    int total = n;
    for (int j = 0; j < n; j++) {
        lp4_comp *rec = &S->comps[j];
        if (rec->r[3] != 0) {
            int p = POS(rec) + 2, blk = p >> 6;
            if (g->flag[blk >> 2] == 0 || cnt[blk] >= 7 || total == 64 || p > 0x3fd)
                return 0;
            cnt[blk]++;
            total++;
        }
    }
    total = n;
    for (int j = 0; j < n; j++) {
        lp4_comp *rec = &S->comps[j];
        if (rec->r[3] != 0) {
            lp4_comp *nw = &S->comps[total];
            nw->r[0] = rec->r[2];
            nw->r[1] = rec->r[3];
            nw->r[2] = nw->r[3] = 0;
            nw->r[4] = POS(rec) + 2;
            nw->r[5] = SF(rec);
            rec->r[2] = rec->r[3] = 0;
            int blk = POS(nw) >> 6;
            g->lst[blk][g->nlst[blk]++] = total;
            total++;
        }
    }
    return total;
}

/* extract_tones: tonal component extraction for the normal allocation path.  Peaks well
 * above their subband's level are coded as 4-value components (quantiser step 3) and removed
 * from the spectrum.  Marks affected subbands with a negative sf and reuses sel[] as the
 * number of coefficients covered by tones.  Returns the bits spent on the tonal block. */
static int extract_tones(int bits_left, const int *P, lp4_alloc_state *S)
{
    lp4_group *g = &S->groups[0];
    memset(g, 0, sizeof *g);
    g->qstep = 3;
    g->nvm1 = 3;
    int hdr = 8, ncomp = 0, bits_tone = 0, nzero = 0;
    for (int i = bits_left < 1100 ? 8 : 0; i < S->nsub; i++) {
        if (ncomp > 64 - (lp4_sb_width[i] >> 4) || 24 * ncomp + hdr > bits_left - 600)
            break;
        int sel = S->sel[i];
        if (sel == 0)
            continue;
        if (P[i] < 0x500)
            sel = 0;
        int start = at3_sb_start[i];
        int band = lp4_sb_width[i] == 32 ? 1 : start >> 8;
        int thr_sf = imax(S->sf[i] - lp4_tone_sf_offset[sel + band * 8], 0);
        int st = start;
        if (ncomp > 0 && POS(&S->comps[ncomp - 1]) + 4 > start) {
            S->sf[i] = -S->sf[i];
            int ov = POS(&S->comps[ncomp - 1]) - start;
            S->sel[i] = ov + 4;
            S->sel[i - 1] += ov;
            st = POS(&S->comps[ncomp - 1]) + 4;
        }
        int peaks[8];
        int np = find_peaks(S->spec, thr_sf, st, at3_sb_start[i + 1], (lp4_sb_width[i] + 8) >> 4, peaks);
        if (np <= 0)
            continue;
        if (S->sf[i] >= 0) {
            S->sf[i] = -S->sf[i];
            S->sel[i] = 0;
        }
        for (int pi = 0; pi < np; pi++) {
            int c = imin(peaks[pi], 0x3fc), nv = 4;
            if (ncomp > 0 && POS(&S->comps[ncomp - 1]) + 4 == c)
                for (;;) {
                    if (fabs((double)S->spec[c + 3]) > fabs((double)S->spec[c - 1]))
                        break;
                    c--;
                    nv--;
                    if (nv <= 1)
                        break;
                }
            int blk = c >> 6;
            if (g->nlst[blk] >= 7)
                continue;
            int ri = ncomp++;
            lp4_comp *rec = &S->comps[ri];
            POS(rec) = c;
            int nb = quant_tone(S->spec, rec, 3, 4);
            if (rec->r[0] == 0 && POS(rec) <= 0x3fc) {
                if (shift_comp(rec, i) < 0) {
                    restore_tone(S->spec, rec, 3, 3);
                    ncomp--;
                    break;
                }
                blk = POS(rec) >> 6;
            }
            bits_tone += nb;
            S->sel[i] += nv;
            if (g->flag[blk >> 2] == 0) {
                g->flag[blk >> 2] = 1;
                hdr += 12;
            }
            g->lst[blk][g->nlst[blk]++] = ri;
            if (rec->r[3] == 0)
                nzero++;
        }
    }
    if (ncomp == 0)
        return 0;
    if (ncomp == 1 && (POS(&S->comps[0]) < 128 || bits_left < 1100)) {
        restore_tone(S->spec, &S->comps[0], 3, 3);
        return 0;
    }
    hdr += S->nbands;
    S->ngroups = 1;
    int len0v = lp4_vlc[0][g->qstep][0].len, len0c = lp4_vlc[1][g->qstep][0].len;
    int vlc = bits_tone, clc = 24 * ncomp;
    int do_shrink = nzero == ncomp;
    if (!do_shrink) {
        int nz = 2 * ncomp - 3 * nzero, x12 = 12 * (ncomp - nzero);
        int est_v = len0v * nz + x12 + bits_tone;
        int est_c = len0c * nz + x12 + 24 * ncomp;
        if (imin(bits_tone, 24 * ncomp) > imin(est_v, est_c)) {
            int r = split_comps(S, ncomp);
            if (r) {
                ncomp = r;
                vlc = len0v * r + est_v;
                clc = len0c * r + est_c;
                do_shrink = 1;
            }
        }
    }
    if (do_shrink) {
        int k = shrink_values(S, ncomp) * ncomp;
        vlc -= len0v * k;
        clc -= len0c * k;
    }
    int best;
    if (vlc > clc) {
        S->sel_mode = 1;
        best = clc;
    } else {
        S->sel_mode = 0;
        best = vlc;
    }
    return hdr + best;
}

/* replace_weakest: make room in a full block list: remove (restore) the component with the
 * lowest sf below `sf` and return its record for reuse, or -1. */
static int replace_weakest(lp4_alloc_state *S, int gi, int blk, int sf)
{
    lp4_group *grp = &S->groups[gi];
    int best = -1;
    for (int j = 0; j < 7; j++)
        if (SF(&S->comps[grp->lst[blk][j]]) < sf) {
            best = j;
            sf = SF(&S->comps[grp->lst[blk][j]]);
        }
    if (best < 0)
        return -1;
    int ri = grp->lst[blk][best];
    for (int j = best; j < grp->nlst[blk] - 1; j++)
        grp->lst[blk][j] = grp->lst[blk][j + 1];
    grp->nlst[blk]--;
    restore_tone(S->spec, &S->comps[ri], grp->nvm1, grp->qstep);
    return ri;
}

/* all_tonal: allocation for frames with plenty of bits.  Every coefficient of magnitude >= 2
 * becomes a tonal component (CLC coded, step 5 for sf <= 30, step 7 above) until the
 * spectrum is exhausted or the bits run out; the remaining subbands get selector 1. */
static int all_tonal(int bits_left, int *P, int *sfidx, lp4_alloc_state *S)
{
    int need = at3_sb_start[S->nsub] >> 8;
    int cost = 22;
    if (S->nbands < need) {
        cost = (need - S->nbands) * 3 + 22;
        S->nbands = need;
    }
    memset(S->groups, 0, sizeof S->groups);
    S->groups[0].qstep = 5;
    S->groups[1].qstep = 7;
    S->groups[0].nvm1 = S->groups[1].nvm1 = 3;
    S->sel_mode = 1;
    S->ngroups = 2;
    int ncomp = 0, limit = bits_left - 200;
    uint32_t thr2 = lp4_sf_bits[18] * 2u;    /* 2.0 */
    int end = at3_sb_start[S->nsub], done = 0;
    while (!done && end > 0) {
        int start_n = ncomp, k = 0;
        while (k < end) {
            if ((f32_bits(S->spec[k]) << 1) >= thr2) {
                int ri = ncomp++;
                POS(&S->comps[ri]) = imin(k, 0x3fc);
                int sf = sf_index4(S->spec + POS(&S->comps[ri]));
                int gi = sf > 30 ? 1 : 0;
                int blk = k >> 6;
                if (S->groups[gi].nlst[blk] == 7) {
                    if (S->groups[1 - gi].nlst[blk] == 7) {
                        ncomp--;
                        ri = replace_weakest(S, gi, blk, sf);
                        if (ri < 0) {
                            k++;
                            continue;
                        }
                    } else {
                        gi = 1 - gi;
                    }
                }
                lp4_group *grp = &S->groups[gi];
                if (grp->flag[blk >> 2] == 0) {
                    grp->flag[blk >> 2] = 1;
                    cost += 12;
                }
                grp->lst[blk][grp->nlst[blk]++] = ri;
                quant_tone(S->spec, &S->comps[ri], grp->qstep, 4);
                cost += 28 + 8 * gi;
                if (cost > limit || ncomp >= 64) {
                    done = 1;
                    break;
                }
                k += 3;
            }
            k++;
        }
        if (start_n == ncomp)
            break;
    }
    for (int i = 0; i < S->nsub; i++) {
        int mx = 0;
        for (int gg = at3_sb_start[i] >> 2; gg < at3_sb_start[i + 1] >> 2; gg++) {
            int v = sf_index4(S->spec + 4 * gg);
            sfidx[gg] = v;
            mx = imax(mx, v);
        }
        S->sf[i] = mx;
        S->sel[i] = mx ? 1 : 0;
        P[i] = mx * 341;
    }
    return cost;
}

/* ---------------------------------------------------------------- allocation */
static int est_bits(int i, const int *sfidx, const lp4_alloc_state *S)
{
    int sel = S->sel[i];
    int thr = S->sf[i] - lp4_est_thr[sel];
    int b = lp4_est_bits[sel] * lp4_sb_width[i] + 60;
    int g0 = at3_sb_start[i] >> 2;
    for (int g = g0; g < g0 + (lp4_sb_width[i] >> 2); g++)
        if (sfidx[g] < thr)
            b -= lp4_est_save[sel];
    return b;
}

static int gain_header_bits(const lp4_unit_info *U, lp4_alloc_state *S)
{
    int nb = (at3_sb_start[S->nsub] + 0xff) >> 8, pts = 0;
    S->nbands = nb;
    for (int b = 0; b < nb; b++)
        pts += U->gain_counts[b];
    return (pts * 3 + nb) * 3;
}

static int drop_empty_bands(const lp4_unit_info *U, lp4_alloc_state *S)
{
    if (U->js && (U->weight_sent & 7) != 7)
        return 0;
    int n = S->nbands, k = n - 1;
    while (k > 0 && U->gain_counts[k] == 0) {
        int any = 0;
        for (int j = 0; j < S->ngroups; j++)
            any |= S->groups[j].flag[k];
        if (any)
            break;
        k--;
    }
    S->nbands = k + 1;
    return (S->ngroups + 3) * (n - (k + 1));
}

/* allocate: choose nsub, selectors, scale factors and tonal components for one unit. */
void lp4_allocate(lp4_alloc_state *S, const lp4_unit_info *U, const float *spec)
{
    int sfidx[256], sums[32], P[32], order[32];
    lp4_alloc_init(S, spec);
    int budget = U->budget, nsub = U->nsub_max;
    for (int g = 0; g < 256; g++)
        sfidx[g] = sf_index4(S->spec + 4 * g);
    int cnt7 = 0, cum = 0, top = 0;
    for (int i = 0; i < 32; i++) {
        int s = 0, mx = 0;
        for (int g = at3_sb_start[i] >> 2; g < at3_sb_start[i + 1] >> 2; g++) {
            int v = sfidx[g];
            s += v;
            mx = imax(mx, v);
            if (v > 7)
                cnt7++;
        }
        sums[i] = s;
        if (top < mx) {
            top = mx;
            if (!U->js) {                    /* a new loudest subband may extend the bandwidth */
                int ok = ((at3_sb_start[i] & 0x80) == 0 || i >= 22) ? cnt7 * 8 < budget
                                                                       : cnt7 * 16 < budget;
                if (ok)
                    nsub = imax(nsub, i <= 26 ? 28 : imin(i + 2, 32));
            }
        } else if (s < lp4_min_sfsum[i] || mx < 3) {
            s = 0;
            mx = 0;
        }
        S->sf[i] = mx;
        P[i] = cnt7;
        cum += s;
        S->sel[i] = cum;
    }
    S->nsub = nsub;
    int shift = budget * 2 < P[nsub - 1] * 12 ? 11 : 10;
    int bits_left = budget - 3 * nsub - gain_header_bits(U, S);
    if (bits_left > P[nsub - 1] * 32) {
        bits_left -= all_tonal(bits_left, P, sfidx, S);
    } else {
        for (int i = 0; i < nsub; i++) {
            int m = S->sf[i], v = 0;
            if (m >= 3) {
                m = m * 256 - cum;
                v = clamp_sel(m, shift);
            }
            S->sel[i] = v;
            P[i] = m;
        }
        S->ngroups = 0;
        if (!U->js)
            bits_left -= extract_tones(bits_left, P, S);
        for (int i = 0; i < nsub; i++) {
            int m = S->sf[i];
            if (m < 0) {                     /* subband lost tonal components: re-evaluate */
                int mx = 0, s = 0;
                for (int g = at3_sb_start[i] >> 2; g < at3_sb_start[i + 1] >> 2; g++) {
                    int v = sf_index4(S->spec + 4 * g);
                    sfidx[g] = v;
                    mx = imax(mx, v);
                    s += v;
                }
                S->sf[i] = mx;
                if (S->sel[i] == 0) {
                    P[i] += 256;
                } else {
                    int k = asr(S->sel[i] - 1, 1);
                    if (k < 0 || k > 4)
                        k = 4;
                    int d = lp4_tone_penalty[k] * (mx + m);
                    if (budget < 1100) {
                        d *= 2;
                    } else if (lp4_sb_width[i] <= 16 && P[i] > 0x1800) {
                        int d2 = 2 * d;
                        d = d2 < P[i] - 0x400 ? d2 : P[i] - 0x400;
                    }
                    if (i > 1)
                        P[i - 1] -= asr(d, 3);
                    if (i < 31)
                        P[i + 1] -= asr(d, 2);
                    P[i] -= d;
                }
                int v = 0;
                if (mx > 7 && s > lp4_min_sfsum[i] && (P[i] >= 0x400 || budget > 0x708))
                    v = clamp_sel(P[i], shift);
                S->sel[i] = v;
            } else if (shift != 10) {
                if (lp4_sb_width[i] * m < sums[i] * 5) {
                    if (m >= 62)
                        continue;
                    S->sf[i] = m + 1;
                }
                if (S->sf[i] < 62)
                    S->sf[i] += 1;
            } else if (lp4_sb_width[i] * m > sums[i] * 8) {
                P[i] += 256;
            }
        }
    }
    /* drop isolated subbands far below their neighbours, estimate the bits of the rest */
    int est = 0, w7 = 0, wsum = 0;
    for (int i = 0; i < nsub; i++) {
        if (S->sel[i] == 0)
            continue;
        if ((i > 1 && P[i] < P[i - 1] - 0xf00) || (i < nsub - 1 && P[i] < P[i + 1] - 0x1400)) {
            S->sel[i] = 0;
        } else {
            wsum += lp4_sb_width[i];
            if (S->sel[i] == 7)
                w7 += lp4_sb_width[i];
            est += est_bits(i, sfidx, S);
        }
    }
    int k = nsub - 1;
    while (k > 0 && S->sel[k] == 0) {
        bits_left += 3;
        k--;
    }
    nsub = k + 1;
    S->nsub = nsub;
    bits_left += drop_empty_bands(U, S);
    int b10 = bits_left * 10, off = 0, width = 0;
    if (b10 <= est)
        width = wsum;
    else if (wsum - w7 != 0)
        width = wsum - w7;
    if (width != 0) {
        int t = bits_left < 600 ? b10 : bits_left * 9;
        int diff = t - est, num = diff * 1024;
        off = abs(num) / (width * 10);
        if ((num < 0) != (width < 0))
            off = -off;
        if (diff < 0 && shift == 10) {
            off += asr(off, 3);
            P[0] -= off;
            P[1] -= asr(off * 188, 8);
            for (k = 2; k < 8; k++)
                P[k] -= asr(off * 97, 8);
            for (k = 8; k < 18; k++)
                P[k] -= asr(off * 74, 8);
            for (k = 0; k < 18; k++)
                if (P[k] + off < 0 && S->sf[k] < 63)
                    S->sf[k]++;
            for (k = 18; k < nsub; k++)
                if (P[k] + off < 0x400)
                    S->sf[k] = imin(63, S->sf[k] - asr(P[k] + off - 0x400, 10));
        }
    }
    for (int i = nsub - 1; i >= 0; i--)
        if (S->sel[i] != 0) {
            S->sel[i] = clamp_sel(P[i] + off, shift);
            b10 -= est_bits(i, sfidx, S);
            P[i] = b10;
        }
    /* exact bit counts, local selector adjustment */
    int used = 0;
    for (int i = 0; i < 32; i++)
        order[i] = -1;
    for (int i = 0; i < nsub; i++) {
        int sel = S->sel[i];
        if (sel == 0)
            continue;
        if (lp4_sb_width[i] > (bits_left - used) * 2 - 16) {
            S->sel[i] = 0;
            continue;
        }
        int d = P[i] - used * 10;
        if (sel < 7 && d > 4000) {
            S->sel[i] = sel + 1;
        } else if (i != 0 && d < 300) {
            if (sel > 1)
                S->sel[i] = sel - 1;
            else if (S->sf[i] != 63)
                S->sf[i]++;
        }
        int b = bitcount(S->sel[i], S->sf[i], lp4_sb_width[i], S->spec + at3_sb_start[i]);
        P[i] = b;
        if (b * 2 == lp4_zero_len2[S->sel[i]] * lp4_sb_width[i] + 12) {   /* everything quantised to 0 */
            S->sel[i] = 0;
        } else {
            used += b;
            order[i] = S->sf[i];
        }
    }
    k = nsub - 1;
    while (k > 0 && S->sel[k] == 0) {
        bits_left += 3;
        k--;
    }
    int slack = bits_left - used;
    nsub = k + 1;
    S->nsub = nsub;
    int spent = 0;
    if (U->prev_count0 + U->gain_counts[0] != 0 && !U->js)
        for (k = 0; k < 5; k++)              /* transient frames: favour the lowest subbands */
            if (order[k] >= 0 && order[k] <= 0x40)
                order[k] = 0x40;
    /* final pass in order of decreasing scale factor: fit the budget, spend leftovers */
    for (;;) {
        int best = -1, idx = 0;
        for (k = 0; k < nsub; k++)
            if (best < order[k]) {
                best = order[k];
                idx = k;
            }
        if (best < 0)
            return;
        int b = P[idx];
        order[idx] = -1;
        int w = lp4_sb_width[idx];
        const float *x = S->spec + at3_sb_start[idx];
        if (spent + b > bits_left) {
            int zthr = lp4_zero_len2[S->sel[idx]] * (w >> 1) + 6;
            k = nsub - 1;
            while (k > idx) {
                if (S->sel[k] != 0 && order[k] == -1)
                    break;
                S->sel[k] = 0;
                bits_left += 3;
                order[k] = -1;
                k--;
            }
            nsub = k + 1;
            S->nsub = nsub;
            while (spent + b > bits_left) {
                S->sf[idx]++;
                if (S->sf[idx] >= 64 || b <= zthr)
                    break;
                b = bitcount(S->sel[idx], S->sf[idx], w, x);
            }
        } else if (S->sel[idx] < 7 && lp4_step_cost[S->sel[idx]] * (w >> 1) < bits_left - spent - 7) {
            int ns = S->sel[idx] + 1;
            if (slack > w * 4 && ns < 7)
                ns++;
            int b2 = bitcount(ns, S->sf[idx], w, x);
            if (b2 - b <= (slack > 0 ? slack : 0)) {
                S->sel[idx] = ns;
                slack += b - b2;
                b = b2;
            }
        }
        spent += b;
        if (spent > bits_left) {
            spent -= b;
            slack += b;
            S->sel[idx] = 0;
        }
    }
}
