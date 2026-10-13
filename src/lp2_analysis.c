/* LP2 encoder analysis front end: QMF band split, gain control, MDCT.
 * Mirrors spec/reference/lp2_analysis.py.  Sony's code is x87 with 53-bit precision: values
 * kept in FPU registers are doubles here, values stored to memory are floats. */
#include <math.h>
#include "internal.h"

/* ================================================================ QMF analysis (qmf_split) */
static void qmf_split(const float *x, int n, float *delay, float *lo, float *hi)
{
    double s[46 + 1024];
    for (int i = 0; i < 46; i++)
        s[i] = delay[i];
    for (int i = 0; i < n; i++)
        s[46 + i] = x[i];
    for (int i = 0; i < 46; i++)
        delay[i] = x[n - 46 + i];
    for (int j = 0; j < n / 2; j++) {
        double P = 0.0, Q = 0.0;
        for (int k = 47; k > 0; k -= 2)
            P = P + (double)lp2_qmf_win[k] * s[2 * j + k];
        for (int k = 46; k >= 0; k -= 2)
            Q = (double)lp2_qmf_win[k] * s[2 * j + k] + Q;
        lo[j] = (float)(Q + P);
        hi[j] = (float)(P - Q);
    }
}

/* ================================================================ gain control */
static int log2_round(double r)
{
    if (r <= 1e-20)
        r = 1e-20;
    double lg2 = log10(2.0);
    double a = (log(r) / log(2.0)) * lg2;
    return (int)(a / lg2 + 0.5);
}

static float absmax(const float *x, int n)
{
    float m = 0.0f;
    for (int i = 0; i < n; i++) {
        float v = fabsf(x[i]);
        if (v > m)
            m = v;
    }
    return m;
}

/* band_gain_points: gain points of one band from its 768-sample history [prev | cur | next]. */
static void band_gain_points(const float *x, const lp2_gain *prev, lp2_gain *out)
{
    float blk[64], grp[16];
    for (int i = 0; i < 64; i++)
        blk[i] = absmax(x + 256 + 8 * i, 8);
    float pre = absmax(x + 224, 32);
    int prev_lev = 0;
    for (int k = 0; k < prev->n; k++)
        prev_lev = imax(prev_lev, prev->lev[k] - 4);
    for (int g = 0; g < 16; g++) {
        grp[g] = blk[4 * g];
        for (int k = 1; k < 4; k++)
            if (blk[4 * g + k] > grp[g])
                grp[g] = blk[4 * g + k];
    }
    const double k185 = (double)1.85f;

    /* candidate list: attack ratios R[0..31] and release ratios R[32..39] */
    float R[40] = {0};
    int cnt = 0;
    double M = prev->peak;
    for (int i = 0; i < 32; i++) {
        if (M < blk[i])
            M = blk[i];
        double v = blk[i + 1];
        if (v > 10.0 && v > 1.5 * M) {
            R[i] = M <= 4.0 ? (float)(v * 0.25) : (float)(v / M);
            cnt++;
        }
    }
    double M2 = blk[63];
    for (int g = 15; g > 7; g--)
        if (M2 < grp[g])
            M2 = grp[g];
    for (int g = 7; g >= 0; g--) {
        if (M2 < grp[g])
            M2 = grp[g];
        double v = g == 0 ? pre : grp[g - 1];
        if (v > 10.0 && v > M2 * k185) {
            R[32 + g] = M2 <= 4.0 ? (float)(0.25 * v) : (float)(v / M2);
            cnt++;
        }
    }
    int sel[40] = {0};                       /* keep the 7 strongest candidates */
    for (int t = 0; t < imin(cnt, 7); t++) {
        double best = 0.0;
        int idx = 0;
        for (int i = 0; i < 40; i++)
            if (best < R[i]) {
                best = R[i];
                idx = i;
            }
        R[idx] = 0.0f;
        sel[idx] = 1;
    }

    /* attacks: level steps, at most 7, total level <= 10 */
    int att_pos[8], att_val[8], natt = 0, tot = prev_lev;
    M = prev->peak;
    for (int i = 0; i < 32; i++) {
        if (M < blk[i])
            M = blk[i];
        double v = blk[i + 1];
        if (!(v > 10.0 && v > 1.5 * M && sel[i]))
            continue;
        int a = log2_round(M <= 4.0 ? v * 0.25 : v / M);
        if (a >= 0) {
            if (a + tot > 10)
                a = 10 - tot;
            tot += a;
            if (a > 0) {
                att_pos[natt] = i;
                att_val[natt++] = a;
            }
        }
        if (natt >= 7)
            break;
    }

    /* releases, scanned backwards in 32-sample groups, total <= min(tot, 4) */
    int cap = imin(tot, 4);
    int rel_pos[8], rel_val[8], nrel = 0;
    if (natt || prev->n) {
        int rs = 0, n_all = natt;
        M2 = blk[63];
        for (int g = 15; g > 7; g--)
            if (M2 < grp[g])
                M2 = grp[g];
        for (int g = 7; g >= 0; g--) {
            if (M2 < grp[g])
                M2 = grp[g];
            double v = g == 0 ? pre : grp[g - 1];
            if (v > 10.0 && v > M2 * k185 && sel[32 + g]) {
                int a = log2_round(M2 <= 4.0 ? v * 0.25 : v / M2);
                if (a >= 0) {
                    if (a + rs > cap)
                        a = cap - rs;
                    rs += a;
                    if (a > 0 && natt < 7) {
                        rel_pos[nrel] = g == 0 ? 1 : 4 * g;
                        rel_val[nrel++] = a;
                        n_all++;
                    }
                }
            }
            if (n_all >= 7)
                break;
        }
    }

    out->peak = absmax(blk, 32);
    /* cumulative levels (from the end), releases capped at 4, attacks at 11 + total release */
    int acc = 0;
    for (int k = nrel - 1; k >= 0; k--) {
        rel_val[k] = imin(rel_val[k] + acc, 4);
        acc = rel_val[k];
    }
    int cap_a = 11 + acc;
    acc = 0;
    for (int k = natt - 1; k >= 0; k--) {
        att_val[k] = imin(att_val[k] + acc, cap_a);
        acc = att_val[k];
    }
    int lv[33] = {0};
    int j = 0;
    for (int k = 0; k < natt; k++)
        while (j <= att_pos[k])
            lv[j++] += att_val[k];
    j = 32;
    for (int k = 0; k < nrel; k++)
        while (j >= rel_pos[k])
            lv[j--] += rel_val[k];
    out->n = 0;
    for (j = 0; j < 32; j++)
        if (lv[j] != lv[j + 1]) {
            out->lev[out->n] = lv[j] - lv[32] + 4;
            out->loc[out->n++] = j;
        }
    while (out->n && out->lev[out->n - 1] == 4)
        out->n--;
}

/* adjust_band0_gain: give band 0 an attack point when band 1 has a strong attack (level spread
 * >= 2) and band 0 has none and is not loud. */
static void adjust_band0_gain(float hist[4][768], const lp2_gain *prev, lp2_gain *cur)
{
    int spread = 0;
    const lp2_gain *p = &cur[1];
    if (p->n) {
        int mn = 4, mi = 0, mx = 0;
        for (int i = 0; i < p->n; i++)
            if (p->lev[i] < mn) {
                mn = p->lev[i];
                mi = i;
            }
        for (int i = 0; i <= mi; i++)
            mx = imax(mx, p->lev[i]);
        spread = mx - mn;
    }
    if (spread < 2 || prev[0].n || cur[0].n)
        return;
    int loc = cur[1].loc[0];
    int ok = !(prev[0].peak > 16384.0);
    for (int k = 256; k < 256 + 8 * (loc + 1); k++)
        if (fabs((double)hist[0][k]) > 16384.0)
            ok = 0;
    if (ok) {
        cur[0].n = 1;
        cur[0].lev[0] = 5;
        cur[0].loc[0] = loc;
    }
}

/* gain_window: 512-sample divisor applied before the MDCT. */
static void gain_window(const lp2_gain *prev, const lp2_gain *cur, float *out)
{
    int G[64] = {0};
    int p = 0;
    for (int k = 0; k < cur->n; k++)
        while (p <= cur->loc[k] + 32)
            G[p++] = cur->lev[k] - 4;
    p = 0;
    for (int k = 0; k < prev->n; k++)
        while (p <= prev->loc[k])
            G[p++] += prev->lev[k] - 4;
    int o = 511, last = 0;
    for (int k = 63; k >= 0; k--) {          /* built backwards in 8-sample blocks */
        int e = -G[k];
        if (e == last) {
            float v = e < 0 ? (float)(1.0 / (1 << -e)) : (float)(1 << e);
            for (int m = 0; m < 8; m++)
                out[o - m] = v;
            o -= 8;
        } else {                             /* 8-step geometric interpolation */
            for (int m = 0; m < 8; m++) {
                int t = (m + 1) * e + (7 - m) * last + 176;
                out[o--] = lp2_pow2_8th[imin(imax(t, 0), 240)];
            }
        }
        last = e;
    }
}

/* ================================================================ MDCT (mdct256) */
static void mdct256(const float *x, int odd, float *out)
{
    double w[512];
    float z[256], tre[128], tim[128], re[128], im[128], y[256];
    for (int n = 0; n < 512; n++)
        w[n] = (double)lp2_mdct_win[n] * (double)x[n];
    for (int n = 0; n < 64; n++)
        z[n] = (float)((-w[384 + 2 * n]) - w[383 - 2 * n]);
    for (int m = 0; m < 128; m++)
        z[64 + m] = (float)(w[2 * m] - w[255 - 2 * m]);
    for (int n = 0; n < 64; n++)
        z[192 + n] = (float)(w[256 + 2 * n] + w[511 - 2 * n]);
    for (int j = 0; j < 128; j++) {
        double c = lp2_pre_cos[j], s = lp2_pre_sin[j], a = z[2 * j], b = z[2 * j + 1];
        tre[j] = (float)(c * a - s * b);
        tim[j] = (float)(s * a + c * b);
    }
    for (int j = 0; j < 128; j++) {
        re[j] = tre[lp2_bitrev7[j]];
        im[j] = tim[lp2_bitrev7[j]];
    }
    int groups = 64, span = 2;
    for (int stage = 0; stage < 7; stage++) { /* radix-2 decimation-in-time butterflies */
        int half = span / 2;
        for (int g = 0; g < groups; g++)
            for (int h = 0; h < half; h++) {
                double tc = lp2_fft_cos[h * groups], ts = lp2_fft_sin[h * groups];
                int top = g * span + h, bot = top + half;
                double rb = re[bot], ib = im[bot];
                double tr = tc * rb - ts * ib;
                double ti = ts * rb + tc * ib;
                double rt = re[top], it = im[top];
                re[bot] = (float)(rt - tr);
                im[bot] = (float)(it - ti);
                re[top] = (float)(tr + rt);
                im[top] = (float)(ti + it);
            }
        groups /= 2;
        span *= 2;
    }
    for (int k = 0; k < 128; k++) {
        double p1 = lp2_post1[k], p2 = lp2_post2[k], p3 = lp2_post3[k], p4 = lp2_post4[k];
        double Rk = re[k], Rm = re[127 - k], Ik = im[k], Im = im[127 - k];
        y[k] = (float)(((Rm * p2 + Ik * p3) + p1 * Rk) + p4 * Im);
        y[255 - k] = (float)(((p3 * Rk - p4 * Rm) - p1 * Ik) + p2 * Im);
    }
    for (int k = 0; k < 256; k++)
        out[odd ? 255 - k : k] = (float)((double)y[k] * 0.0078125);
}

/* ================================================================ per-channel front end */
void lp2_analysis_init(lp2_analysis *a)
{
    memset(a, 0, sizeof *a);
}

/* ChannelAnalysis.process.  After the call cur_gain holds this frame's gain points and prev_gain the
 * previous frame's. */
void lp2_analysis_process(lp2_analysis *a, const float *pcm, float *spec_gc, float *spec_plain)
{
    float lo[512], hi[512], band[4][256], buf[512], w[512];
    lp2_gain *prev = a->prev_gain, *cur = a->cur_gain;

    memcpy(prev, cur, sizeof a->cur_gain);
    qmf_split(pcm, 1024, a->qmf_delay[0], lo, hi);
    qmf_split(lo, 512, a->qmf_delay[1], band[0], band[1]);
    qmf_split(hi, 512, a->qmf_delay[2], band[3], band[2]);
    for (int b = 0; b < 4; b++) {
        memmove(a->hist[b], a->hist[b] + 256, 512 * sizeof(float));
        memcpy(a->hist[b] + 512, band[b], 256 * sizeof(float));
    }
    for (int b = 0; b < 4; b++)                                  /* gain control */
        band_gain_points(a->hist[b], &prev[b], &cur[b]);
    adjust_band0_gain(a->hist, prev, cur);

    int any = 0;
    for (int b = 0; b < 4; b++) {                                /* MDCT */
        const float *mid = a->hist[b] + 256;
        if (!cur[b].n && !prev[b].n) {
            memcpy(buf, a->overlap_gc[b], 256 * sizeof(float));
            memcpy(buf + 256, mid, 256 * sizeof(float));
        } else {
            any = 1;
            gain_window(&prev[b], &cur[b], w);
            for (int i = 0; i < 256; i++)
                buf[i] = (float)((double)a->overlap_gc[b][i] / (double)w[i]);
            for (int i = 0; i < 256; i++)
                buf[256 + i] = (float)((double)mid[i] / (double)w[256 + i]);
        }
        mdct256(buf, b & 1, spec_gc + 256 * b);
        memcpy(a->overlap_gc[b], mid, 256 * sizeof(float));
    }
    if (any) {
        for (int b = 0; b < 4; b++) {
            const float *mid = a->hist[b] + 256;
            memcpy(buf, a->overlap_plain[b], 256 * sizeof(float));
            memcpy(buf + 256, mid, 256 * sizeof(float));
            mdct256(buf, b & 1, spec_plain + 256 * b);
            memcpy(a->overlap_plain[b], mid, 256 * sizeof(float));
        }
    } else {
        memcpy(spec_plain, spec_gc, 1024 * sizeof(float));
        for (int b = 0; b < 4; b++)
            memcpy(a->overlap_plain[b], a->hist[b] + 256, 256 * sizeof(float));
    }
}
