/* Joint-stereo stage of the LP4 encoder (x87 code: doubles with float stores): JointStereo.analyse,
 * the per-band conversion and split_bits.  Mirrors spec/reference/lp4_stereo.py. */
#include <math.h>
#include <stdlib.h>
#include "internal.h"

void lp4_stereo_init(lp4_stereo *js, int limit)
{
    memset(js, 0, sizeof *js);
    js->limit = limit;
    for (int b = 0; b < 4; b++) {
        js->dec[b] = -3;
        js->r_prev[b] = js->r[b] = 1.0;
    }
    js->w = 15;
}

/* JointStereo.analyse: band energies (sum of |x|), per-band decisions, weight index. */
void lp4_stereo_analyse(lp4_stereo *js, float X0[256][4], float X1[256][4])
{
    double acc0 = (double)0.1f, acc1 = acc0;
    for (int b = 0; b < 4; b++) {
        double s0 = 0.0, s1 = 0.0;
        for (int t = 0; t < 256; t++) {
            s0 += fabs((double)X0[t][b]);
            s1 += fabs((double)X1[t][b]);
        }
        double e0 = (float)s0, e1 = (float)s1;
        js->energy[b][0] = e0;
        js->energy[b][1] = e1;
        int prev = js->dec[b];
        double k0 = lp4_js_k_side, k1 = lp4_js_k_side;
        if (b < js->limit) {                 /* hysteresis: easier to stay in a one-sided mode */
            k0 = k1 = lp4_js_k_matrix;
            if (prev == 0)
                k0 = lp4_js_k_hyst;
            else if (prev == 1)
                k1 = lp4_js_k_hyst;
        }
        int d;
        if (k0 * s0 < s1)
            d = 0;
        else if (k1 * s1 < s0)
            d = 1;
        else if (s1 + s0 < lp4_js_quiet)
            d = -3;
        else
            d = 3;
        js->dec_prev[b] = abs(prev);
        js->dec[b] = d;
        if (b >= js->limit && (d == 3 || d == -3)) {
            acc0 = acc0 * 16.0 + e0;
            acc1 = acc1 * 16.0 + e1;
        }
    }
    js->w_prev = js->w;
    double tot = (float)(acc1 + acc0), mn;
    int side;
    if (acc0 <= acc1) {
        mn = acc0;
        side = 0;
    } else {
        mn = (float)acc1;
        side = 8;
    }
    int w = (int)(f32_bits((float)(mn / tot * 14.0 + 12582912.0)) & 7) + side;
    if (w == 0 || w == 8)
        w++;
    js->w = w;
}

/* Apply the matrix (b < limit) or the intensity weighting (b >= limit) to band b in place. */
void lp4_stereo_convert(lp4_stereo *js, int b, float X0[256][4], float X1[256][4])
{
    int cur = abs(js->dec[b]);
    if (b >= js->limit) {
        int wc = js->w, wp = js->w_prev;
        if (cur == 0) {
            wc = (wc ^ 8) + 8;
            wp = (wp ^ 8) + 8;
        } else if (cur == 1) {
            wc += 8;
            wp += 8;
        } else if (cur == 3) {
            wc &= 7;
            wp &= 7;
        }
        double tc = lp4_chwt_table[wc], tp = lp4_chwt_table[wp];
        double step = (tp - tc) * 0.125;
        for (int i = 0; i < 8; i++) {
            int k = i - 8;
            X0[i][b] = (float)(((double)X1[i][b] + (double)X0[i][b]) / (tc - k * step));
        }
        double inv = 1.0 / tc;
        for (int i = 8; i < 256; i++)
            X0[i][b] = (float)(((double)X1[i][b] + (double)X0[i][b]) * inv);
        return;
    }
    double e0 = js->energy[b][0], e1 = js->energy[b][1];
    int prv = js->dec_prev[b];
    double c0 = lp4_chsw_table[cur], c1 = lp4_chsw_table[cur + 1];
    double p0 = lp4_chsw_table[prv], p1 = lp4_chsw_table[prv + 1];
    double r = 1.0;
    if (cur == 3) {                          /* side gain from the energy imbalance */
        if (!(e1 <= e0 * 3.94) && !(e0 * 256.0 <= e1))
            r = (float)((e1 - e0) * lp4_js_ratio_gain / (e1 + e0));
        else if (!(e0 <= e1 * 3.94) && !(e1 * 256.0 <= e0))
            r = (float)((e0 - e1) * lp4_js_ratio_gain / (e1 + e0));
        else if (e1 == 0.0 || e0 == 0.0)
            r = 1.0;
        else if (e1 * 256.0 <= e0 || e0 * 256.0 <= e1)
            r = (float)(5.0 / 3.0);
    }
    double q = js->r_prev[b];
    if (!(r - q < 0.05))                     /* limit the change to 0.05 per frame */
        r = (float)(q + 0.05);
    if (!(q - r < 0.05))
        r = (float)(q - 0.05);
    js->r[b] = r;
    double c1r = c1 * r;
    double dc = (float)((p0 - c0) * 0.125);
    double dw = (p1 * q - c1r) * 0.125;
    js->r_prev[b] = r;
    for (int i = 0; i < 8; i++) {
        int k = i - 8;
        double x0 = X0[i][b], x1 = X1[i][b];
        double m = (x1 + x0) * 0.5;
        X1[i][b] = (float)((x0 - (c0 - dc * k) * m) / (c1r - k * dw));
        X0[i][b] = (float)m;
    }
    double inv = 1.0 / c1r;
    for (int i = 8; i < 256; i++) {
        double x0 = X0[i][b], x1 = X1[i][b];
        double m = (x1 + x0) * 0.5;
        X1[i][b] = (float)((x0 - m * c0) * inv);
        X0[i][b] = (float)m;
    }
}

/* split_bits: bits granted to the first sound unit in joint-stereo mode. */
int lp4_split_bits(const float *spec0, const float *spec1, int nsub0, int nsub1,
                   int bytes_per_ch, int base0, int dec_prev0, int gaincount1_b0)
{
    int total = bytes_per_ch * 16 - 0x3b;
    double sum0 = 0.0, sum1 = 0.0;
    for (int i = imin(nsub0 - 1, 18); i >= 0; i--) {
        double mx = 0.0;
        for (int k = at3_sb_start[i]; k < at3_sb_start[i + 1]; k++)
            if (!(mx >= fabs((double)spec0[k])))
                mx = fabs((double)spec0[k]);
        sum0 += lp4_sb_width[i] * mx;
    }
    sum0 = (float)sum0;
    for (int i = 0; i < nsub1; i++) {
        double mx = 0.0;
        for (int k = at3_sb_start[i]; k < at3_sb_start[i + 1]; k++)
            if (!(mx >= fabs((double)spec1[k])))
                mx = fabs((double)spec1[k]);
        sum1 += lp4_sb_width[i] * mx;
    }
    double t = (float)(sum1 + sum1);
    if (dec_prev0 != 3)
        t = (float)(t + t);
    if (t < 1.0 || (dec_prev0 != 3 && t < 4.0))
        return total;
    int spare = total - base0;
    if (gaincount1_b0 == 0 && ((t * 3.0 < sum0 && dec_prev0 == 3) || t * 6.0 < sum0)) {
        int x;
        if (!(t * 12.0 > sum0)) {
            x = spare;
        } else {
            if (!(t * 10.0 >= sum0))
                spare = (int)(spare * 1.5);
            else if (!(t * 5.0 >= sum0))
                spare = (int)(spare * 1.25);
            x = (int)((sum0 - t) / (sum0 + sum0) * spare);
        }
        return base0 + (x & ~7);
    }
    if (t >= sum0 || dec_prev0 != 3)
        return base0;
    int x = (int)((sum0 - t) / (sum0 * 10.0) * spare);
    return base0 + (x & ~7);
}
