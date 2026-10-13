/* Gain control (transient pre-processing) of the LP4 encoder (x87 code: doubles with float
 * stores).  Mirrors spec/reference/lp4_gain.py. */
#include "internal.h"

/* floor(log2(ratio * sqrt2)) from the float32 exponent field, i.e. round(log2(ratio)). */
static int log2_level(double ratio)
{
    return (int)(f32_bits((float)(ratio * lp4_gc_sqrt2)) >> 23) - 127;
}

/* detect: gain points for one band from A = 64 block peaks [previous frame 32 | new 32].
 * mode: band index for ordinary units, or -1 / 0 / 5 for the joint-stereo unit. */
static int detect(lp4_band_gain *G, const double *A, int mode)
{
    double bm[8];
    for (int j = 0; j < 8; j++) {
        bm[j] = A[4 * j];
        for (int k = 1; k < 4; k++)
            if (A[4 * j + k] > bm[j])
                bm[j] = A[4 * j + k];
    }
    double prev_last = G->last_block;
    G->last_block = bm[7];
    int end = mode <= 0 ? 64 : (8 - mode) * 8;
    double ref = bm[7];                      /* look-ahead reference level */
    for (int i = 32; i < end; i++)
        if (A[i] > ref)
            ref = A[i];
    if (ref < lp4_gc_floor)
        ref = lp4_gc_floor;
    double thr = (float)(lp4_gc_ratio * ref);
    int budget = 4;                          /* total attenuation available (levels) */
    int slot = 7;                            /* decay points are stored from slot 6 down */
    G->loc[7] = 32;
    double blk[9];
    blk[0] = prev_last;
    for (int j = 0; j < 8; j++)
        blk[j + 1] = bm[j];
    /* part 1: scan backwards for decays (blocks louder than 1.85x what follows them) */
    for (int e = 8; e >= 0; e--) {
        double v = blk[e];
        if (v < ref)
            continue;
        if (v > lp4_gc_min && v > thr) {
            int loc = 4 * e;
            if (e != 0) {
                if (A[4 * e] >= thr || A[4 * e - 1] >= thr)
                    loc = 4 * e;
                else if (A[4 * e - 2] >= thr)
                    loc = 4 * e - 1;
                else
                    loc = 4 * e - 2;
            }
            G->loc[--slot] = loc;
            int ex = imin(log2_level(v / ref), budget);
            budget -= ex;
            G->lev[slot] = -ex;              /* relative step, made absolute below */
            if (budget <= 0 || slot == 5)
                break;
        }
        ref = v;
        thr = (float)(v * lp4_gc_ratio);
    }
    /* part 2: scan forwards for attacks before the first decay point */
    int n = 0;
    int room = imin(15, 0x83 - (int)(f32_bits((float)(G->peak * lp4_gc_expscale)) >> 23)) - budget;
    if (room > 0) {
        double fac = mode == -1 ? 2.0 : (double)1.6f;
        ref = G->peak;
        if (ref < A[0])
            ref = A[0];
        thr = fac * ref;
        for (int i = 0; i < G->loc[slot]; i++) {
            double a = A[i + 1];
            if (a < ref)
                continue;
            if (a > lp4_gc_min && a > thr) {
                G->loc[n] = i;
                int ex = log2_level(a / ref);
                if (n >= 1 && G->loc[n - 1] == i - 1 && ex >= G->lev[n - 1]) {
                    n--;                     /* merge with an adjacent attack */
                    G->loc[n] = i;
                    room += G->lev[n];
                    ex += G->lev[n];
                }
                ex = imin(ex, room);
                room -= ex;
                G->lev[n++] = ex;
                if (n == slot || room <= 0)
                    break;
            }
            ref = a;
            thr = a * fac;
        }
    }
    while (slot < 7) {                       /* append the decay points */
        G->lev[n] = G->lev[slot];
        G->loc[n++] = G->loc[slot++];
    }
    int acc = 4;                             /* cumulative levels, last point -> 4 */
    for (int k = n - 1; k >= 0; k--) {
        acc += G->lev[k];
        G->lev[k] = acc;
    }
    return n;
}

/* apply_gain: multiply a band signal by the gain curve: constant 2**(lev-4) up to and
 * including sample 8*loc, then a 7-sample geometric ramp towards the next level. */
static void apply_gain(const lp4_band_gain *G, int n, float prev[256][4], int b)
{
    int lev[9];
    for (int i = 0; i < n; i++)
        lev[i] = G->lev[i];
    lev[n] = 4;
    int pos = 0;
    for (int i = 0; i <= n; i++) {
        int c = lev[i];
        if (c == 4) {
            if (i == n)
                return;
            pos = G->loc[i] * 8 + 1;
        } else {
            int end = G->loc[i] * 8;
            float p = lp4_gain_pow2[c & 15];
            while (pos < end) {
                for (int k = 0; k < 8; k++)
                    prev[pos + k][b] = prev[pos + k][b] * p;
                pos += 8;
            }
            prev[pos][b] = prev[pos][b] * p;
            pos++;
        }
        int d = lev[i + 1] - c;
        int v = c * 8 + d;
        for (int j = 0; j < 7; j++) {        /* g = 2**(v/8 - 4) */
            float g = f32_from_bits(((uint32_t)v << 20) + lp4_gain_ramp_bits[v & 7]);
            prev[pos][b] = prev[pos][b] * g;
            pos++;
            v += d;
        }
    }
}

/* gain_band: gain control of band b of one unit.  Modifies prev[.][b] (the previous frame's
 * samples, about to be transformed) and bands[b]; returns the point count. */
int lp4_gain_band(lp4_band_gain *bands, int b, const float *env_new, float prev[256][4],
                  int js_unit, const int *dec_prev, const int *dec_cur, lp4_unit *unit)
{
    lp4_band_gain *G = &bands[b];
    double A[64];
    double peak = lp4_gc_floor;
    for (int i = 0; i < 32; i++) {
        A[i] = G->env[i];
        A[32 + i] = env_new[i];
        if (G->env[i] > peak)
            peak = G->env[i];
    }
    for (int i = 0; i < 32; i++)
        G->env[i] = env_new[i];
    int mode = b;
    if (js_unit) {
        if (dec_prev[b] != dec_cur[b])
            mode = 5;
        else if (dec_prev[b] != 3)
            mode = 0;
        else
            mode = -1;
    }
    int n = detect(G, A, mode);
    G->peak = peak;
    if (b == 0)
        unit->prev_count0 = G->count;
    if (n == 0) {
        /* band 0 may copy a strong attack found in band 1 (bass smearing protection) */
        const lp4_band_gain *g1 = &bands[1];
        if (b == 0 && !js_unit && G->count == 0 && g1->count != 0 && peak <= lp4_gc_band0_max) {
            int mx = 4, mn = 4;
            for (int k = 0; k < g1->count; k++) {
                mx = imax(mx, g1->lev[k]);
                mn = imin(mn, g1->lev[k]);
            }
            if (mx - mn < 2)
                return 0;
            int l0 = g1->loc[0];
            for (int k = 0; k <= l0; k++)
                if (A[k] > lp4_gc_band0_max)
                    return 0;
            n = 1;
            G->loc[0] = l0;
            G->lev[0] = 5;
        } else {
            return 0;
        }
    }
    apply_gain(G, n, prev, b);
    return n;
}
