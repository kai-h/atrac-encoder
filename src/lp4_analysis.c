/* Signal analysis of the LP4 encoder: QMF band split, block envelope and the band MDCT.
 * Mirrors spec/reference/lp4_analysis.py.  All three stages are SSE single-precision code in
 * Sony's encoder: every operation rounds to float, in the order given here.
 * Layout: a QMF frame is X[256][4] (time index, band). */
#include <math.h>
#include "internal.h"

/* qmf_analysis: 1024 PCM samples -> 4 bands of 256 samples; state[138] carried over. */
void lp4_qmf_analysis(const float *x, float *state, float X[256][4])
{
    float xh[46 + 1024];
    float s1all[23 + 256][4];
    memcpy(xh, state + 92, 46 * sizeof(float));
    memcpy(xh + 46, x, 1024 * sizeof(float));
    memcpy(s1all, state, 92 * sizeof(float));

    /* stage 1: SA uses offsets +4, 0, -4, ..., -40; SB the same shifted by +2 */
    for (int t = 0; t < 256; t++) {
        const float *v = xh + 40 + 4 * t;    /* vec(o)[t][j] = v[o + j] */
        float sa[4], sb[4];
        for (int j = 0; j < 4; j++) {
            float s = lp4_qmf1[11][j] * v[4 + j] + lp4_qmf1[0][j] * v[-40 + j];
            for (int k = 9; k >= 0; k--)
                s = s + lp4_qmf1[10 - k][j] * v[-4 * k + j];
            sa[j] = s;
            s = lp4_qmf1[11][j] * v[6 + j] + lp4_qmf1[0][j] * v[2 - 40 + j];
            for (int k = 9; k >= 0; k--)
                s = s + lp4_qmf1[10 - k][j] * v[2 - 4 * k + j];
            sb[j] = s;
        }
        float u0 = sa[2] + sa[0], u1 = sb[2] + sb[0], u2 = sa[3] + sa[1], u3 = sb[3] + sb[1];
        float *s1 = s1all[23 + t];
        s1[0] = u2 + u0;
        s1[1] = u3 + u1;
        s1[2] = u2 - u0;
        s1[3] = u3 - u1;
    }
    /* stage 2: 24-tap filter over stage-1 vectors, then the final butterfly */
    for (int t = 0; t < 256; t++) {
        float acc[4];
        for (int j = 0; j < 4; j++) {
            float a = lp4_qmf2[23][j] * s1all[23 + t][j];
            for (int k = 22; k >= 0; k--)
                a = a + lp4_qmf2[22 - k][j] * s1all[22 - k + t][j];
            acc[j] = a;
        }
        X[t][0] = acc[0] + acc[1];
        X[t][1] = acc[1] - acc[0];
        X[t][2] = acc[3] - acc[2];
        X[t][3] = acc[2] + acc[3];
    }
    memcpy(state, s1all[256], 92 * sizeof(float));
    memcpy(state + 92, xh + 1024, 46 * sizeof(float));
}

/* envelope: peak magnitude of each block of 8 samples per band. */
void lp4_envelope(float X[256][4], float env[4][32])
{
    for (int b = 0; b < 4; b++)
        for (int i = 0; i < 32; i++) {
            float m = fabsf(X[8 * i][b]);
            for (int k = 1; k < 8; k++)
                if (fabsf(X[8 * i + k][b]) > m)
                    m = fabsf(X[8 * i + k][b]);
            env[b][i] = m;
        }
}

/* mdct: forward MDCT of all four bands.  prev = the previous frame's 256 band samples (after
 * gain modification), state = folded first-half data of the previous call (updated).
 * spec receives 4 x 256 coefficients, odd bands reversed in frequency. */
void lp4_mdct(float prev[256][4], float state[128][4], const int first_levels[4], float *spec)
{
    float W[256][4], hi2[64][4], lo2[64][4];
    float g[4];
    for (int b = 0; b < 4; b++)
        g[b] = lp4_gain_pow2[first_levels[b] & 15];

    for (int m = 0; m < 128; m++) {
        int k = lp4_mdct_perm[m];
        const float *tw = lp4_mdct_tw[k];
        for (int b = 0; b < 4; b++) {
            float xa = prev[128 + k][b], xb = prev[127 - k][b];
            float s2g = state[m][b] * g[b];
            float t = tw[2] * (xb * tw[0] + xa);
            float d = tw[3] * (t - s2g);
            state[m][b] = tw[1] * (xa * tw[0] - xb);
            W[2 * m][b] = (t + s2g) + d;
            W[2 * m + 1][b] = d;
        }
    }
    /* first radix stage (groups of 4 vectors) */
    for (int gi = 0; gi < 64; gi++) {
        float T = lp4_mdct_fft_tw[2 * gi];
        for (int b = 0; b < 4; b++) {
            float a = W[4 * gi][b], bb = W[4 * gi + 1][b], c = W[4 * gi + 2][b], e = W[4 * gi + 3][b];
            float x3 = T * (bb - e), x0 = T * (a - c);
            W[4 * gi][b] = (a + c) + x3;
            W[4 * gi + 1][b] = (bb + e) + x0;
            W[4 * gi + 2][b] = x0;
            W[4 * gi + 3][b] = x3;
        }
    }
    /* further stages, block size 2n for n = 4 .. 64 */
    for (int n = 4; n != 128; n *= 2) {
        int nb = 256 / (2 * n);
        for (int bi = 0; bi < nb; bi++) {
            float twd = lp4_mdct_fft_tw[n / 2 - 1 + n * bi];
            float (*V)[4] = W + bi * 2 * n;
            for (int r = 0; r < n; r++)
                for (int b = 0; b < 4; b++) {
                    hi2[r][b] = twd * (V[r][b] - V[n + r][b]);
                    lo2[r][b] = V[r][b] + V[n + r][b];
                }
            for (int r = 0; r < n; r++)
                for (int b = 0; b < 4; b++) {
                    V[r][b] = hi2[n - 1 - r][b] + lo2[r][b];
                    V[n + r][b] = hi2[r][b];
                }
        }
    }
    for (int i = 0; i < 128; i++)
        for (int b = 0; b < 4; b++) {
            float lo = W[i][b], hi = W[128 + i][b];
            W[128 + i][b] = lp4_mdct_sqrt_half * (lo - hi);
            W[i][b] = lo + hi;
        }
    for (int b = 0; b < 4; b++) {
        float *s = spec + 256 * b;
        for (int i = 0; i < 128; i++) {
            float wp = W[i][b] + W[255 - i][b];
            if (b % 2 == 0) {
                s[i] = wp;
                s[255 - i] = W[255 - i][b];
            } else {                         /* odd bands come out frequency-reversed */
                s[i] = W[255 - i][b];
                s[255 - i] = wp;
            }
        }
    }
}
