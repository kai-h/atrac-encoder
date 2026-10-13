/* Internal declarations shared by the ATRAC3 encoder and decoder sources.
 *
 * Each source file mirrors one Python module of spec/reference, and the comments name the
 * corresponding function there.  Arithmetic follows the models exactly: `float` where the
 * model rounds to float32, `double` where it keeps a double (x87 register) value.  This only
 * holds when the compiler neither fuses multiply-adds nor uses excess precision, so the
 * build must use -ffp-contract=off (GCC/Clang) or /fp:precise (MSVC) and SSE2 on x86.
 */
#ifndef AT3_INTERNAL_H
#define AT3_INTERNAL_H

#include <stdint.h>
#include <string.h>
#include "atrac3.h"

#if defined(__FLT_EVAL_METHOD__) && __FLT_EVAL_METHOD__ != 0
#error "ATRAC3 needs FLT_EVAL_METHOD 0 (SSE2 on x86): build with -msse2 -mfpmath=sse"
#endif
/* The same as -ffp-contract=off, for builds that compile these sources without our flags. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#elif defined(_MSC_VER)
#pragma fp_contract(off)
#endif

#define AT3_FRAME_SAMPLES 1024

/* ---------------------------------------------------------------- helpers */
static inline uint32_t f32_bits(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }
static inline float f32_from_bits(uint32_t u) { float f; memcpy(&f, &u, 4); return f; }
/* Python's floor shift (v >> s) for negative v, portably. */
static inline int asr(int v, int s) { return v >= 0 ? v >> s : ~(~v >> s); }
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

/* MSB-first bit writer into a zero-filled buffer.  Writes past the end are dropped, so an
 * oversized unit can never corrupt memory (the models raise an error instead). */
typedef struct {
    uint8_t *buf;
    int cap_bits;
    int pos;
} bitwriter;

static inline void bw_put(bitwriter *bw, uint32_t v, int n)
{
    for (int k = n - 1; k >= 0; k--) {
        if (((v >> k) & 1) && bw->pos < bw->cap_bits)
            bw->buf[bw->pos >> 3] |= (uint8_t)(0x80 >> (bw->pos & 7));
        bw->pos++;
    }
}

/* ---------------------------------------------------------------- tables (tables.c) */
typedef struct { uint8_t code, len; } at3_code;

extern const int at3_sb_start[33];
extern const float at3_sf[64];
extern const at3_code at3_codes[2][8][64];   /* decoder-derived code tables */
extern const at3_code lp4_vlc[2][8][64];     /* Sony's LP4 encoder tables (incl. filler entries) */

extern const float lp2_t1[13], lp2_qmf_win[48], lp2_mdct_win[512];
extern const float lp2_pre_cos[128], lp2_pre_sin[128], lp2_fft_cos[64], lp2_fft_sin[64];
extern const int lp2_bitrev7[128];
extern const float lp2_post1[128], lp2_post2[128], lp2_post3[128], lp2_post4[128];
extern const float lp2_pow2_8th[241];

extern const int lp4_is_joint[6], lp4_bitrate[6], lp4_block_align[6];
extern const int lp4_channel_bytes0[6], lp4_channel_bytes1[6];
extern const int lp4_channel_hertz0[6], lp4_channel_hertz1[6];
extern const int lp4_sb_width[32];
extern const float lp4_qmf1[12][4], lp4_qmf2[24][4];
extern const int lp4_mdct_perm[128];
extern const float lp4_mdct_tw[128][4], lp4_mdct_fft_tw[127], lp4_mdct_sqrt_half;
extern const float lp4_gain_pow2[16];
extern const uint32_t lp4_gain_ramp_bits[8];
extern const double lp4_chsw_table[5], lp4_chwt_table[24];
extern const uint32_t lp4_sf_bits[64], lp4_qscale_bits[22];
extern const int lp4_qmask[8], lp4_min_sfsum[32], lp4_zero_len2[8], lp4_step_cost[8];
extern const int lp4_est_save[8], lp4_est_bits[8], lp4_est_thr[8], lp4_tone_penalty[5];
extern const int lp4_tone_sf_offset[32], lp4_sel1_len[16];
extern const double lp4_js_k_side, lp4_js_k_matrix, lp4_js_k_hyst, lp4_js_quiet;
extern const double lp4_js_ratio_gain, lp4_gc_floor, lp4_gc_min, lp4_gc_ratio, lp4_gc_sqrt2;
extern const double lp4_gc_expscale, lp4_gc_band0_max, lp4_tone_best_init;

extern const double dec_sf[64], dec_qmf_win[48], dec_window[512];
extern const double dec_weights[16][2], dec_inc_pow[31][8];

/* ================================================================ LP2 encoder ("algorithm 1") */
#define LP2_UNIT_BYTES 192
#define LP2_UNIT_BITS (LP2_UNIT_BYTES * 8)

typedef struct {             /* GainInfo: gain points of one QMF band for one frame */
    int n;
    int lev[32], loc[32];    /* at most 7 in valid streams */
    float peak;
} lp2_gain;

typedef struct {             /* ChannelAnalysis */
    float qmf_delay[3][46];
    float hist[4][768];
    float overlap_gc[4][256], overlap_plain[4][256];
    lp2_gain prev_gain[4], cur_gain[4];
} lp2_analysis;

typedef struct {             /* Component */
    int pos, sfi, len_code, sel, mode;
    int mant[8];
} lp2_comp;

typedef struct {             /* ToneGroup */
    int sel, len_code, mode;
    int flags[4];
    int ncell[16];
    lp2_comp cells[16][7];
} lp2_group;

#define LP2_MAX_GROUPS 4

typedef struct {             /* SoundUnit */
    int nbands;
    const lp2_gain *gain;    /* 4 bands */
    int ngroups;
    lp2_group groups[LP2_MAX_GROUPS];
    int tone_mode, nsub, mode;
    int sel[32], sfi[32];
    int q[1024];
    int ncomp, bits;
} lp2_unit;

typedef struct {
    int channels;
    lp2_analysis ana[2];
    lp2_unit unit;
    float spec_gc[1024], spec_plain[1024];
    int counts[3];
    int calls;
} lp2_encoder;

void lp2_analysis_init(lp2_analysis *a);
void lp2_analysis_process(lp2_analysis *a, const float *pcm, float *spec_gc, float *spec_plain);
void lp2_encode_unit(lp2_unit *u, const float *spec_gc, const float *spec_plain,
                       const lp2_gain *gain, const lp2_gain *prev_gain);
int lp2_mantissa_bits(int mode, int sel, const int *q, int n);
int lp2_gain_bits(const lp2_gain *gain, int nbands);
int lp2_packed_bits(const lp2_unit *u);
void lp2_pack_unit(const lp2_unit *u, uint8_t *out);

void lp2_encoder_init(lp2_encoder *e, int channels);
/* One library call: returns 1 if a frame was written to out (channels * 192 bytes). */
int lp2_encoder_call(lp2_encoder *e, const int16_t *block, int nsamples, uint8_t *out, int *end);

/* ================================================================ LP4 encoder ("algorithm 0") */
#define LP4_QMF_STATE_LEN 138

typedef struct {             /* BandGain */
    int loc[8];
    int lev[9];
    double env[32];
    double peak, last_block;
    int count;
} lp4_band_gain;

typedef struct {             /* Unit */
    int budget, js, nsub_max;    /* js: nonzero for the joint-stereo unit (unit 0's base budget) */
    float qmf_state[LP4_QMF_STATE_LEN];
    float prev_time[256][4];
    float mdct_state[128][4];
    lp4_band_gain gain[4];
    int prev_count0;
} lp4_unit;

typedef struct {             /* JointStereo */
    int limit;
    int dec_prev[4], dec[4];
    int w_prev, w;
    double r_prev[4], r[4];
    double energy[4][2];     /* this frame's band energies (L, R), float32 values */
} lp4_stereo;

typedef struct {             /* Comp: q[0..3], pos, sf */
    int r[6];
} lp4_comp;

typedef struct {             /* ToneGroup */
    int flag[4];
    int qstep, nvm1;
    int nlst[16];
    int lst[16][7];          /* indices into lp4_alloc_state.comps */
} lp4_group;

typedef struct {             /* AllocState */
    int sel[32], sf[32];
    int nsub, nbands, sel_mode, ngroups;
    lp4_group groups[2];
    lp4_comp comps[80];
    float spec[1024];
} lp4_alloc_state;

typedef struct {             /* UnitInfo */
    int budget, js, nsub_max;
    int gain_counts[4];
    int prev_count0, weight_sent;
} lp4_unit_info;

typedef struct {
    int frame_bytes, joint, bytes_per_ch;
    lp4_unit units[2];
    lp4_stereo stereo;
    int counter;
    float X[2][256][4];
    float spec[2][1024];
    lp4_alloc_state S;
    uint8_t buf[512 + 8];
} lp4_encoder;

void lp4_qmf_analysis(const float *x, float *state, float X[256][4]);
void lp4_envelope(float X[256][4], float env[4][32]);
void lp4_mdct(float prev[256][4], float state[128][4], const int first_levels[4], float *spec);

int lp4_gain_band(lp4_band_gain *bands, int b, const float *env_new, float prev[256][4],
                  int js_unit, const int *dec_prev, const int *dec_cur, lp4_unit *unit);

void lp4_stereo_init(lp4_stereo *js, int limit);
void lp4_stereo_analyse(lp4_stereo *js, float X0[256][4], float X1[256][4]);
void lp4_stereo_convert(lp4_stereo *js, int b, float X0[256][4], float X1[256][4]);
int lp4_split_bits(const float *spec0, const float *spec1, int nsub0, int nsub1,
                   int bytes_per_ch, int base0, int dec_prev0, int gaincount1_b0);

void lp4_alloc_init(lp4_alloc_state *S, const float *spec);
void lp4_allocate(lp4_alloc_state *S, const lp4_unit_info *U, const float *spec);
uint32_t lp4_quantise(float x, int sel, int sf);

int lp4_encoder_init(lp4_encoder *e, int kbps);
/* One library call (nbytes = 4 * samples, 4096 for a full block).  Returns 1 and writes
 * frame_bytes to out when a frame is produced, 0 when not, -1 when the flush is complete. */
int lp4_encoder_call(lp4_encoder *e, int nbytes, const int16_t *pcm, uint8_t *out);

#endif
