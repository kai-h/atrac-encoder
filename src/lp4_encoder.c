/* LP4 encoder ("algorithm 0"): sound-unit packing and the library call.  Mirrors
 * spec/reference/lp4_bitstream.py and the Encoder class of atrac3_encoder_lp4.py. */
#include "internal.h"

/* ---------------------------------------------------------------- packing */
static void pack_tones(bitwriter *bw, const lp4_alloc_state *S)
{
    bw_put(bw, S->sel_mode, 2);
    for (int gi = 0; gi < S->ngroups; gi++) {
        const lp4_group *g = &S->groups[gi];
        int v = g->flag[0];
        for (int k = 1; k < S->nbands; k++)
            v = v * 2 + g->flag[k];
        bw_put(bw, v, S->nbands);
        bw_put(bw, g->qstep + g->nvm1 * 8, 6);
        const at3_code *tab = lp4_vlc[S->sel_mode][g->qstep];
        for (int blk = 0; blk < S->nbands * 4; blk++) {
            if (!g->flag[blk >> 2])
                continue;
            bw_put(bw, g->nlst[blk], 3);
            for (int j = 0; j < g->nlst[blk]; j++) {
                const lp4_comp *rec = &S->comps[g->lst[blk][j]];
                bw_put(bw, rec->r[5], 6);
                bw_put(bw, rec->r[4] & 0x3f, 6);
                for (int k = 0; k <= g->nvm1; k++)
                    bw_put(bw, tab[rec->r[k] & 63].code, tab[rec->r[k] & 63].len);
            }
        }
    }
}

/* pack_unit: write one sound unit starting at byte byte_off of buf (zero filled, cap bytes).
 * js: second unit of a joint-stereo frame, preceded by the weight index w and the four
 * matrix selectors mat.  Returns the bit position after the unit. */
static int pack_unit(uint8_t *buf, int cap, int byte_off, const lp4_alloc_state *S,
                     const lp4_band_gain *gain, int js, int w, const int *mat)
{
    int head;
    if (byte_off + 2 > cap)                  /* unit 0 overran the frame (the model raises) */
        return byte_off * 8;
    if (js) {
        buf[byte_off++] = (uint8_t)(((w * 4 + mat[0]) * 4 + mat[1]) & 0xff);
        head = (mat[3] + mat[2] * 4) * 16 + 0x0c;     /* 2 matrix selectors + 2-bit id 3 */
    } else {
        head = 0xa0;                                  /* 6-bit id 0x28 */
    }
    buf[byte_off] = (uint8_t)((S->nbands + head - 1) & 0xff);
    bitwriter bw = {buf, cap * 8, byte_off * 8 + 8};
    for (int b = 0; b < S->nbands; b++) {
        bw_put(&bw, gain[b].count, 3);
        for (int k = 0; k < gain[b].count; k++)
            bw_put(&bw, gain[b].lev[k] * 32 + gain[b].loc[k], 9);
    }
    bw_put(&bw, S->ngroups, 5);
    if (S->ngroups)
        pack_tones(&bw, S);
    bw_put(&bw, S->nsub * 2 - 2, 6);                  /* (nsub - 1) and coding mode 0 = VLC */
    for (int i = 0; i < S->nsub; i++)
        bw_put(&bw, S->sel[i], 3);
    for (int i = 0; i < S->nsub; i++)
        if (S->sel[i])
            bw_put(&bw, S->sf[i], 6);
    for (int i = 0; i < S->nsub; i++) {
        int sel = S->sel[i];
        if (!sel)
            continue;
        const float *x = S->spec + at3_sb_start[i];
        int n = at3_sb_start[i + 1] - at3_sb_start[i];
        const at3_code *tab = lp4_vlc[0][sel];
        if (sel == 1) {
            for (int j = 0; j < n; j += 2) {
                uint32_t q0 = lp4_quantise(x[j], sel, S->sf[i]) & 3;
                uint32_t q1 = lp4_quantise(x[j + 1], sel, S->sf[i]) & 3;
                bw_put(&bw, tab[q0 * 4 + q1].code, tab[q0 * 4 + q1].len);
            }
        } else {
            for (int j = 0; j < n; j++) {
                uint32_t idx = lp4_quantise(x[j], sel, S->sf[i]) & (uint32_t)lp4_qmask[sel];
                bw_put(&bw, tab[idx].code, tab[idx].len);
            }
        }
    }
    return bw.pos;
}

/* ---------------------------------------------------------------- encoder */
static void bandwidth_subbands(int hz, int *nsub, int *idx)
{
    *idx = (hz << 11) / 44100;
    int i = 0;
    while (i < 0x1f && !(*idx <= at3_sb_start[i + 1]))
        i++;
    *nsub = i + 1;
}

int lp4_encoder_init(lp4_encoder *e, int kbps)
{
    int row = -1;
    for (int r = 0; r < 6; r++)
        if (lp4_bitrate[r] == kbps)
            row = r;
    if (row < 0)
        return -1;
    memset(e, 0, sizeof *e);
    e->frame_bytes = lp4_block_align[row];
    e->joint = lp4_is_joint[row];
    int cb[2] = {lp4_channel_bytes0[row], lp4_channel_bytes1[row]};
    int hz[2] = {lp4_channel_hertz0[row], lp4_channel_hertz1[row]};
    e->bytes_per_ch = (cb[0] + cb[1]) / 2;
    int limit = -1;
    for (int c = 0; c < 2; c++) {
        lp4_unit *u = &e->units[c];
        int js = e->joint && c == 1, nsub, idx;
        u->budget = cb[c] * 8 - (js ? 8 : 0) - 19;
        bandwidth_subbands(hz[c], &nsub, &idx);
        if (js) {
            js = e->units[0].budget;         /* Sony stores unit 0's budget as the flag */
            limit = imax(1, (idx + 0xff) >> 8);
        }
        u->js = js;
        u->nsub_max = nsub;
    }
    if (e->joint)
        lp4_stereo_init(&e->stereo, limit);
    return 0;
}

/* Encoder.encode_frame: analyse 1024 new samples per channel, emit one frame. */
static void encode_frame(lp4_encoder *e, float x[2][AT3_FRAME_SAMPLES], uint8_t *out)
{
    lp4_stereo *js = e->joint ? &e->stereo : NULL;
    float env[4][32], prev[256][4];

    for (int c = 0; c < 2; c++)
        lp4_qmf_analysis(x[c], e->units[c].qmf_state, e->X[c]);
    if (js) {
        lp4_stereo_analyse(js, e->X[0], e->X[1]);
        for (int b = 0; b < 4; b++)
            lp4_stereo_convert(js, b, e->X[0], e->X[1]);
    }
    for (int c = 0; c < 2; c++) {
        lp4_unit *u = &e->units[c];
        int fl[4];
        lp4_envelope(e->X[c], env);
        memcpy(prev, u->prev_time, sizeof prev);
        for (int b = 3; b >= 0; b--) {
            int n = lp4_gain_band(u->gain, b, env[b], prev, u->js != 0, js ? js->dec_prev : NULL,
                                  js ? js->dec : NULL, u);
            u->gain[b].count = n;
            u->gain[b].lev[n] = 4;
        }
        for (int b = 0; b < 4; b++)
            fl[b] = u->gain[b].lev[0];
        lp4_mdct(prev, u->mdct_state, fl, e->spec[c]);
        memcpy(u->prev_time, e->X[c], sizeof u->prev_time);
    }
    if (js) {
        lp4_unit *u0 = &e->units[0], *u1 = &e->units[1];
        u0->budget = lp4_split_bits(e->spec[0], e->spec[1], u0->nsub_max, u1->nsub_max,
                                    e->bytes_per_ch, u1->js, js->dec_prev[0], u1->gain[0].count);
    }
    int cap = e->frame_bytes + 8;
    memset(e->buf, 0, cap);
    int offset = 0;
    for (int c = 0; c < 2; c++) {
        lp4_unit *u = &e->units[c];
        lp4_alloc_state *S = &e->S;
        lp4_unit_info info = {u->budget, u->js, u->nsub_max,
                              {u->gain[0].count, u->gain[1].count, u->gain[2].count, u->gain[3].count},
                              u->prev_count0, js ? js->w_prev : 7};
        if (u->budget < 40) {
            lp4_alloc_init(S, e->spec[c]);
            S->nbands = S->nsub = 1;
            S->ngroups = 0;
            u->gain[0].count = 0;
        } else {
            lp4_allocate(S, &info, e->spec[c]);
        }
        int is_js = js && u->js;
        int end = pack_unit(e->buf, cap, offset, S, u->gain, is_js, is_js ? js->w_prev : 0,
                            is_js ? js->dec_prev : NULL);
        if (js) {
            if (c == 0) {
                offset = (end + 7) >> 3;
                e->units[1].budget = (e->bytes_per_ch * 2 - offset) * 8 - 0x1b;
            } else {                         /* unit 1 is stored byte-reversed from the end */
                for (int i = offset, j = e->frame_bytes - 1; i < j; i++, j--) {
                    uint8_t t = e->buf[i];
                    e->buf[i] = e->buf[j];
                    e->buf[j] = t;
                }
            }
        } else {
            offset += e->bytes_per_ch;
        }
    }
    memcpy(out, e->buf, e->frame_bytes);
}

/* a3e0_encode / a3e0_flush: one library call.  A block shorter than 1024 samples marks the
 * end of the input; afterwards calls with nbytes 0 flush the pending frames. */
int lp4_encoder_call(lp4_encoder *e, int nbytes, const int16_t *pcm, uint8_t *out)
{
    float x[2][AT3_FRAME_SAMPLES];
    e->counter++;
    if (nbytes != 4096 && e->counter > 0)
        e->counter = nbytes < 0xdd8 ? -3 : -4;
    if (e->counter == 0)
        return -1;
    memset(x, 0, sizeof x);
    for (int i = 0; i < nbytes / 4; i++) {
        x[0][i] = pcm[2 * i];
        x[1][i] = pcm[2 * i + 1];
    }
    encode_frame(e, x, out);
    return e->counter == 1 ? 0 : 1;
}
