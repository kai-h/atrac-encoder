/* Public encoder: mode selection, at3tool's framing and the .at3 header.
 * Mirrors encode_pcm / encode_file (atrac3_encoder_lp2.py) and encode_wav
 * (atrac3_encoder_lp4.py); see spec/ATRAC3.md sections 2, 3 and 6.4. */
#include <stdlib.h>
#include "internal.h"

#define PRIME_ZEROS 955                      /* 1024 - 69: zeros the tool feeds before the source */

struct at3_encoder {
    int channels;            /* input channels */
    int codec_channels;      /* channels handed to the library (52 kbps mono runs as stereo) */
    int frame_bytes;         /* output bytes per frame */
    int joint;
    int16_t block[AT3_FRAME_SAMPLES * 2];
    int fill;
    int frames_seen;         /* frames returned by the library; the tool drops the first */
    int done;
    uint8_t frame[512];
    lp2_encoder *lp2;
    lp4_encoder *lp4;
};

at3_encoder *at3_encoder_create(int kbps, int channels)
{
    at3_encoder *e = calloc(1, sizeof *e);
    if (!e)
        return NULL;
    e->channels = channels;
    if ((kbps == 132 && channels == 2) || (kbps == 66 && channels == 1)) {
        /* algorithm 1; mono is channel 0 of the stereo encoder (channels are independent) */
        e->codec_channels = channels;
        e->frame_bytes = LP2_UNIT_BYTES * channels;
        e->lp2 = malloc(sizeof *e->lp2);
        if (e->lp2)
            lp2_encoder_init(e->lp2, channels);
    } else if ((kbps == 66 || kbps == 105) && channels == 2) {
        e->codec_channels = 2;
        e->lp4 = malloc(sizeof *e->lp4);
        if (e->lp4)
            lp4_encoder_init(e->lp4, kbps);
        e->frame_bytes = e->lp4 ? e->lp4->frame_bytes : 0;
        e->joint = kbps == 66;
    } else if (kbps == 52 && channels == 1) {
        /* the tool duplicates mono input into the 105 kbps stereo encoder, keeps unit 0 */
        e->codec_channels = 2;
        e->frame_bytes = 152;
        e->lp4 = malloc(sizeof *e->lp4);
        if (e->lp4)
            lp4_encoder_init(e->lp4, 105);
    }
    if (!e->lp2 && !e->lp4) {
        free(e);
        return NULL;
    }
    e->fill = PRIME_ZEROS;                   /* block[] starts zeroed */
    return e;
}

void at3_encoder_destroy(at3_encoder *e)
{
    if (!e)
        return;
    free(e->lp2);
    free(e->lp4);
    free(e);
}

int at3_encoder_frame_bytes(const at3_encoder *e) { return e->frame_bytes; }
int at3_encoder_joint_stereo(const at3_encoder *e) { return e->joint; }

size_t at3_encoder_max_output(const at3_encoder *e, size_t nsamples)
{
    return ((e->fill + nsamples) / AT3_FRAME_SAMPLES) * (size_t)e->frame_bytes;
}

/* One library call on the pending block (n samples).  Returns the bytes written to out. */
static size_t library_call(at3_encoder *e, int n, uint8_t *out)
{
    int r, end = 0;
    if (e->lp2) {
        r = lp2_encoder_call(e->lp2, e->block, n, e->frame, &end);
    } else {
        r = lp4_encoder_call(e->lp4, n * 4, e->block, e->frame);
        end = r < 0;
    }
    if (end)
        e->done = 1;
    if (r != 1 || e->frames_seen++ == 0)
        return 0;
    memcpy(out, e->frame, e->frame_bytes);
    return e->frame_bytes;
}

size_t at3_encoder_write(at3_encoder *e, const int16_t *pcm, size_t nsamples, uint8_t *out)
{
    size_t written = 0;
    if (e->done)
        return 0;
    for (size_t s = 0; s < nsamples; s++) {
        for (int c = 0; c < e->codec_channels; c++)
            e->block[e->fill * e->codec_channels + c] = pcm[s * e->channels + (e->channels == 1 ? 0 : c)];
        if (++e->fill == AT3_FRAME_SAMPLES) {
            written += library_call(e, AT3_FRAME_SAMPLES, out + written);
            e->fill = 0;
        }
    }
    return written;
}

size_t at3_encoder_finish(at3_encoder *e, uint8_t *out)
{
    size_t written = 0;
    int n = e->fill;
    while (!e->done) {                       /* the short block, then flush calls until done */
        written += library_call(e, n, out + written);
        n = 0;
    }
    e->fill = 0;
    return written;
}

/* ---------------------------------------------------------------- .at3 header */
static uint8_t *put16(uint8_t *p, unsigned v) { p[0] = v & 0xff; p[1] = (v >> 8) & 0xff; return p + 2; }
static uint8_t *put32(uint8_t *p, uint32_t v) { put16(p, v & 0xffff); put16(p + 2, v >> 16); return p + 4; }
static uint8_t *put4cc(uint8_t *p, const char *s) { memcpy(p, s, 4); return p + 4; }

void at3_write_header(uint8_t *out, int channels, int frame_bytes, int joint_stereo,
                      uint32_t total_samples, uint32_t data_bytes)
{
    uint8_t *p = out;
    p = put4cc(p, "RIFF");
    p = put32(p, AT3_HEADER_BYTES - 8 + data_bytes);
    p = put4cc(p, "WAVE");
    p = put4cc(p, "fmt ");
    p = put32(p, 32);
    p = put16(p, 0x270);
    p = put16(p, channels);
    p = put32(p, 44100);
    p = put32(p, (uint32_t)(44100.0 * frame_bytes / 1024.0 + 0.5));
    p = put16(p, frame_bytes);
    p = put16(p, 0);
    p = put16(p, 14);
    p = put16(p, 1);
    p = put32(p, channels * 2048);
    p = put16(p, joint_stereo);
    p = put16(p, joint_stereo);
    p = put16(p, 1);
    p = put16(p, 0);
    p = put4cc(p, "fact");
    p = put32(p, 12);
    p = put32(p, total_samples);
    p = put32(p, 1024);
    p = put32(p, 1024);
    p = put4cc(p, "data");
    put32(p, data_bytes);
}
