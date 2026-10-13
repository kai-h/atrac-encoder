/* Whole-file helpers shared by the command-line tool and the WebAssembly module:
 * 16-bit PCM WAV <-> .at3 (RIFF/WAVE, format 0x270), entirely in memory. */
#include <stdlib.h>
#include "internal.h"

static uint32_t rd16(const uint8_t *p) { return p[0] | (uint32_t)p[1] << 8; }
static uint32_t rd32(const uint8_t *p) { return rd16(p) | rd16(p + 2) << 16; }

typedef struct {
    const uint8_t *fmt, *fact, *data;
    uint32_t fmt_len, fact_len, data_len;
    int has_smpl;
} riff;

static const char *parse_riff(const uint8_t *b, size_t len, riff *r)
{
    memset(r, 0, sizeof *r);
    if (len < 12 || memcmp(b, "RIFF", 4) || memcmp(b + 8, "WAVE", 4))
        return "not a RIFF/WAVE file";
    size_t pos = 12;
    while (pos + 8 <= len) {
        const uint8_t *id = b + pos;
        uint32_t size = rd32(b + pos + 4);
        size_t avail = len - pos - 8;
        uint32_t n = size < avail ? size : (uint32_t)avail;     /* tolerate truncated files */
        const uint8_t *body = b + pos + 8;
        if (!memcmp(id, "fmt ", 4)) {
            r->fmt = body;
            r->fmt_len = n;
        } else if (!memcmp(id, "fact", 4)) {
            r->fact = body;
            r->fact_len = n;
        } else if (!memcmp(id, "data", 4) && !r->data) {
            r->data = body;
            r->data_len = n;
        } else if (!memcmp(id, "smpl", 4)) {
            r->has_smpl = 1;
        }
        if (size >= avail)      /* chunk runs to the end of the file */
            break;
        pos += 8 + (size_t)size + (size & 1);
    }
    if (!r->fmt || r->fmt_len < 16)
        return "missing fmt chunk";
    if (!r->data)
        return "missing data chunk";
    return NULL;
}

void at3_free(void *p)
{
    free(p);
}

const char *at3_encode_wav(const uint8_t *wav, size_t wav_len, int kbps,
                           uint8_t **out, size_t *out_len, int *warnings)
{
    riff r;
    const char *err = parse_riff(wav, wav_len, &r);
    if (err)
        return err;
    int tag = rd16(r.fmt), channels = rd16(r.fmt + 2), bits = rd16(r.fmt + 14);
    if (tag == 0xfffe && r.fmt_len >= 26)    /* WAVE_FORMAT_EXTENSIBLE: subformat GUID */
        tag = rd16(r.fmt + 24);
    if (tag != 1 || bits != 16)
        return "input must be 16-bit PCM";
    if (rd32(r.fmt + 4) != 44100)
        return "input must be 44.1 kHz";
    if (channels != 1 && channels != 2)
        return "input must be mono or stereo";
    if (kbps == 0)
        kbps = channels == 2 ? 132 : 66;
    at3_encoder *enc = at3_encoder_create(kbps, channels);
    if (!enc)
        return channels == 2 ? "stereo input supports 132 (LP2), 105 or 66 (LP4) kbps"
                             : "mono input supports 66 or 52 kbps";
    size_t n = r.data_len / (2 * channels);
    if (warnings)
        *warnings = (n <= 6144 ? AT3_WARN_SHORT_INPUT : 0) | (r.has_smpl ? AT3_WARN_LOOP_IGNORED : 0);
    int fb = at3_encoder_frame_bytes(enc);
    size_t cap = AT3_HEADER_BYTES + ((n + 955) / AT3_FRAME_SAMPLES + 5) * (size_t)fb;
    uint8_t *buf = malloc(cap);
    int16_t *pcm = malloc(4096 * sizeof(int16_t));
    if (!buf || !pcm) {
        free(buf);
        free(pcm);
        at3_encoder_destroy(enc);
        return "out of memory";
    }
    size_t pos = AT3_HEADER_BYTES;
    for (size_t s = 0; s < n; s += 2048) {   /* little-endian samples, any alignment */
        size_t m = n - s < 2048 ? n - s : 2048;
        const uint8_t *p = r.data + s * 2 * channels;
        for (size_t i = 0; i < m * channels; i++)
            pcm[i] = (int16_t)rd16(p + 2 * i);
        pos += at3_encoder_write(enc, pcm, m, buf + pos);
    }
    pos += at3_encoder_finish(enc, buf + pos);
    at3_write_header(buf, channels, fb, at3_encoder_joint_stereo(enc), (uint32_t)n,
                     (uint32_t)(pos - AT3_HEADER_BYTES));
    at3_encoder_destroy(enc);
    free(pcm);
    *out = buf;
    *out_len = pos;
    return NULL;
}

static uint8_t *put16(uint8_t *p, unsigned v) { p[0] = v & 0xff; p[1] = (v >> 8) & 0xff; return p + 2; }
static uint8_t *put32(uint8_t *p, uint32_t v) { put16(p, v & 0xffff); put16(p + 2, v >> 16); return p + 4; }

const char *at3_decode_at3(const uint8_t *at3, size_t at3_len, uint8_t **out, size_t *out_len)
{
    riff r;
    const char *err = parse_riff(at3, at3_len, &r);
    if (err)
        return err;
    if (rd16(r.fmt) != 0x270)
        return "not an ATRAC3 file (format 0x270)";
    int channels = rd16(r.fmt + 2), align = rd16(r.fmt + 12);
    int js = r.fmt_len >= 26 ? rd16(r.fmt + 24) : 0;
    at3_decoder *dec = at3_decoder_create(channels, align, js);
    if (!dec)
        return "unsupported ATRAC3 parameters";
    uint64_t nframes = r.data_len / align;
    /* decoded sample i is source sample i - (K + 69); fact = {total, [K,] K} */
    uint64_t skip = 1024 + 69, total = nframes * 1024;
    if (r.fact && r.fact_len >= 8) {
        skip = rd32(r.fact + r.fact_len - 4) + (uint64_t)69;
        total = rd32(r.fact);
    }
    uint64_t avail = nframes * 1024 > skip ? nframes * 1024 - skip : 0;
    if (total > avail)
        total = avail;
    uint64_t data_bytes = total * 2 * channels;
    if (data_bytes > 0xffffffffu - 44 || data_bytes > SIZE_MAX - 44) {
        at3_decoder_destroy(dec);
        return "decoded audio would exceed 4 GB";
    }
    uint8_t *buf = malloc((size_t)(44 + data_bytes));
    int16_t *pcm = malloc(2048 * sizeof(int16_t));
    if (!buf || !pcm) {
        free(buf);
        free(pcm);
        at3_decoder_destroy(dec);
        return "out of memory";
    }
    uint8_t *p = buf;
    memcpy(p, "RIFF", 4);
    p = put32(p + 4, (uint32_t)(36 + data_bytes));
    memcpy(p, "WAVEfmt ", 8);
    p = put32(p + 8, 16);
    p = put16(p, 1);
    p = put16(p, channels);
    p = put32(p, 44100);
    p = put32(p, 44100 * 2 * channels);
    p = put16(p, 2 * channels);
    p = put16(p, 16);
    memcpy(p, "data", 4);
    p = put32(p + 4, (uint32_t)data_bytes);
    uint64_t at = 0;                         /* decoded sample index */
    for (size_t f = 0; f < nframes && at < skip + total; f++, at += 1024) {
        if (at3_decoder_decode_frame(dec, r.data + f * align, pcm)) {
            free(buf);
            free(pcm);
            at3_decoder_destroy(dec);
            return "malformed ATRAC3 frame";
        }
        for (int i = 0; i < 1024; i++) {
            uint64_t s = at + i;
            if (s < skip || s >= skip + total)
                continue;
            for (int c = 0; c < channels; c++)
                p = put16(p, (uint16_t)pcm[i * channels + c]);
        }
    }
    at3_decoder_destroy(dec);
    free(pcm);
    *out = buf;
    *out_len = (size_t)(44 + data_bytes);
    return NULL;
}
