/* ATRAC3 (MiniDisc LP2/LP4) encoder and decoder.
 *
 * C implementation of the reference models in spec/reference.  The encoder produces output
 * byte-identical to Sony's at3tool 3.0.0.0 (`at3tool -e -br N in.wav out.at3`).
 *
 * Supported modes (kbps, input channels):
 *   132, 2   MiniDisc LP2           (Sony "algorithm 1")
 *    66, 2   MiniDisc LP4, joint stereo (algorithm 0)
 *   105, 2   105 kbps stereo        (algorithm 0)
 *    66, 1   66 kbps mono           (algorithm 1)
 *    52, 1   52 kbps mono           (algorithm 0)
 * Input is 16-bit PCM at 44.1 kHz.
 */
#ifndef ATRAC3_H
#define ATRAC3_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AT3_HEADER_BYTES 80         /* RIFF + fmt (32) + fact (12) + data chunk header */

/* ---------------------------------------------------------------- streaming encoder */
typedef struct at3_encoder at3_encoder;

/* Returns NULL for an unsupported (kbps, channels) combination or when out of memory. */
at3_encoder *at3_encoder_create(int kbps, int channels);
void at3_encoder_destroy(at3_encoder *enc);
int at3_encoder_frame_bytes(const at3_encoder *enc);
int at3_encoder_joint_stereo(const at3_encoder *enc);
/* Output space needed by at3_encoder_write for nsamples samples per channel. */
size_t at3_encoder_max_output(const at3_encoder *enc, size_t nsamples);
/* Feed interleaved PCM (nsamples per channel); writes whole frames, returns bytes written.
 * Framing is at3tool's: the stream is primed with 955 zero samples and its first frame is
 * dropped, so decoded sample i corresponds to source sample i - 1093. */
size_t at3_encoder_write(at3_encoder *enc, const int16_t *pcm, size_t nsamples, uint8_t *out);
/* End of input: writes the remaining frames (at most 4 * frame_bytes), returns bytes. */
size_t at3_encoder_finish(at3_encoder *enc, uint8_t *out);

/* .at3 (RIFF/WAVE, format 0x270) header as written by at3tool, AT3_HEADER_BYTES long. */
void at3_write_header(uint8_t *out, int channels, int frame_bytes, int joint_stereo,
                      uint32_t total_samples, uint32_t data_bytes);

/* ---------------------------------------------------------------- frame decoder */
typedef struct at3_decoder at3_decoder;

at3_decoder *at3_decoder_create(int channels, int frame_bytes, int joint_stereo);
void at3_decoder_destroy(at3_decoder *dec);
/* Decode one frame into 1024 interleaved samples per channel.  Returns 0, or -1 for a
 * malformed frame. */
int at3_decoder_decode_frame(at3_decoder *dec, const uint8_t *frame, int16_t *pcm);

/* ---------------------------------------------------------------- whole files in memory */
#define AT3_WARN_SHORT_INPUT 1      /* <= 6144 samples: at3tool refuses such input */
#define AT3_WARN_LOOP_IGNORED 2     /* input has a smpl (loop) chunk; loops are not supported */

/* Encode a 16-bit 44.1 kHz mono/stereo WAV file image into an .at3 file image.  kbps 0
 * selects 132 for stereo and 66 for mono.  On success returns NULL and sets *out (free with
 * at3_free); on failure returns a static error message.  warnings may be NULL. */
const char *at3_encode_wav(const uint8_t *wav, size_t wav_len, int kbps,
                           uint8_t **out, size_t *out_len, int *warnings);
/* Decode an .at3 file image into a 16-bit WAV file image, sample-aligned with the source
 * (the codec delay is removed and the length taken from the fact chunk, like at3tool -d). */
const char *at3_decode_at3(const uint8_t *at3, size_t at3_len, uint8_t **out, size_t *out_len);
void at3_free(void *p);

#ifdef __cplusplus
}
#endif
#endif
