/* Emscripten exports for web/atrac3.mjs: whole-file encode/decode with the result kept
 * here until the next call, so JavaScript needs no pointer-to-pointer handling. */
#include <emscripten.h>
#include "atrac3.h"

static uint8_t *result;
static size_t result_len;

EMSCRIPTEN_KEEPALIVE const char *at3w_encode(const uint8_t *wav, size_t len, int kbps)
{
    at3_free(result);
    result = NULL;
    result_len = 0;
    return at3_encode_wav(wav, len, kbps, &result, &result_len, NULL);
}

EMSCRIPTEN_KEEPALIVE const char *at3w_decode(const uint8_t *at3, size_t len)
{
    at3_free(result);
    result = NULL;
    result_len = 0;
    return at3_decode_at3(at3, len, &result, &result_len);
}

EMSCRIPTEN_KEEPALIVE uint8_t *at3w_result(void) { return result; }
EMSCRIPTEN_KEEPALIVE size_t at3w_result_len(void) { return result_len; }
