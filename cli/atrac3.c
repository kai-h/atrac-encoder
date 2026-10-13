/* atrac3: command-line ATRAC3 encoder/decoder with at3tool's interface. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atrac3.h"

#ifndef AT3_VERSION
#define AT3_VERSION "dev"
#endif

static const char USAGE[] =
    "atrac3 " AT3_VERSION " - ATRAC3 (MiniDisc LP2/LP4) encoder, byte-identical to Sony at3tool\n"
    "\n"
    "usage: atrac3 -e [-br KBPS] in.wav out.at3     encode 16-bit 44.1 kHz WAV\n"
    "       atrac3 -d in.at3 out.wav                decode to 16-bit WAV\n"
    "\n"
    "  -br 132   LP2, stereo (default for stereo input)\n"
    "  -br 66    LP4 joint stereo with stereo input; 66 kbps mono with mono input (default)\n"
    "  -br 105   105 kbps stereo\n"
    "  -br 52    52 kbps mono\n";

static unsigned char *read_file(const char *path, size_t *len)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        return NULL;
    size_t cap = 1 << 20, n = 0;
    unsigned char *buf = malloc(cap);
    while (buf) {
        n += fread(buf + n, 1, cap - n, f);
        if (n < cap)
            break;
        unsigned char *nb = realloc(buf, cap *= 2);
        if (!nb)
            free(buf);
        buf = nb;
    }
    if (ferror(f)) {
        free(buf);
        buf = NULL;
    }
    fclose(f);
    *len = n;
    return buf;
}

static int write_file(const char *path, const unsigned char *data, size_t len)
{
    FILE *f = fopen(path, "wb");
    if (!f)
        return -1;
    size_t n = fwrite(data, 1, len, f);
    return (fclose(f) == 0 && n == len) ? 0 : -1;
}

int main(int argc, char **argv)
{
    int mode = 0, kbps = 0;
    const char *paths[2];
    int npaths = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "-d")) {
            mode = argv[i][1];
        } else if (!strcmp(argv[i], "-br") && i + 1 < argc) {
            kbps = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version")) {
            puts("atrac3 " AT3_VERSION);
            return 0;
        } else if (argv[i][0] == '-' && argv[i][1]) {
            fputs(USAGE, stderr);
            return 2;
        } else if (npaths < 2) {
            paths[npaths++] = argv[i];
        } else {
            fputs(USAGE, stderr);
            return 2;
        }
    }
    if (!mode || npaths != 2) {
        fputs(USAGE, stderr);
        return 2;
    }
    size_t in_len, out_len;
    unsigned char *in = read_file(paths[0], &in_len), *out = NULL;
    if (!in) {
        fprintf(stderr, "atrac3: cannot read %s\n", paths[0]);
        return 1;
    }
    const char *err;
    int warnings = 0;
    if (mode == 'e')
        err = at3_encode_wav(in, in_len, kbps, &out, &out_len, &warnings);
    else
        err = at3_decode_at3(in, in_len, &out, &out_len);
    free(in);
    if (err) {
        fprintf(stderr, "atrac3: %s: %s\n", paths[0], err);
        return 1;
    }
    if (warnings & AT3_WARN_SHORT_INPUT)
        fprintf(stderr, "atrac3: note: input has 6144 samples or fewer (at3tool refuses such input)\n");
    if (warnings & AT3_WARN_LOOP_IGNORED)
        fprintf(stderr, "atrac3: note: loop points (smpl chunk) are ignored\n");
    if (write_file(paths[1], out, out_len)) {
        fprintf(stderr, "atrac3: cannot write %s\n", paths[1]);
        at3_free(out);
        return 1;
    }
    at3_free(out);
    return 0;
}
