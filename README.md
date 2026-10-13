# ATRAC3 (MiniDisc LP2/LP4) specification

A specification of Sony's ATRAC3 audio codec, as used by MiniDisc LP2 (132 kbps) and LP4 (66 kbps), written for the digital preservation of MiniDisc audio. It was reconstructed by observing Sony's reference encoder/decoder (`at3tool` 3.0.0.0) and verifying an independent implementation against it.

- **[spec/](spec/)** contains the specification: the bitstream, the decoder, both of Sony's encoders and all tables.
- **[spec/reference/](spec/reference/)** contains the Python reference models.
  - The decoder matches Sony's within ±1 LSB.
  - The LP2 and LP4 encoders produce output byte-identical to Sony's.
  - They need only NumPy.
- **[tools/](tools/)** contains the verification and analysis tooling: a CPU-emulator harness for the reference tool, comparison scripts, signal generators and a hardware test-kit generator.

- **[src/](src/)** is a C99 implementation of the decoder and both encoders, built as the `atrac3` command-line tool and as a WebAssembly module ([web/](web/)). Its encoder output is byte-identical to the reference models, and so to Sony's, and it runs several hundred times faster than real time. See [Source code](#source-code).

```bash
python spec/reference/atrac3_encoder_lp2.py in.wav out.at3        # LP2 (stereo) / 66 kbps mono
python spec/reference/atrac3_encoder_lp4.py in.wav out.at3 66     # LP4 joint stereo
python spec/reference/atrac3_decoder.py in.at3 out.wav
```

## Command-line tool

`atrac3` takes the same arguments as Sony's `at3tool`. It reads 16-bit 44.1 kHz WAV files.

```bash
atrac3 -e -br 132 in.wav out.at3    # LP2 (stereo input; the default)
atrac3 -e -br 66 in.wav out.at3     # LP4 joint stereo, or 66 kbps mono with mono input
atrac3 -e -br 105 in.wav out.at3    # 105 kbps stereo
atrac3 -e -br 52 in.wav out.at3     # 52 kbps mono
atrac3 -d in.at3 out.wav            # decode, sample-aligned with the source
```

Prebuilt binaries for Linux (x86_64 and arm64, statically linked), macOS (universal), Windows (x86_64) and WebAssembly are attached to each [release](https://github.com/kai-h/atrac-encoder/releases). The macOS binary isn't signed, so remove the quarantine flag once with `xattr -d com.apple.quarantine atrac3`.

To build it yourself you need CMake and a C compiler:

```bash
cmake -S . -B build
cmake --build build
python3 tests/run_tests.py          # byte-exact regression test, standard library only
```

The WebAssembly module needs [Emscripten](https://emscripten.org):

```bash
emcmake cmake -S . -B build-wasm
cmake --build build-wasm
python3 tests/run_tests.py --wasm build-wasm
```

`build-wasm` then holds `atrac3.mjs` (the JavaScript API, described under [WebAssembly module](#webassembly-module)), the compiled module and `index.html`, a small page that encodes and decodes files in the browser. Serve the folder over HTTP to use it.

**Floating-point arithmetic must not be fused or reordered.** Byte-exact output depends on every operation rounding exactly as written. The CMake build disables multiply-add contraction (`-ffp-contract=off`, or `/fp:precise` with MSVC), and the sources carry the equivalent pragmas for other build systems. Don't build with `-ffast-math`, and on 32-bit x86 use SSE2 rather than the x87 FPU. The test suite runs on every platform in CI and catches any deviation. `tools/verify_c.py` compares the C build directly with the Python models.

## Source code

The C implementation is plain C99 with no dependencies beyond the C standard library and libm. CMake builds it as one static library, `libatrac3`, which both the command-line tool and the WebAssembly module link against, so the two produce identical output.

### Library

The public API is [src/atrac3.h](src/atrac3.h), in three parts:

- **Streaming encoder**: `at3_encoder_create`, `at3_encoder_write` and `at3_encoder_finish` take interleaved 16-bit PCM in blocks of any size and return whole ATRAC3 frames, with `at3_write_header` writing the 80-byte `.at3` header. The framing matches `at3tool`'s, so decoded sample *i* corresponds to source sample *i* − 1093.
- **Frame decoder**: `at3_decoder_create` and `at3_decoder_decode_frame` decode one frame to 1024 samples per channel.
- **Whole files in memory**: `at3_encode_wav` and `at3_decode_at3` convert a WAV file image to an `.at3` file image and back. Both front ends use these.

Each encoder and decoder source file mirrors one Python module in [spec/reference/](spec/reference/), and its comments name the corresponding Python function, so the C code can be read side by side with the specification.

| File | Contents |
|---|---|
| `atrac3.h` | Public API |
| `internal.h` | Declarations shared by the library sources, and the floating-point pragmas |
| `encoder.c` | Mode selection, `at3tool`'s framing and the `.at3` header |
| `lp2_analysis.c` | LP2 encoder: QMF band split, gain control and MDCT |
| `lp2_alloc.c` | LP2 encoder: tonal component extraction, bit allocation and quantisation |
| `lp2_bitstream.c` | LP2 encoder: bit counting and frame packing |
| `lp4_analysis.c` | LP4 encoder: QMF band split, block envelope and MDCT |
| `lp4_gain.c` | LP4 encoder: gain control (transient pre-processing) |
| `lp4_stereo.c` | LP4 encoder: joint stereo |
| `lp4_alloc.c` | LP4 encoder: tonal component extraction, bit allocation and quantisation |
| `lp4_encoder.c` | LP4 encoder: frame packing |
| `decoder.c` | Frame decoder |
| `files.c` | WAV and `.at3` file reading and writing |
| `tables.c` | Codec tables, generated from the reference models by `tools/gen_c_tables.py` |

Sony uses two separate ATRAC3 encoders. The LP2 encoder also handles 66 kbps mono, and the LP4 encoder also handles 105 kbps stereo and 52 kbps mono.

### Command-line tool

[cli/atrac3.c](cli/atrac3.c) parses `at3tool`'s arguments, reads the input file into memory, converts it with `at3_encode_wav` or `at3_decode_at3` and writes the result. It prints a note when the input has 6144 samples or fewer (which `at3tool` refuses) or carries loop points, which are ignored.

### WebAssembly module

Emscripten compiles the same library, plus two small files from [web/](web/):

- [web/glue.c](web/glue.c) exports the whole-file functions to JavaScript, keeping each result in the module until the next call so JavaScript never handles pointers to pointers.
- [web/atrac3.mjs](web/atrac3.mjs) is the JavaScript API, an ES module that works in browsers, web workers and Node. It copies the input into the module's memory, calls the export and returns a copy of the result, throwing an `Error` with the library's message on failure.

```js
import { createAtrac3 } from './atrac3.mjs';
const atrac3 = await createAtrac3();
const at3 = atrac3.encodeWav(wavBytes, 132);   // Uint8Array .at3 file
const wav = atrac3.decodeAt3(at3);             // Uint8Array 16-bit WAV file
```

`encodeWav` takes the same bitrates as the command-line tool, and 0 selects 132 kbps for stereo and 66 kbps for mono. [web/index.html](web/index.html) is the demo page. [tests/wasm_cli.mjs](tests/wasm_cli.mjs) wraps the module in the command-line tool's arguments, which lets `tests/run_tests.py` run the same tests against both builds.

## Not included

This repository contains no Sony software or files. The numeric tables the codec needs (code tables, filter coefficients and the encoders' tuning constants) were read from Sony's reference tool, so that the output matches it. The verification tools expect the reference binary at `private/psp_at3tool.exe` (gitignored), which you must supply yourself. The specification and reference models stand on their own without it.

ATRAC1 (MiniDisc SP) and ATRAC3plus are not covered.

## Licence

The specification, reference models, C implementation and tools are released under the [MIT licence](LICENSE).

## Trademarks

ATRAC, MiniDisc, MDLP and PSP are trademarks of Sony Group Corporation or its affiliates. This project is not affiliated with or endorsed by Sony.
