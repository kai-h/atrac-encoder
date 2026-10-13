# ATRAC3 codec specification

**Derived from:** Sony Computer Entertainment "SCEI ATRAC3plus Codec TOOL" version 3.0.0.0 (`psp_at3tool.exe`, built 18 August 2009, library version 1.02.00)
**Scope:** ATRAC3 at 44.1 kHz as used by MiniDisc LP2 (132 kbps) and LP4 (66 kbps), plus the 105 kbps stereo and 52/66 kbps mono modes the tool also supports
**Status:** complete for all five ATRAC3 modes of the tool. The decoder is verified against Sony's to within float rounding. Both Sony encoders are reproduced byte for byte by the reference models (§7, §8).

This document describes the codec's behaviour in our own words. It was written by analysing the reference binary and observing it under emulation, and by differentially testing our own implementation against it. Every statement marked **[verified]** has been checked numerically against the reference. Statements marked **[inferred]** come from code reading only, because the reference encoder never produced that case. Tables live in [Appendix A](ATRAC3-tables.md). Executable reference models live in [`reference/`](reference/).

---

## 1. Overview

ATRAC3 is a hybrid sub-band/transform codec. The signal chain for each channel is as follows.

```
PCM 44.1 kHz ──► 4-band QMF analysis tree ──► per band: gain control ─► 256-line MDCT (512-sample window, 50 % overlap)
              (bands 0–5.5, 5.5–11, 11–16.5, 16.5–22 kHz, each 11025 Hz sample rate)
                                               │
          1024 spectral lines per channel per frame (4 bands × 256)
                                               │
       tonal component extraction (optional) + 32 block-floating subbands
                                               │
       scale factor (6 bit) + quantiser selector (3 bit) per subband, mantissas VLC- or CLC-coded
```

A **frame** covers 1024 PCM samples per channel (23.22 ms) and has a fixed byte size. Each channel's data in a frame is a **sound unit**. In joint-stereo mode, two sound units (a sum signal and a difference signal) share one frame, and their channel weighting and matrixing parameters travel in the second unit.

The decoder is deterministic, and its output is fixed by this specification up to floating-point rounding. Sony's decoder itself is not bit-exact against an exact-arithmetic decoder; differences of ±1 LSB occur in about 0.2 % of samples. The encoder is not normative. This document describes Sony's two encoder implementations so that their behaviour (and quality) can be reproduced.

## 2. Operating modes

| Mode | kbps | Channels | Frame bytes | Sound units | Joint stereo | Sony encoder used | Codec-info word |
|---|---|---|---|---|---|---|---|
| MD **LP2** | 132 | 2 | 384 | 2 × 192 (independent) | no | algorithm 1 | `0x002030` |
| PSP 105 | 105 | 2 | 304 | 2 × 152 (independent) | no | algorithm 0 | `0x002026` |
| MD **LP4** | 66 | 2 | 192 | 1 frame shared by 2 units | **yes** | algorithm 0 | `0x022018` |
| 66 mono | 66 | 1 | 192 | 1 × 192 | no | algorithm 1 | `0x002030` |
| 52 mono | 52 | 1 | 152 | 1 × 152 | no | algorithm 0 | `0x002026` |

Bit rate = frame_bytes × 8 × 44100 / 1024. The table in the binary sits at `0x451A50`, with 9 dwords per entry: codec family, kbps, channels, samples per frame, sample rate, block align, codec-info word, encoder algorithm and the "mono via stereo" flag. **[verified]**

**Mono modes** are implemented as half of a stereo independent-channel frame. **[verified]**
- **Encoder:** it duplicates the mono input into both channels, encodes with the stereo configuration (304 or 384-byte frames), and keeps only the first sound unit (152 or 192 bytes) of each frame.
- **Decoder:** it copies the single unit into both halves of a stereo frame, decodes, and outputs channel 0. A mono stream is therefore exactly one independent-mode sound unit per frame.

### 2.1 Codec-info word

This is a 24-bit word in the same layout as Sony's OMA/OpenMG headers. Bits 24–31 select the codec family (0 = ATRAC3, 1 = ATRAC3plus, 0x20–0x22 = ATRAC Advanced Lossless). **[verified]**

| bits | meaning (ATRAC3) |
|---|---|
| 0–9 | frame size in bytes / 8 (stereo frame) |
| 13–15 | sample rate: 0 = 32000, **1 = 44100**, 2 = 48000 (the Sony encoder and decoder only accept 44100) |
| 17 | joint stereo (1) or independent channels (0) |

## 3. File container (`.at3`, RIFF/WAVE)

The tool reads and writes little-endian RIFF/WAVE files. The chunk order on output is `fmt `, `fact`, optional `smpl`, then `data`. **[verified]**

### 3.1 `fmt ` chunk (32 bytes)

| offset | size | value |
|---|---|---|
| 0 | 2 | `wFormatTag` = 0x0270 (WAVE_FORMAT_SONY_SCX / ATRAC3) |
| 2 | 2 | channels (1 or 2) |
| 4 | 4 | 44100 |
| 8 | 4 | `nAvgBytesPerSec` = round(44100 × block_align / 1024 + 0.5), i.e. 16538 / 13092 / 8269 for 384 / 304 / 192 bytes |
| 12 | 2 | `nBlockAlign` = frame bytes (384, 304, 192, 152) |
| 14 | 2 | `wBitsPerSample` = 0 |
| 16 | 2 | `cbSize` = 14 |
| 18 | 2 | 1 |
| 20 | 4 | channels × 2048 (0x1000 for stereo, 0x800 for mono) |
| 24 | 2 | joint-stereo flag (codec-info bit 17) |
| 26 | 2 | joint-stereo flag (repeated) |
| 28 | 2 | 1 |
| 30 | 2 | 0 |

The reader requires the word at offset 18 to be 1, the two flags at offsets 24 and 26 to be equal, and the words at offsets 28 and 30 to be 1 and 0 respectively. **[verified, from the reader code]**

### 3.2 `fact` chunk (12 bytes)

`{ total_samples, K, K }`. Here `total_samples` is the number of PCM samples per channel in the source, and K is the **decoder skip value**. Without loop points K = 1024. With loop points see §9.2.

Files written by at3tool versions before 2.2.0 use an 8-byte `fact` chunk, `{ total_samples, K }`; several of the PSP SDK's sample files are in this form. A reader must accept both, taking K from the last dword. For ATRAC3 the meaning of K is the same in both forms. (For ATRAC3plus the old form needs 184 added; that is out of scope here.) **[from the reader code; the SDK sample files use both forms]**

### 3.3 `smpl` chunk (60 bytes, only when loop points are given)

This is the standard sampler chunk with a single loop:
- `dwSamplePeriod` = round(10⁹/44100 + 0.5) = 22676,
- `dwMIDIUnityNote` = 60,
- `cSampleLoops` = 1,
- `cbSamplerData` = 24 (a Sony quirk: this holds the size of one loop record),
- the loop record is `{0, 0, loop_start + K, loop_end + K, 0, 0}`.

### 3.4 `data` chunk

This holds the concatenated frames, a whole number of `nBlockAlign`-byte frames.

### 3.5 Timing

Decoded sample *i* (counting from the first sample of the first frame in the file) corresponds to source sample *i − (K + 69)*. A decoder that wants the original timing therefore discards the first `K + 69` output samples and emits `total_samples` samples. **[verified: the Sony decoder's output is sample-aligned with the source]**

> FFmpeg's `atrac3` decoder ignores this delay (it outputs from sample 0) and its output is also **inverted in polarity** relative to Sony's. Do not use it as a timing or sign reference.

## 4. Bitstream syntax

All fields are unsigned and read MSB-first. Bit position 0 is the MSB of the first byte of the sound unit.

### 4.1 Frame layout

**Independent channels (LP2, 105 kbps, mono).** The frame splits into `channels` equal parts, and part *c* (frame_bytes/channels bytes) holds the sound unit of channel *c*. Unused trailing bits in a unit are padding. (Sony's decoder ignores them, and Sony's encoders write zero bits.)

**Joint stereo (LP4).** The frame holds two sound units:
- **Unit 0** (the sum or "M" signal) starts at the first bit of the frame and is read forwards.
- **Unit 1** is stored **byte-reversed from the end of the frame**: its first byte is the frame's last byte, its second byte is the frame's second-to-last byte, and so on. Unit 1 starts with a 16-bit joint-stereo header (§4.6).
- **Gap bytes:** Sony's encoders fill the gap between the end of unit 0 and the (reversed) end of unit 1 with zero bytes. **[verified]**
- **Decoders:** after reversing, a decoder must skip any leading `0xF8` bytes before reading unit 1. Sony's decoder does this, presumably for streams from other encoders. It also flags a frame as damaged when unit 0's bits run into the area occupied by unit 1. **[inferred, from decoder code]**

### 4.2 Sound unit

```
sound_unit(js_second):
    if js_second:  js_header()            // 12 bits, see §4.6
                   unit_id      : 2  == 3
    else:          unit_id      : 6  == 0x28 (binary 101000)
    num_qmf_bands_m1 : 2                  // QMF bands carrying data: 1..4
    for b in 0 .. num_qmf_bands_m1:  gain_control(b)
    tonal_components(num_qmf_bands_m1)
    spectrum()
```

QMF bands above `num_qmf_bands_m1` have no gain points. Tonal components can't be placed in them either. (Sony's decoder uses the count only for those two things; spectral data in `spectrum()` may still extend above it.)

### 4.3 Gain control data

```
gain_control(b):
    num_points : 3                        // 0..7
    repeat num_points:
        level    : 4                      // gain 2^(4 - level)
        location : 5                      // sample 8·location in the 256-sample band block
```

Locations must be strictly increasing within a band. Sony's decoder treats a violation as a frame error. **[verified from decoder]**

### 4.4 Tonal components

```
tonal_components(nb_m1):
    num_groups : 5
    if num_groups == 0: return
    mode_select : 2                       // 0: all VLC, 1: all CLC, 3: per group; 2 is invalid
    repeat num_groups:
        band_flag[0..nb_m1] : 1 each      // which QMF bands carry components of this group
        values_m1   : 3                   // values per component = values_m1 + 1 (1..8)
        selector    : 3                   // quantiser selector, must be >= 2
        coding_mode : (mode_select == 3) ? 1 bit : (mode_select)
        for b in 0..nb_m1 where band_flag[b]:
            for g in 0..3:                // four 64-line groups per QMF band
                count : 3
                repeat count:
                    sf_index : 6
                    position : 6          // first line = 256·b + 64·g + position
                    n = min(values_m1 + 1, 1024 - first_line)
                    mantissas(selector, coding_mode, n)
```

At most 64 components may appear per sound unit (decoder error otherwise). Component values are **added** to the spectrum at lines `first_line .. first_line + n − 1`. They use the same dequantisation as §4.5. **[verified]**

### 4.5 Spectrum (block floating subbands)

```
spectrum():
    num_subbands_m1 : 5                   // subbands 0 .. num_subbands_m1 are described
    coding_mode     : 1                   // 0 = VLC, 1 = CLC
    selector[i]     : 3   for i in 0..num_subbands_m1
    sf_index[i]     : 6   for each i with selector[i] != 0
    for i in 0..num_subbands_m1 with selector[i] != 0:
        mantissas(selector[i], coding_mode, width of subband i)   // Appendix A.1
```

Lines not covered by a coded subband are zero (before tonal components are added).

**Mantissa coding** (Appendix A.4):
- **Selector 1:** each codeword carries two consecutive mantissas (a subband width is always even).
- **Selectors 2–7:** each codeword carries one mantissa.
- **CLC codes:** a fixed-width two's-complement field per value, `2·2`, 3, 3, 4, 4, 5 and 6 bits for selectors 1–7. Codewords that represent an out-of-range value are invalid.
- **VLC codes:** these are Huffman codes. For selectors 6 and 7 the maximum-magnitude values (±15 and ±31) get short 4-bit codes, because the scale factor is chosen from each block's peak.

**Dequantisation** (the same formula for subbands and tonal components): **[verified]**

```
coef = mantissa × SF[sf_index] / MAX_QUANT[selector]
SF[i] = 2^((i − 15)/3),   MAX_QUANT = {—, 1.5, 2.5, 3.5, 4.5, 7.5, 15.5, 31.5}
```

The resulting coefficients are already in the scale used by §5; no other gain constant is involved.

### 4.6 Joint-stereo header (start of unit 1)

```
js_header():
    weight_flag   : 1
    weight_level  : 3                     // weight index W = 8·weight_flag + weight_level
    matrix_sel[b] : 2   for b = 0..3      // per QMF band; value 2 is invalid
```

Then the `unit_id` (2 bits, value 3) follows. In byte terms, the first two bytes of the reversed unit are: high nibble of byte 0 = W, then four 2-bit matrix selectors, then `11`, then `num_qmf_bands_m1`.

## 5. Decoding process

For each frame, each sound unit is parsed (§4), synthesised into four 256-sample QMF-band signals (§5.1–5.3), joint-stereo-processed if applicable (§5.4), and recombined to 1024 PCM samples (§5.5–5.6). All arithmetic can be done in single precision; the reference model uses double.

### 5.1 Inverse MDCT per QMF band **[verified]**

For QMF band b (0..3), take the coefficients X[k] = spectrum[256·b + k], k = 0..255. **For odd bands (1 and 3) reverse the order first:** X[k] ← X[255 − k]. (The QMF analysis spectrally inverts odd bands.) Then:

```
y[n] = w[n] · Σ_{k=0}^{255} X[k] · cos( π/256 · (n + 128.5) · (k + 0.5) ),   n = 0..511
```

The transform is unnormalised and w is the window in Appendix A.5. A band whose spectrum is entirely zero yields y = 0. Sony's decoder skips the IMDCT for bands above the highest band with non-zero data in the current or previous frame, which is equivalent.

### 5.2 Gain compensation and overlap-add **[verified]**

Each QMF band keeps the second half of the previous frame's windowed IMDCT (`overlap[256]`) and the previous frame's gain points (`prev_points`). With `cur_points` being this frame's gain points for the band:

```
scale = (cur_points not empty) ? 2^(4 − cur_points[0].level) : 1
g[n]  = gain_curve(prev_points)[n]                      // n = 0..255, defined below
out[n] = (y[n] · scale + overlap[n]) · g[n]             // n = 0..255
overlap = y[256..511];  prev_points = cur_points
```

`gain_curve(P)` is piecewise constant with 8-sample exponential ramps. Start with `pos = 0` and append a terminal level of 4 (unity) after the last point:

```
for i, (L_i, loc_i) in P:
    start = 8 · loc_i
    g[pos .. start−1]   = 2^(4 − L_i)
    g[start + k]        = 2^(4 − L_i) · 2^(k · (L_i − L_{i+1}) / 8),  k = 0..7   (L_{last+1} = 4)
    pos = start + 8
g[pos .. 255] = 1
```

Gain points therefore describe the **previous** frame's second half. The encoder attenuated a loud part of the signal before the MDCT; the decoder restores it here.

### 5.3 Uncoded bands

QMF bands with no spectral data (all coefficients zero) still go through §5.2. Their output is the overlap from the previous frame times the gain curve, or zero if that was also silent.

### 5.4 Joint-stereo reconstruction (LP4) **[verified]**

This operates on the four QMF-band time signals of both units, after §5.2 and before §5.5. Let L_b[n] be unit 0's band-b signal and R_b[n] unit 1's.

The parameters are applied **one frame late**. When frame *t* is decoded, the parameters parsed from frame *t − 1* are "current" and those from frame *t − 2* are "previous". The initial state is matrix selector 3 for all bands, no weighting, and "last coded band" = 4.

The per-band coefficients are:
- `(c0, c1) = MATRIX[matrix_sel_cur[b]]`, with MATRIX = {0: (0, 1), 1: (1, 1), 3: (0.5, 0.5)}.
- `(w0, w1) = WEIGHTS[W_cur]` if `b > last_coded_band_cur`, else `(2, 2)`. Here `last_coded_band_cur` is the `num_qmf_bands_m1` of unit 1 in that same frame, and WEIGHTS is in Appendix A.8. The same is done with the "previous" parameters to obtain `(p0, p1)` and `(pw0, pw1)`.

Then, for each band b:

```
interp = (pw0 != w0) or (matrix_sel_prev[b] != matrix_sel_cur[b])
for n in 0..255:
    if interp and n < 8:
        t  = (n − 8) / 8
        a0 = c0 + t·(c0 − p0);   a1 = c1 + t·(c1 − p1)
        v0 = w0 + t·(w0 − pw0);  v1 = w1 + t·(w1 − pw1)
    else:
        a0, a1, v0, v1 = c0, c1, w0, w1
    s      = a0·L_b[n] + a1·R_b[n]
    R_b[n] = (L_b[n] − s) · v1
    L_b[n] = s · v0
```

With matrix 3 and no weighting this gives L = M + S and R = M − S, where unit 0 = M = (L+R)/2 and unit 1 = S = (L−R)/2. Weighting applies only to bands that unit 1 does not code, which gives intensity-stereo-style panning of the sum signal.

> FFmpeg's decoder differs here in two ways. It applies weighting to bands 1–3 regardless of how many bands unit 1 codes. Its non-interpolated matrix-1 formula is also inconsistent with its interpolated one. This matters in practice, because Sony's LP4 encoder uses matrix selector 1 whenever channel 0 dominates a band (§8). Sony's decoder, and the formula above, give `R = −2·S` in steady state. **[verified on LP4 streams using selector 1]**

### 5.5 QMF synthesis **[verified]**

This is a three-stage tree of identical two-band synthesis filters with the 48-tap prototype h (Appendix A.6, used as 2·h):

```
low  = qmf_synth(lo = band0, hi = band1)      // 256+256 -> 512 samples
high = qmf_synth(lo = band3, hi = band2)      // note: band 3 is the "low" input
pcm  = qmf_synth(lo = low,   hi = high)       // 512+512 -> 1024 samples

qmf_synth(lo[N], hi[N]) with persistent delay d[46]:
    t[0..45] = d;  for i in 0..N−1: t[46+2i] = lo[i] + hi[i];  t[47+2i] = lo[i] − hi[i]
    for j in 0..N−1:
        out[2j]   = Σ_{k=0}^{23} t[2j + 2k + 1] · 2·h[2k + 1]
        out[2j+1] = Σ_{k=0}^{23} t[2j + 2k]     · 2·h[2k]
    d = t[2N .. 2N+45]
```

Each of the three filter instances keeps its own delay line, per channel.

### 5.6 Output conversion **[verified]**

`sample = round_half_even(clamp(x, −32767, +32767))`. Note the symmetric clamp at −32767. Sony's decoder never produces −32768.

## 6. Encoder: common framework

Sony's two encoders are separate implementations. The first is written in portable C with x87 floating point (algorithm 1). The second uses SSE kernels (algorithm 0). Their decision logic differs completely (§7, §8), but they share the signal model below. Any ATRAC3 encoder must produce spectra on the scale the decoder expects (§5), and this section states that scale. **[verified in both chapters]**

### 6.1 Analysis filterbank

The encoder uses the time-reversed counterpart of §5.5, with the same 48-tap prototype h (Appendix A.6). One two-band split of a signal s (with 46 samples of history) is:

```
for j = 0, 1, 2, ...:
    P = Σ_{k odd}  h[k] · s[2j + k]
    Q = Σ_{k even} h[k] · s[2j + k]
    low[j]  = Q + P
    high[j] = P − Q
```

The tree is `(L, H) = split(x)`, then `(band0, band1) = split(L)` and `(band3, band2) = split(H)`. Note the swapped outputs in the high half, which mirror the decoder.

### 6.2 Forward MDCT and window

Each QMF band is transformed with 50 % overlap. The 512-sample input is the previous 256-sample block followed by the current one. With the analysis window `wa[n] = sin²(π(n + 0.5)/512)`:

```
X[k] = (1/128) · Σ_{n=0}^{511} wa[n] · x[n] · cos( π/256 · (n + 128.5) · (k + 0.5) ),   k = 0..255
```

For odd QMF bands (1 and 3) the 256 outputs are stored in reverse order. With 16-bit input in integer units, this gives exactly the coefficient scale of §4.5 and §5.1.

The decoder's synthesis window (Appendix A.5) equals `ws[n] = wa[n] / (wa[n]² + wa[(n + 256) mod 512]²)`; Appendix A.5 writes the same thing with `wi = 2·wa`. Because `sin⁴ + cos⁴ ≠ 1`, the window is not power-complementary, and this normalisation makes analysis followed by synthesis a perfect-reconstruction pair. **[verified numerically: reconstruction error 2·10⁻¹⁰ on random input]** Algorithm 1 computes the transform with an FFT. Algorithm 0 carries internal scale factors in its QMF (about 1.17·10⁹) and its MDCT window that cancel out. Both reach the formula above to within float32 rounding.

### 6.3 Gain control in the encoder

Gain points (§4.3) describe level changes within a 256-sample band block. Before the MDCT, the encoder divides each sample by the gain curve that the decoder will multiply back in (§5.2). The block being transformed is the one whose second half the points describe. A loud attack late in a block is therefore attenuated to the level of the quiet part before transformation, which confines pre-echo to the 8-sample transition. Each chapter gives its own attack and release detector.

### 6.4 Framing, priming and flushing (both algorithms)

- **Library calls:** each call passes up to 1024 samples per channel; a short block is zero-padded. Both libraries return **no frame for their first call**. A short block marks the end of the input. Further empty "flush" calls then produce frames until the library reports completion: 2 more frames if the short block held fewer than 886 samples, otherwise 3.
- **at3tool:** the input stream is 955 zero samples followed by the source (generally `Z = (K − 69) mod 1024`, §9.2). The tool **discards the first frame** the library returns whenever `K − 69 < 1024`, which is always the case without loop points (K = 1024). For N source samples let `T = N + 955` and `r = T mod 1024`. The file then holds `⌈T/1024⌉` frames if `0 < r < 886`, and one frame more otherwise.
- **End-to-end timing:** the file timing is that of §3.5 (decoded sample `i` ↔ source sample `i − (K + 69)`) for both encoders. Sony's own decoder, and the reference decoder here, reproduce the source sample-aligned. **[verified]**

## 7. Encoder algorithm 1 (LP2, 132 kbps; 66 kbps mono)

**Complete: see [ATRAC3-encoder-LP2.md](ATRAC3-encoder-LP2.md).** The reference model [`reference/atrac3_encoder_lp2.py`](reference/atrac3_encoder_lp2.py) produces files **byte-identical** to `at3tool -e -br 132`, and to `-br 66` with mono input. It was tested on 9,780 frames: synthetic signals, edge cases and PSP SDK music. **[verified]**

In brief: each channel is coded independently into a 192-byte unit. The pipeline is:
- QMF analysis.
- Gain control: attack/release detection on the current block plus one block of look-ahead, at most 7 points per band.
- A gain-windowed MDCT, with spectrum A coded and spectrum B (no gain window) used for decisions.
- Tonal extraction: either around one dominant isolated peak and its QMF alias, or wherever the spectrum lies more than 32 scale-factor steps (64 dB) above its mean.
- Initial word lengths from the scale-factor profile with a simple masking slope.
- A 15-step global offset search to meet the 1536-bit budget, dead-zone widening and word-length reduction if over budget, then five refinement passes that spend leftover bits.
- A final scale-factor correction.

Bandwidth is fixed at 15.1 kHz (29 subbands).

## 8. Encoder algorithm 0 (LP4, 66 kbps joint stereo; 105 kbps; 52 kbps mono)

**Complete: see [ATRAC3-encoder-LP4.md](ATRAC3-encoder-LP4.md).** The reference model [`reference/atrac3_encoder_lp4.py`](reference/atrac3_encoder_lp4.py) produces files **byte-identical** to `at3tool -e -br 66` and `-br 105`. Its first unit also reproduces `-br 52` mono. Testing covered 17 synthetic signals, 32,452 crafted sound units injected directly into the allocator, and over 4 minutes of PSP SDK music. **[verified]**

In brief, at 66 kbps the frame is split unevenly: unit 0 carries the sum signal (nominally up to 12.4 kHz, which a loud frame can extend), and unit 1 carries a side signal only up to 2.7 kHz.
- **Joint stereo, per QMF band:**
  - Band 0 is matrixed. Its matrix selector comes from the band's energy balance, with hysteresis.
  - The side signal is additionally divided by an untransmitted gain r, between 1 and about 1.68 and changing by at most 0.05 per frame. This narrows the stereo image when one channel dominates.
  - Bands 1–3 carry only the sum, with one intensity weight per frame.
- **Gain control:** a backwards decay search, then a forwards attack search, up to 7 points.
- **Bit split:** a peak-weighted comparison of the two spectra sets unit 0's budget. Unit 1 gets the rest of the 192 bytes.
- **Allocation:** a scale-factor survey, an all-tonal path for quiet frames, priority-driven word lengths, and tone extraction with VLC/CLC choice and tone splitting. After exact bit counting there is a final fit in order of scale factor.

## 9. Reference tool behaviour (`at3tool`)

### 9.1 Command line

```
at3tool -e [-br N] [-loop S E | -wholeloop] in.wav out.at3     encode (16-bit PCM WAV, 44.1 kHz)
at3tool -d [-repeat N] in.at3 out.wav                          decode to 16-bit PCM WAV
```

- **`-br`:** defaults to 128 kbps, which is an ATRAC3plus rate. For ATRAC3 use `-br 132`, `105` or `66` with stereo input, or `66` or `52` with mono input. ATRAC3plus rates are 32–352 kbps.
- **Supported combinations:** these are listed in §2. The bit rate together with the input channel count selects the table entry.
- **Input length:** input must be at least 6144 samples (0x1800).
- **Loops:** must satisfy `0 ≤ S < S + 6143 ≤ E < N`. A `smpl` chunk in the input is used if no `-loop` is given.
- **`-repeat N`:** the decoder plays the loop section N times (default 2).

### 9.2 Loop handling (PSP-specific, not used for MiniDisc)

With loop end E, the skip value is `K = 2047 − ((1024 + 69 + E) mod 1024)`. This puts sample E+1 of the decoded stream on a frame boundary. The encoder:
- prepends `Z = (K − 69) mod 1024` zeros,
- keeps the first frame the library returns when `K − 69 ≥ 1024`, and discards it otherwise (§6.4),
- after reaching E+1, rewinds the input to the loop start and encodes 2 further frames (8 for ATRAC3plus) so that a looping decoder has seamless overlap data.

**[verified for E = 100000 and whole-file loops]**

## 10. Verification

| Component | Method | Result |
|---|---|---|
| Bitstream parser | Spectra of every sound unit compared with the reference decoder's internal buffers, 12 streams (LP2 and LP4, 6 signal types) | identical to float32 precision |
| Synthesis (IMDCT, window, gain, overlap, JS, QMF) | Per-frame PCM before integer conversion | residual ≤ −97 dB in every frame (float rounding) |
| PCM output and delay | 16-bit WAV output of the reference CLI, all five modes (LP2, LP4, 105, 66 mono, 52 mono) | identical sample count and alignment; ±1 LSB in ~0.2 % of samples |
| Real music | 7 Sony-encoded ATRAC3 files shipped with the PSP SDK (5 × 132 kbps stereo, 2 × 52 kbps mono; ~23 000 frames, about 10 minutes, produced by several at3tool versions), frame-by-frame against the reference decoder | worst frame −90 dB, typically −110 dB |
| Joint-stereo matrix selector 1 | LP4 streams of one-sided and panned stereo (303 band-frames with selector 1) | residual ≤ −108 dB |
| Encoders | See §7 and §8: byte-identical output files | 100 % of frames |
| Not exercised by the Sony encoders (code reading only) | CLC spectrum coding, tonal `mode_select = 3`, `0xF8` padding | [inferred] |

Reference models: `reference/atrac3_decoder.py` (decoder), `reference/atrac3_encoder_lp2.py` and `reference/atrac3_encoder_lp4.py` (encoders) and `reference/at3_tables.py` (tables). Verification harnesses: `tools/verify_decoder.py`, `tools/verify_lp2.py` and `tools/verify_lp4.py`. These need the reference binary and the emulator in `tools/at3emu.py`.
