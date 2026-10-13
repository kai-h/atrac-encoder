# ATRAC3 encoder "algorithm 0" (MiniDisc LP4 66 kbps JS, 105 kbps; 52 kbps mono)

**Chapter 8 of the [ATRAC3 specification](ATRAC3.md).** Written from a byte-exact reverse-engineered model. Every section is
**verified by test** (the executable model [`reference/atrac3_encoder_lp4.py`](reference/atrac3_encoder_lp4.py) reproduces the Sony
encoder byte for byte) unless it is marked *inferred* or *unreachable*.

Executable companions:
- [`reference/atrac3_encoder_lp4.py`](reference/atrac3_encoder_lp4.py) – top level (`Encoder`, `encode_wav`),
  `lp4_analysis.py` (QMF, block peaks, MDCT), `lp4_stereo.py` (joint stereo, bit split),
  `lp4_gain.py` (gain control), `lp4_alloc.py` (allocation, tones, quantisation, bit counting),
  `lp4_bitstream.py` (packing), `lp4_fp.py` (arithmetic helpers),
  `lp4_tables.py` (all constant tables, generated from the binary by `tools/extract_lp4_tables.py`).
- `tools/verify_lp4.py` – comparison against the emulated Sony encoder.

Names in brackets, e.g. [allocate], are the corresponding functions and tables of the reference
model; addresses refer to `psp_at3tool.exe` and are for traceability only.

The bitstream syntax is the one in [ATRAC3.md](ATRAC3.md) §4 (decoded by [`reference/atrac3_decoder.py`](reference/atrac3_decoder.py)); this document only
describes how the encoder chooses what to write.

---------------------------------------------------------------------------------------------------

## 1. Overview

The encoder processes 1024 stereo samples per call and produces one frame (192 bytes at 66 kbps,
304 bytes at 105 kbps) containing two sound units.

```
PCM int16 L/R ──► QMF 4 bands ──► [joint stereo: matrix / intensity weights] ──► block peaks
                     (per ch)        (66 kbps only)                                 │
                                                                                    ▼
             frame ◄── packing ◄── allocation + tones ◄── [bit split] ◄── MDCT ◄── gain control
```

Per channel ("unit") and frame:

1. **QMF analysis** splits the 1024 new samples into four bands of 256 samples (0–5.5, 5.5–11,
   11–16.5, 16.5–22 kHz).
2. **Joint stereo** (66 kbps only, [JointStereo.analyse]): per band a decision between "left dominant",
   "right dominant" and "comparable"; band 0 is converted into a sum and a scaled difference signal,
   bands 1–3 carry only the sum, scaled by a weight that the decoder uses to rebuild L and R.
3. **Block peaks** of the new samples (32 blocks of 8 per band).
4. **Gain control**: using the peaks of the previous and the new samples, gain points are chosen for
   the *previous* frame's band samples, which are then multiplied by the gain curve.
5. **MDCT** of the previous frame's samples (second half of the window) together with the folded
   frame before that (first half, kept as state).
6. **Bit split** between the two units (66 kbps only, [split_bits]).
7. **Bit allocation** [allocate]: bandwidth, quantiser selectors, scale factors, tonal components.
8. **Packing**; at 66 kbps unit 2 is written backwards from the end of the frame.

Timing: the spectrum coded in a frame describes the band samples that arrived one call earlier
(the new samples are only look-ahead for gain control and for the stereo analysis). Joint-stereo
decisions and weights are transmitted one frame after they are made, which is exactly when the
decoder needs them. The tool's lead-in of 955 zero samples and its dropped first frame (section 3)
are related to these delays (*inferred*); measured end to end, `spec/reference/atrac3_decoder.py`
applied to a tool-encoded file reproduces the input delayed by 1093 = 1024 + 69 samples, at the
original level.

## 2. Configuration

The mode table at 0x48cf74 [IS_JOINT, CHANNEL_MODE, BITRATE, BLOCK_ALIGN,
CHANNEL_BYTES0/1, CHANNEL_HERTZ0/1] has six rows; the command line tool reaches row 0
(66 kbps) and row 2 (105 kbps) of algorithm 0:

| | 66 kbps (LP4) | 105 kbps |
|---|---|---|
| frame bytes | 192 | 304 |
| joint stereo | yes | no |
| nominal bytes unit 0 / unit 1 | 144 / 48 | 152 / 152 |
| bandwidth unit 0 / unit 1 | 12400 / 2700 Hz | 13700 / 13700 Hz |
| bytes per channel `bpc` (frame bytes / 2) | 96 | 152 |

Per unit c the encoder derives (0x434090, 0x434580 [Encoder.__init__]):

- static bit budget `budget_c = 8·bytes_c − 19`, minus a further 8 for the joint-stereo unit
  (1133 and 357 at 66 kbps; 1197 and 1197 at 105 kbps). 19 bits are the fixed fields of a unit
  (id, band count, tone group count, subband count, coding mode); the extra byte is the
  joint-stereo header. The 66 kbps unit-1 value is replaced every frame (section 9).
- bandwidth index `idx = floor(hz · 2048 / 44100)` and the subband limit `nsub_max` = 1 + the
  smallest `i` (0 ≤ i < 31) with `idx ≤ SB_START[i+1]`: 26 / 13 at 66 kbps, 27 / 27 at 105 kbps.
- for the joint-stereo unit: the number of *matrixed* QMF bands
  `limit = max(1, floor((idx + 255) / 256))` = 1 at 66 kbps. Without joint stereo, limit = −1 and the
  whole stereo stage is skipped.
- the joint-stereo unit is flagged; Sony's flag field holds unit 0's static budget (1133), which
  is reused as the base value of the bit split.

**52 kbps mono** is this encoder in the 105 kbps configuration, fed with the mono signal in both channels. The tool writes only the first 152-byte unit of each frame. **[verified: 131/131 frames identical to `at3tool -e -br 52`, mono input]**

Initial joint-stereo state: current decisions −3 for all bands, transmitted selectors 0, current
weight 15, transmitted weight 0, side gains 1.0 (current and previous), energies 0. All other state
(QMF, MDCT, gain envelopes, previous samples) starts at zero.

## 3. Call protocol, framing and container

### 3.1 Library calls (a3e0_encode / a3e0_flush, [Encoder.encode] / [Encoder.flush])

Each call passes up to 4096 bytes of interleaved little-endian int16 stereo PCM. A counter `cnt`
(initially 0) controls output:

```
on every call (flush calls pass 0 bytes):
    cnt = cnt + 1
    if bytes != 4096 and cnt > 0:            # first short block = end of input
        cnt = -3 if bytes < 3544 else -4      # 3544 bytes = 886 samples
    if cnt == 0: return "finished, no frame"
    convert the samples to float (channel 0 = even samples), zero-pad to 1024, encode one frame
    return no frame if cnt == 1 else the frame
```

So the first call primes the encoder and returns nothing; a short final block is followed by two
(fewer than 886 samples) or three further frames from flush calls; if the input ends exactly on a
block boundary the first flush call behaves like an empty short block (three flush frames).

### 3.2 at3tool framing (cli_encode_loop 0x403290)

`at3tool -e in.wav out.at3 -br 66|105` builds the first block from 955 zero samples followed by the
first 69 input samples (1024 minus the codec delay of 69), then feeds blocks of 1024 samples, the
last one short, and flushes. Because its lead-in is shorter than one block the tool discards the
first frame the library returns. WAV files carrying a `smpl` loop take a different path of the tool
(loop-aligned extra blocks); that path is not described here (*not modelled*).

Container (verified byte-exact for files without `smpl`): `RIFF` / `WAVE`, a 32-byte `fmt ` chunk
{tag 0x0270, channels 2, 44100 Hz, bytes/s = floor(44100·frame/1024 + 0.5), block align = frame bytes,
bits 0, cbSize 14, then 16-bit 1, 32-bit 0x1000, 16-bit js, 16-bit js, 16-bit 1, 16-bit 0}
(js = 1 at 66 kbps, 0 at 105 kbps), a 12-byte `fact` chunk {number of input samples, 1024, 1024},
and the `data` chunk with the frames.

### 3.3 Frame layout

At 105 kbps unit 0 occupies bytes 0..151 and unit 1 bytes 152..303, each starting at its first byte
(bits are written MSB first; unused bits and bytes are zero). At 66 kbps unit 0 is written from
byte 0; if it ends at bit position `e`, it occupies `n0 = ceil(e/8)` bytes. Unit 1 is then written
forwards starting at byte `n0`, after which the byte range `[n0, 192)` is reversed, so unit 1 reads
backwards from the last byte and the gap between the units is zero bytes. (The 0xF8 padding that
the decoder tolerates is not produced by this encoder.)

## 4. Arithmetic conventions

Exact reproduction needs the precision of each operation:

- **SSE stages** (QMF, block peaks, MDCT, bit counting, final quantisation): IEEE single precision,
  each operation rounded to nearest-even, operations in the order given.
- **x87 stages** (joint stereo, gain control, bit split, tone quantisation): the FPU works with
  53-bit mantissas, i.e. IEEE double arithmetic. Values are rounded to float32 only when stored
  into float variables; this document marks such stores with `f32(·)`. Loads of float32 values
  are exact. Conversions to integer by the C runtime truncate toward zero (`trunc`).
- **Rounding by magic number**: adding `M = 1.5·2^23 = 12582912` to a float32 value `y` with
  `|y| < 2^22` and storing the float32 sum leaves `0x4B400000 + round_half_even(y)` in its bit pattern.
  The encoder uses the low bits (masked) or the low 16 bits (sign-extended) of that pattern as the
  quantised integer. In SSE code `y` itself is a float32 product (`f32(f32(x·s) + M)`); in x87 code
  `y + M` is computed in double and rounded once to float32 (`f32(x·s + M)`).
- **Scale factor index of a group** (`sf_index4`, 0x438760): let `u` be the largest of the four
  values `(bits(x_i) << 1) mod 2^32` (the magnitude bits), `e = u >> 24`, `m = u & 0xFFFFFF`;
  `i = 3e − 364`, plus 1 if `m > 0x965FE9`, minus 1 if `m < 0x428A30`; if `i` is outside 0..63 the
  result is 0. This is the smallest index with `2^((i−15)/3)` above the peak (the thresholds are the
  mantissas of 2^(1/3) and 2^(2/3)); magnitudes of 2^16 and more wrap to index 0.
- **Quantiser scale** `qscale(sf, sel)` = MAX_QUANT[sel] / 2^((sf−15)/3) as float32, obtained as
  `QSCALE[3·sel − sf mod 3]` with `floor(sf/3)` subtracted from its exponent
  (QSCALE[3·sel − r] = f32(MAX_QUANT[sel]·2^(5 − r/3)), MAX_QUANT = 0, 1.5, 2.5, 3.5, 4.5, 7.5, 15.5,
  31.5). `floor(sf/3)` is computed as `(sf·0x2B0000) >> 23`, which is equal for sf < 128.

## 5. QMF analysis (0x435c50, inlined in [Encoder.encode_frame])

**Mathematically** (verified to float32 rounding, relative error ≈ 2·10⁻⁷): the 4-band split is the
classic ATRAC tree of two-band QMF analyses using the 48-tap prototype `W` of the decoder
(`QMF_WIN` in `atrac3_decoder.py`), scaled by 1.169770·10⁹ (the square of the ratio between Sony's
lattice coefficients and `W`) and with no extra delay. One two-band analysis of a sequence with
46 samples of history `p[0..45]` and new samples `p[46..]` is

```
for j = 0, 2, 4, ...:
    a = Σ_{i<24} W[2i]   · p[47 + j − 2i]
    b = Σ_{i<24} W[2i+1] · p[46 + j − 2i]
    low[j/2] = a + b ;  high[j/2] = a − b
```

and the tree is: (L, H) = split(x); (band0, band1) = split(L); (band3, band2) = split(H).

**Exactly** (what the encoder computes): the state is 138 floats, `s[0..91]` = the last 23 stage-1
output vectors (4 floats each), `s[92..137]` = the last 46 input samples. Let `h` be the 46 history
samples followed by the 1024 new ones (`h[k]`, k = 0..1069), and for frame step `i = 0..255` define
4-vectors `V_i(o) = (h[40+4i+o], …, h[43+4i+o])`. With the coefficient vectors `Q1[0..11]`
(0x48bfd0) and `Q2[0..23]` (0x48c090) from `lp4_tables.py`, all products/sums in float32 and
summed strictly left to right:

```
SA = Q1[11]·V(4)  + Q1[0]·V(−40) + Q1[1]·V(−36) + … + Q1[10]·V(0)
SB = Q1[11]·V(6)  + Q1[0]·V(−38) + Q1[1]·V(−34) + … + Q1[10]·V(2)
     (i.e. term k = 0..10 is Q1[10−k]·V(−4k), resp. V(2−4k); the sum starts with the Q1[11] term,
      then adds k = 10, 9, …, 0)
u  = (SA2+SA0, SB2+SB0, SA3+SA1, SB3+SB1)
g_i = (u0+u2, u1+u3, u2−u0, u3−u1)                                   # stage-1 output vector
w  = Q2[23]·g_i + Q2[0]·g_{i−23} + Q2[1]·g_{i−22} + … + Q2[22]·g_{i−1}   (left to right)
band samples of step i:  band0 = w0 + w1,  band1 = w1 − w0,  band2 = w3 − w2,  band3 = w2 + w3
```

(SAj denotes lane j of the vector SA; each lane of a vector product is an independent float32
multiply.) Stage-1 vectors with negative index come from the state.
After the frame the state holds the last 23 stage-1 vectors and the last 46 input samples.
The output is kept as an interleaved array X[t][b] (t = 0..255, b = 0..3).

## 6. Joint stereo (66 kbps only)

### 6.1 Decision [JointStereo.analyse] (0x439630)

State per band: `cur[b]`, the decision of the previous call, and `sent[b]`. For each band
b = 0..3 (in order):

```
s0 = Σ_t |X0[t][b]|,  s1 = Σ_t |X1[t][b]|          (double accumulation over t = 0..255)
e0, e1 = f32(s0), f32(s1)                           (stored; history of 3 frames is kept)
k0 = k1 = 19.697716 (f32)
if b < limit:  k0 = k1 = 2304;  if cur[b] == 0: k0 = 96;  elif cur[b] == 1: k1 = 96   # hysteresis
if   k0·s0 < s1:            d = 0       # channel 1 dominant
elif k1·s1 < s0:            d = 1       # channel 0 dominant
elif s0 + s1 < 4.791992e12: d = −3      # comparable, quiet
else:                       d = 3       # comparable
sent[b] = |cur[b]|;  cur[b] = d         # 'sent' = matrix selector transmitted in this frame
if b ≥ limit and |d| == 3:  A0 = 16·A0 + e0;  A1 = 16·A1 + e1      (doubles, A0 = A1 = f32(0.1) at the start of each call)
```

Weight index (one per frame, used by bands ≥ limit): the previous current weight becomes the
transmitted weight; then

```
T = f32(A1 + A0)
if A0 ≤ A1: (m, side) = (A0, 0) else (m, side) = (f32(A1), 8)
w = (bits(f32(m / T · 14 + M)) & 7) + side ;  if w ∈ {0, 8}: w = w + 1     (*the fix-up is unreachable*)
```

i.e. `w = side + round(14·min/sum)`, where `side = 8` means channel 1 is the weaker one.

### 6.2 Conversion of band b (0x437150)

Tables: `CHSW = [0, 2, 2, 1, 1]` [CHSW_TABLE] gives the pair `(c0, c1) = (CHSW[d], CHSW[d+1])`
for |d| = 0, 1, 3: (0, 2), (2, 2), (1, 1). `CHWT[0..23]` [CHWT_TABLE] holds the intensity divisors.
`k` runs over −8..−1 for the first 8 samples (a linear crossfade from the previous frame's
parameters to the current ones); samples 8..255 use the current parameters.

**Bands b < limit (matrixed)**, with `d = |cur[b]|` (this call's decision) and `p = sent[b]`:

```
r = 1.0
if d == 3:                                   # side gain from the energy balance of this frame
    if  e1 > 3.94·e0 and e0·256 > e1:  r = f32((e1 − e0)·1.6802721 / (e1 + e0))
    elif e0 > 3.94·e1 and e1·256 > e0: r = f32((e0 − e1)·1.6802721 / (e1 + e0))
    elif e1 == 0 or e0 == 0:           r = 1.0
    elif e1·256 ≤ e0 or e0·256 ≤ e1:   r = f32(5/3)
q = previous r of this band;  limit the change: if r − q ≥ 0.05: r = f32(q + 0.05)
                                                 if q − r ≥ 0.05: r = f32(q − 0.05)
(c0, c1) = CHSW pair of d; (p0, p1) = CHSW pair of p
c1r = c1·r;  Dc = f32((p0 − c0)/8);  Dw = (p1·q − c1r)/8
for samples t = 0..7 (k = t − 8):
    m = (x0 + x1)·0.5
    x1' = (x0 − (c0 − Dc·k)·m) / (c1r − k·Dw)
    x0' = m
for t = 8..255:  m = (x0 + x1)·0.5;  x1' = (x0 − m·c0) · (1/c1r);  x0' = m
```

(3.94 is a double constant, 256, 1.6802721 and 0.5 float constants, all arithmetic in double,
results stored as float32; `1/c1r` is a double.) Unit 0 then carries the mid signal
A = (L+R)/2 and unit 1 the difference B = (L − c0·A)/(c1·r). With d = 3: B = (L−R)/(2r).

**Bands b ≥ limit (intensity)**: only channel 0 is changed. With `d = |cur[b]|`, the current and the
previous weight indices `wc`, `wp` are mapped to table indices: for d = 0 → `(w XOR 8) + 8`,
d = 1 → `w + 8`, d = 3 → `w AND 7` (d = 2 does not occur):

```
tc, tp = CHWT[wc'], CHWT[wp'];  D = (tp − tc)/8
t = 0..7:   x0' = (x1 + x0) / (tc − k·D)
t = 8..255: x0' = (x1 + x0) · (1/tc)
```

CHWT[8..15] = 2l/7 and CHWT[16..23] = 2√(2 − (l/7)²) (l = 0..7) are the decoder's weights; d = 0
divides by the weight of channel 1, d = 1 by the weight of channel 0. CHWT[0..7] (√2 … 2) is used for
d = ±3.

The transmitted parameters of a frame are `sent[0..3]` (matrix selectors) and the transmitted weight
index; both describe the samples converted in the previous call, which are the samples completed
by this frame's MDCT (section 1).

## 7. Block peaks and gain control

### 7.1 Block peaks (0x435ab0)

`P[b][g] = max_{t=8g..8g+7} |X[t][b]|`, g = 0..31, computed on the band samples after the stereo
conversion.

### 7.2 Per band state

Per unit and band: the 32 peaks of the previous frame (`env`), `peak` (see below), `last` (the peak of
the last 32 samples of the previous analysis), the gain points `loc[k]`, `lev[k]` and their count.

### 7.3 Detection (0x439bb0, called by 0x439820 for b = 3, 2, 1, 0)

Inputs: `A[0..63]` = previous peaks (32) followed by the new peaks (32); `mode`:
for an ordinary unit `mode = b`; for the joint-stereo unit `mode = 5` if the band's decision changed
(`sent[b] ≠ cur[b]`, section 6.1), else `0` if `sent[b] ≠ 3`, else `−1`. Also `peak_old` = the band's `peak` from
the previous call. Constants: FLOOR = 1.8718718e10, MIN = 4.6796792e10, R = 1.85 (all float32).

```
blk[j] = max(A[4j..4j+3]) for j = 0..7;   blk[−1] = last;   last = blk[7]
end = 64 if mode ≤ 0 else 8·(8 − mode)
ref = max(FLOOR, blk[7], A[32..end−1]);  thr = f32(R·ref)
budget = 4; slots = []                                     # decays, filled from the end
for j = 7 down to −1:                                      # backwards over 32-sample blocks
    v = blk[j]
    if v < ref: continue
    if v > MIN and v > thr:                                # level drops after block j
        e = 4(j+1)
        loc = e if (j == −1 or A[e] ≥ thr or A[e−1] ≥ thr) else (e−1 if A[e−2] ≥ thr else e−2)
        x = min(L(v/ref), budget);  budget −= x;  prepend (loc, −x) to slots
        if budget ≤ 0 or len(slots) == 2: stop
    ref = v;  thr = f32(v·R)
```

where `L(z) = exponent of f32(z·√2) = floor(log2(√2·z))` (√2 as float32 1.4142135; `z` a double).
(The block j = 7 can never trigger because `ref ≥ blk[7]`.)

Then the attack search over the part of the previous frame before the first decay point
(`stop = loc` of the first slot, or 32):

```
room = min(15, 131 − biasedexp(f32(peak_old · 9.2225066e−15))) − budget
points = []
if room > 0:
    f = 2.0 if mode == −1 else f32(1.6)
    ref = max(peak_old, A[0]);  thr = f·ref                     (thr is a double)
    for i = 0 .. stop−1:
        a = A[i+1]
        if a < ref: continue
        if a > MIN and a > thr:
            x = L(a/ref)
            if points and last point is at i−1 and x ≥ its level:    # merge with it
                remove it, room += its level, x += its level
            x = min(x, room);  room −= x;  append (i, x)
            if len(points) == len(free slots before the decays) or room ≤ 0: stop
        ref = a;  thr = a·f
```

("free slots" = 7 − number of decays; the loop also ends when the attack list reaches that size.)
The final list is the attacks followed by the decays; levels are made absolute by summing from the
end: `lev[last] = 4 + x_last`, `lev[k] = lev[k+1] + x_k` (attack steps are positive, decay steps
negative), and `lev[count] = 4` is appended. Locations are in units of 8 samples (0..31).

After detection `peak = max(FLOOR, previous env)` is stored and `env` = the new peaks.

**Band-0 copy rule.** If band 0 (ordinary unit) found no point, had no points in the previous frame,
band 1 has points, and `peak ≤ 7.6671869e13`: if the levels of band 1 (together with 4) span at least
2, and none of `A[0..loc1[0]]` exceeds 7.667e13 (*always true here*), band 0 gets one point
`(loc1[0], 5)`.

### 7.4 Application

The previous frame's 256 band samples are multiplied in place (float32 products) by the gain curve:
for each point k (with `lev[count] = 4` as terminator), samples from the end of the previous ramp
up to and including sample `8·loc[k]` are multiplied by `2^(lev[k]−4)` (skipped when the level is 4),
then the 7 samples `8·loc[k]+j`, j = 1..7, by `2^((8·lev[k] + j·(lev[k+1]−lev[k]))/8 − 4)` (exactly the
float32 of that power; Sony builds it from a 3-bit fraction table). The scaled samples are used only
for the transform; the stored copy for the next frame is the unscaled band signal.

The count, locations and levels are transmitted in this frame (4-bit level, 5-bit location).

## 8. MDCT (0x435380, inlined in [Encoder.encode_frame])

**Mathematically** (verified to float precision): with `v[0..255]` = the gain-scaled previous-frame
samples of the call before (first half of the window) multiplied by `2^(lev0−4)` where `lev0` is the
level of this frame's first gain point of the band (4 if none), and `v[256..511]` = this call's
gain-scaled previous-frame samples,

```
X[k] = Σ_{n<512} w[n]·v[n]·cos(π/256·(n + 0.5 + 128)·(k + 0.5)),   k = 0..255
w[n] = w[511−n] = K·(1 + sin((n + 0.5)·π/256 − π/2)) for n < 256,   K = 8.34727·10⁻¹³
```

and the stored spectrum of bands 1 and 3 is reversed (`X[255−k]` at position k). Together with the
QMF scale this puts the coefficients on the decoder's scale (the decoder multiplies quantised values
by 2^((sf−15)/3)/MAX_QUANT and applies an unnormalised IMDCT).

**Exactly** (all float32, per band lane; `cur[0..255]` = this call's previous-frame samples,
`state[0..127]` from the previous call, `g = 2^(lev0−4)`):

1. Pre-rotation, for j = 0..127 with `p = PERM[j]`:
   `a = cur[128+p]`, `b = cur[127−p]`, `s = state[j]·g`,
   `state'[j] = R1[p]·(a·R0[p] − b)`, `t = R2[p]·(b·R0[p] + a)`, `d = R3[p]·(t − s)`,
   `Y[2j] = (t + s) + d`, `Y[2j+1] = d`.
   PERM is the order f(128) with f(1) = (0), f(2N) = (x, 2N−1−x for each x of f(N)), e.g. 0, 127, 63,
   64, 31, 96, …. R0[p] = f32(W(128+p)/W(127−p)) with W(n) = 1 + sin(((n+0.5)/256 − 0.5)π) and
   R3[p] = f32(1/(2cos(π(2p+1)/512))) (both verified exact); R1, R2 contain the window, the overall
   scale and further rotation factors and are taken from the table (0x48c280, `MDCT_TW`).
2. First butterfly stage over groups of four `(y0, y1, y2, y3)` = Y[4j..4j+3], j = 0..63, with
   `c = C(2, j, 64)`: `e3 = c·(y1 − y3)`, `e0 = c·(y0 − y2)`, result
   `((y0 + y2) + e3, (y1 + y3) + e0, e0, e3)`.
3. Stages n = 4, 8, 16, 32, 64: Y is divided into blocks of 2n; for block k (k = 0 .. 256/(2n) − 1) with
   `c = C(n, k, 256/(2n))`, `lo = Y[base .. base+n)`, `hi = Y[base+n .. base+2n)`:
   `hi' = c·(lo − hi)`, `lo' = (lo + hi) + reverse(hi')` (element m of `lo'` adds element n−1−m of `hi'`).
4. Final stage: `hi' = 0.70710677·(Y[0..127] − Y[128..255])`, `lo' = Y[0..127] + Y[128..255]`.
5. Output: for i = 0..127, `Z = Y[i] + Y[255−i]`; even bands: `spec[i] = Z`, `spec[255−i] = Y[255−i]`;
   odd bands: `spec[i] = Y[255−i]`, `spec[255−i] = Z`.

Here `C(n, k, B) = f32(1/(2cos(π·n·(2·f(B)[k] + 1)/512)))` (verified exact; the table at 0x48ca80
holds these 127 values). The new state replaces the old one.

The spectrum of a unit is stored band-major: coefficient `256b + k`.

## 9. Bit split between the units [split_bits] (0x4378c0, 66 kbps only)

With A = spectrum of unit 0, B = spectrum of unit 1, `nsub0` = 26, `nsub1` = 13, `bpc` = 96,
`base0` = 1133, `d0` = sent[0] (the band-0 selector transmitted this frame) and `g1` = number of
gain points of unit 1 in band 0:

```
S0 = f32( Σ_{i = min(nsub0−1, 18) down to 0} SB_WIDTH[i] · max_{k in subband i} |A[k]| )   (double sum)
S1 = Σ_{i=0}^{nsub1−1} SB_WIDTH[i] · max |B[k]|;   T = f32(2·S1);  if d0 ≠ 3: T = f32(2·T)
total = 16·bpc − 59 = 1477
if T < 1 or (d0 ≠ 3 and T < 4):            budget0 = total
elif g1 == 0 and ((3T < S0 and d0 == 3) or 6T < S0):
    spare = total − base0
    if not (12T > S0): x = spare
    else:
        if not (10T ≥ S0): spare = trunc(1.5·spare)
        elif not (5T ≥ S0): spare = trunc(1.25·spare)
        x = trunc((S0 − T)/(2·S0) · spare)
    budget0 = base0 + (x AND ~7)
elif T ≥ S0 or d0 ≠ 3:                       budget0 = base0
else:                                        budget0 = base0 + (trunc((S0 − T)/(10·S0) · (total − base0)) AND ~7)
```

After unit 0 has been packed into `n0` bytes, unit 1 gets `budget1 = (2·bpc − n0)·8 − 27`.
A unit whose budget is below 40 bits is written as an empty unit: 1 QMF band, its band-0 gain count is
forced to 0 (this also affects the next frame's band-0 state), no tones, 1 subband with selector 0.

## 10. Bit allocation [allocate] (0x437b40)

Inputs per unit: the spectrum `x[0..1023]` (the allocator works on a copy that tone extraction
modifies), `budget`, `nsub_max`, whether it is the joint-stereo unit (`js`), this frame's gain point
counts per band, the previous frame's band-0 count, and (js unit) the transmitted weight.
Outputs: `nsub`, selectors `sel[i]`, scale factors `sf[i]`, `nbands`, tonal groups, residual spectrum.
SB_START/SB_WIDTH are the 32 standard subbands; "group" g = coefficients 4g..4g+3.

### 10.1 Scale factor survey

`gsf[g] = sf_index4(x[4g..4g+3])` for g = 0..255. Then for i = 0..31 (running variables `cnt7` =
number of groups with gsf > 7 so far, `cum` = running total, `top` = largest subband maximum so far):

```
sum_i = Σ gsf over the subband;  mx_i = max gsf;  cnt7 += #(gsf > 7)
if top < mx_i:
    top = mx_i
    if not js:                              # a new loudest subband may extend the bandwidth
        ok = cnt7·16 < budget if (SB_START[i] AND 0x80) and i < 22 else cnt7·8 < budget
        if ok: nsub = max(nsub, 28 if i ≤ 26 else min(i+2, 32))
elif sum_i < MIN_SFSUM[i] or mx_i < 3:  sum_i = mx_i = 0
sf[i] = mx_i;  C[i] = cnt7;  cum += sum_i
```

(nsub starts at nsub_max.) `shift = 11 if 2·budget < 12·C[nsub−1] else 10`.

`nbands = floor((SB_START[nsub] + 255)/256)`; gain bits `G = 3·(3·Σ_{b<nbands} count_b + nbands)`;
`left = budget − 3·nsub − G`.

### 10.2 Path A: quiet frames (all-tonal, 0x4387c0)

Taken when `left > 32·C[nsub−1]` (few loud groups). Two tone groups are set up: group 0 with step 5,
group 1 with step 7, both 4 values per tone, tonal coding mode 1 (CLC). `cost = 22` (plus 3 per band
if `nbands` had to grow to `floor(SB_START[nsub]/256)`, *unreachable*). Repeated passes over
k = 0 .. SB_START[nsub]−1 until a pass adds no tone:

```
if |x[k]| ≥ 2.0:
    pos = min(k, 1020);  s = sf_index4(x[pos..pos+3]);  grp = 1 if s > 30 else 0;  blk = k >> 6
    if group grp's block blk already has 7 tones:
        if the other group's block blk also has 7: replace: among the 7 tones of grp in blk take the
            one with the smallest sf below s (first in list order); if none: skip k;
            otherwise add it back into x (10.6), remove it from the list and re-use its record
            (with its old position) instead of a new tone
        else grp = other group
    first tone of grp in the 256-coefficient band blk>>2: cost += 12 (band flag set)
    append the tone to grp's block list, quantise it at its position (10.5) with grp's step
    cost += 28 + 8·grp;  if cost > left − 200 or 64 tones: stop everything
    k += 3
```

Afterwards for every subband i < nsub: recompute gsf of its groups on the residual, `sf[i] = max`,
`sel[i] = 1 if sf[i] ≠ 0 else 0`, `P[i] = 341·sf[i]`. `left −= cost`. Continue with 10.7.

### 10.3 Path B: preliminary selectors

```
for i < nsub:  m = sf[i];  sel[i] = 0
    if m ≥ 3: m = 256·m − cum;  sel[i] = clamp(m >> shift, 1, 7)
    P[i] = m                      (a priority value, arithmetic shifts throughout)
```

If not js: tone extraction (10.4), `left −= its cost`. Then for each i < nsub:

- if `sf[i] < 0` (the subband lost tones): recompute gsf of its groups and `mx`, `sum` on the
  residual; `sf[i] = mx`; if `sel[i]` (here: number of coefficients taken by tones) is 0,
  `P[i] += 256`, else `k = min((sel[i]−1) >> 1, 4)`, `D = TONE_PENALTY[k]·(mx + old_sf)` (old_sf is
  negative), doubled when `budget < 1100` (*unreachable at 66/105 kbps*); otherwise, if
  `SB_WIDTH[i] ≤ 16 and P[i] > 0x1800`, `D = min(2D, P[i] − 0x400)`; then `P[i−1] −= D>>3` (i > 1),
  `P[i+1] −= D>>2` (i < 31), `P[i] −= D`. Finally `sel[i] = clamp(P[i] >> shift, 1, 7)` if
  `mx > 7 and sum > MIN_SFSUM[i] and (P[i] ≥ 0x400 or budget > 0x708)`, else 0.
- elif shift = 11: if `SB_WIDTH[i]·m < 5·sum_i` then (if m ≥ 62 skip this subband) `sf[i] = m+1`;
  then if `sf[i] < 62`: `sf[i] += 1` (coarser quantisation of dense subbands).
- else if `SB_WIDTH[i]·m > 8·sum_i`: `P[i] += 256` (peaky subbands get priority).

### 10.4 Tone extraction (0x438e60)

One tone group: step 3, 4 values per tone. `hdr = 8`, counters `n` (tones), `bits`, `nzero`
(tones whose 4th value is 0). Subbands i from 0 (from 8 if `left < 1100`) while i < nsub, stopping
when `n > 64 − SB_WIDTH[i]/16` or `24n + hdr > left − 600`; subbands with `sel[i] = 0` are skipped:

```
cls = 0 if P[i] < 0x500 else sel[i];  band = 1 if SB_WIDTH[i] == 32 else SB_START[i] >> 8
t = max(0, sf[i] − TONE_SF_OFFSET[cls + 8·band]);  start = SB_START[i]
if n > 0 and pos(last tone) + 4 > start:                 # previous tone overlaps this subband
    sf[i] = −sf[i];  o = pos(last) − start;  sel[i] = o + 4;  sel[i−1] += o;  start = pos(last) + 4
peaks = positions k from start (at least once, then while k < SB_START[i+1]) with
        |x[k]| ≥ SFVAL[t]; after a hit k skips 3 extra; more than (SB_WIDTH[i]+8)/16 peaks → none
if peaks: if sf[i] ≥ 0: sf[i] = −sf[i]; sel[i] = 0
for each peak p:
    c = min(p, 1020);  nv = 4
    if n > 0 and pos(last) + 4 == c:                       # adjacent: slide back while useful
        repeat: if |x[c+3]| > |x[c−1]|: break;  c −= 1;  nv −= 1  until nv ≤ 1
    if the 64-block c>>6 already holds 7 tones: skip peak
    new tone at c, quantise (10.5) → cost b
    if its first value is 0 and c ≤ 1020: shift it right past its leading zeros (e = index of
        the first nonzero of values 1..3, 4 if none); if c+e ≥ SB_START[i+1]: undo the tone
        (add back, 10.6), drop it and leave this subband; if c+e > 1020: keep as is;
        else move: values shift left by e, zeros appended, pos += e
    bits += b;  sel[i] += nv;  new 256-band → hdr += 12;  append to block list;  if 4th value 0: nzero += 1
```

After all subbands: no tone → cost 0. One tone at position < 128, or one tone with `left < 1100`:
add it back, cost 0 (the negative sf and changed sel stay). Otherwise `hdr += nbands`, one group is
coded and the representation is optimised:

```
V = bits (VLC), Cc = 24n (CLC); z0 = code length of value 0: VLC 1, CLC 3
if nzero ≠ n:
    nz = 2n − 3·nzero;  x12 = 12(n − nzero)
    Ve = z0v·nz + x12 + bits;  Ce = z0c·nz + x12 + 24n
    if min(bits, 24n) > min(Ve, Ce) and split succeeds (n' = new count):
        V = z0v·n' + Ve;  C = z0c·n' + Ce;  n = n';  shrink = yes
else shrink = yes
if shrink: k = (values removed by shrink)·n;  V −= z0v·k;  C −= z0c·k
mode = CLC if V > C else VLC;   cost = hdr + min(V, C)
```

*Split* (0x439400): every tone whose 4th value is nonzero is split into itself with values (q0, q1, 0,
0) and a new tone at pos+2 with (q2, q3, 0, 0) and the same sf, appended after the others. It fails
(nothing changed) if for any such tone the band of pos+2 has no tones yet, its 64-block would exceed
7 tones, the total would reach 64, or pos+2 > 1021.

*Shrink* (0x439360): with `v` = values per tone − 1 (initially 3), repeat: for every tone (from last
to first) whose value v is nonzero, move it one position right (values shift left, value v stays);
mark "stop" if value v−1 is nonzero and (value 0 is nonzero or pos mod 64 == 63); then `v −= 1`;
until "stop". The number of removed values is 3 − v.

### 10.5 Tone quantisation (0x438af0, [best_sf_tone] 0x438bb0)

For a tone at `pos` with step `s` (x87):

```
best = 1e10;  lo = sf_index4(x[pos..pos+3]);  hi = min(lo + 3, 64) − 1
for f = hi down to lo:
    sc = qscale(f, s);  err = 0 (float32);  acc = 0
    for the 4 values: y = sc·x;  q = low16(f32(y + M));  acc |= q;  if |y − q| > err: err = f32(|y − q|)
    if sc·best > err:
        best = err/sc;  sf = f
        while acc is even and sf < 61: acc >>= 1 (arithmetic); sf += 3      # drop a factor 2
sc = qscale(sf, s);  inv = f32(1/sc)
for the 4 values: q = low16(f32(sc·x + M));  x = f32(x − q·inv);  value = q AND QMASK[s]
cost = 12 + Σ VLC length of value (mode 0 table of step s)
```

### 10.6 Adding a tone back (0x438dc0)

`inv = 1/qscale(sf, s)` (double, not rounded); each stored value is sign-extended from the CLC width
of step s (3, 4 or 6 bits for steps 3, 5, 7) and `x = f32(q·inv + x)`.

### 10.7 Pruning and priority offset (both paths)

```
for i < nsub with sel[i] ≠ 0:
    if (i > 1 and P[i] < P[i−1] − 0xF00) or (i < nsub−1 and P[i] < P[i+1] − 0x1400): sel[i] = 0
    else: W += SB_WIDTH[i];  if sel[i] == 7: W7 += SB_WIDTH[i];  E += est(i)
drop trailing subbands with sel 0 (left += 3 each, nsub shrinks)
left += drop_bands()
L10 = 10·left
if L10 ≤ E: width = W  elif W − W7 ≠ 0: width = W − W7  else no offset (off = 0)
if width defined:
    t = L10 if left < 600 else 9·left;  diff = t − E;  off = trunc(1024·diff / (10·width))
    if diff < 0 and shift == 10:
        off += off >> 3;  P[0] −= off;  P[1] −= (188·off) >> 8
        P[2..7] −= (97·off) >> 8;  P[8..17] −= (74·off) >> 8
        for k < 18: if P[k] + off < 0 and sf[k] < 63: sf[k] += 1
        for 18 ≤ k < nsub: if P[k] + off < 1024: sf[k] = min(63, sf[k] − ((P[k] + off − 1024) >> 10))
for i = nsub−1 down to 0 with sel[i] ≠ 0:
    sel[i] = clamp((P[i] + off) >> shift, 1, 7);  L10 −= est(i);  P[i] = L10
```

`est(i)` (0x4385c0, tenths of a bit) = `EST_BITS[sel]·SB_WIDTH[i] + 60 − EST_SAVE[sel]·#{groups g of i
with gsf[g] < sf[i] − EST_THR[sel]}`.

`drop_bands()` (0x438690): for the joint-stereo unit while the transmitted weight index AND 7 ≠ 7,
nothing. Otherwise starting at the top band, while band k > 0 has no gain points and no tone group
uses it, `k −= 1`; `nbands = k + 1`; returns `(groups + 3)·(bands removed)`.

### 10.8 Exact counting

```
used = 0
for i < nsub with sel[i] ≠ 0:
    if SB_WIDTH[i] > 2·(left − used) − 16: sel[i] = 0; continue
    d = P[i] − 10·used
    if sel[i] < 7 and d > 4000: sel[i] += 1
    elif i ≠ 0 and d < 300: if sel[i] > 1: sel[i] −= 1  elif sf[i] ≠ 63: sf[i] += 1
    b = bits(i);  P[i] = b
    if b == (code length of 0)·SB_WIDTH[i] + 6 (all values quantise to 0): sel[i] = 0
    else: used += b;  order[i] = sf[i]
drop trailing zero-selector subbands (left += 3 each);  slack = left − used
if not js and (previous band-0 gain count + this one) ≠ 0:
    for k < 5: if 0 ≤ order[k] ≤ 64: order[k] = 64      # transient frames: low subbands first
```

`bits(i)` (0x436740 [bitcount]) = 6 + the VLC lengths of the subband's values
quantised with (sel, sf) (section 11).

### 10.9 Final fit (in order of decreasing `order`, ties → lowest index)

```
spent = 0
repeat: pick i with the largest order ≥ 0 (none → done);  order[i] = −1;  b = P[i]
    if spent + b > left:
        clear unprocessed subbands above i from the top down (*unreachable*, see 14), nsub follows
        while spent + b > left:
            sf[i] += 1;  if sf[i] ≥ 64 or b ≤ zl·(SB_WIDTH[i]>>1) + 6: break
            b = bits(i)                                  # zl = 2 × length of the 0 code
    elif sel[i] < 7 and STEP_COST[sel[i]]·(SB_WIDTH[i]>>1) < left − spent − 7:
        ns = sel[i] + 1;  if slack > 4·SB_WIDTH[i] and ns < 7: ns += 1
        b2 = bits(i) with ns;  if b2 − b ≤ max(slack, 0): sel[i] = ns; slack += b − b2; b = b2
    spent += b;  if spent > left: spent −= b; slack += b; sel[i] = 0
```

Subbands dropped here stay inside `nsub` with selector 0.

## 11. Quantisation and code lengths (SSE)

For subband i with (sel, sf): `y = f32(x·qscale(sf, sel))`, `v = f32(y + M)`, `q = bits(v) − 0x4B400000`.
For sel ≥ 2 the code is VLC table `(0, sel)` at index `q AND QMASK[sel]` (QMASK = 7, 7, 15, 15, 31,
63 for sel 2..7). For sel 1 pairs (q0, q1) use index `(q0 AND 3)·4 + (q1 AND 3)` of table (0, 1);
for counting, each 4 values cost `SEL1_LEN[pattern]` (pattern bit j = value j nonzero), which equals
the sum of the two pair codes. Values outside the selector's range alias to other table entries
(usually the zero code); the allocator's scale factors normally prevent that.

## 12. Packing (0x437490)

```
[js unit]  byte: weight(4) · sel0(2) · sel1(2)
           byte: sel2(2) · sel3(2) · 11 (id) · nbands−1(2)
[unit 0]   byte: 101000 (id 0x28) · nbands−1(2)
for b < nbands: count(3), then per point level(4) location(5)
tone groups(5); if > 0: mode(2) and per group: band flags (nbands bits, band 0 first),
    values−1 (3), step (3), then for each 64-block of the flagged bands: count(3) and per tone
    sf(6), pos mod 64 (6), the values with table (mode, step)
nsub−1 (5), coding mode 0 (1), sel (3) × nsub, sf (6) for each sel ≠ 0, then the spectral codes of
each coded subband in order
```

The weight and selectors are the transmitted values (section 6). Section 3.3 gives the frame layout.

## 13. Tables

Small tables (all in `lp4_tables.py` with addresses):

```
SB_START  0 8 16 24 32 40 48 56 64 80 96 112 128 144 160 176 192 224 256 288 320 352 384 416 448
          480 512 576 640 704 768 896 1024
SB_WIDTH  8×8, 8×16, 10×32, 4×64, 2×128 (in subband order)
MIN_SFSUM 7 5 5 4 4 4 4 3 3 3 3 4 5 5 5 6 6 7 7 8 10 13 17 22 28 35 49 74 109 155 250 441
by selector 0..7:
  ZERO_LEN2   441* 1 2 2 2 4 6 6        STEP_COST 1 2 2 2 4 6 6 –
  EST_BITS    –  15 20 25 29 35 45 55   EST_SAVE  – 40 40 60 76 60 60 100
  EST_THR     –   3  5  7  9 12 15 18   (*unused entry)
TONE_PENALTY  −225 −266 −307 −317 −1024
TONE_SF_OFFSET[cls + 8·band]: band0 11 9 9 9 9 9 10 13 | band1 9 8 8 8 8 9 10 13 |
                              band2 8 7 7 7 8 9 10 13 | band3 7 7 7 7 8 9 9 9
SEL1_LEN  2 5 4 6 5 8 7 9 4 7 6 8 6 9 8 10
CHSW 0 2 2 1 1;  CHWT[0..7] 1.4142135 1.5498369 1.6707656 1.7762830 1.8650550 1.9348577 1.9820011 2.0
SFVAL[i] = 2^((i−15)/3) (float32, identical to the decoder's SF table)
```

Large tables reproduced in `lp4_tables.py`: QMF coefficient vectors Q1 (12×4) and Q2 (24×4) – all
entries are taken from the 24 lattice coefficients of Sony's QMF (1.0, 6.296898 and the 22 values of
`QMF_HALF_SONY` in `at3_tables.py`) –, MDCT rotation constants (128×4, of which columns 0 and 3
have the closed forms of section 8), the 127 DCT factors (closed form in section 8), and the VLC/CLC
code tables per (mode, selector) as (code, length) by masked index. The code tables agree with the
decoder's tables in `at3_tables.py` for every in-range value; masked indices that correspond to no
legal value hold the zero code.

Scalar constants: stereo 19.697716, 2304, 96, 4.791992e12, 1.6802721, 3.94 (double), 0.05 (double),
5/3; gain 1.8718718e10, 4.6796792e10, 1.85, 1.6, 2.0, √2, 9.2225066e−15, 7.6671869e13; tones 1e10, 2.0.

## 14. Edge cases and unreachable code

- Peak magnitudes ≥ 2^16 give scale factor index 0 (wrap) in `sf_index4`; with 16-bit input the
  largest index observed is 60, so this needs crafted spectra (tested by injection).
- Tone positions are clamped to 1020 so that 4 values fit; tones near the top can therefore keep a
  leading zero value (section 10.4).
- *Unreachable* in these two modes (kept literally in the model): the weight fix-up for 0/8 (the
  decision thresholds bound the energy ratio, so min/sum ≥ 1/20.7); budget < 1100 for a unit that
  extracts tones; nbands growth in path A; the band-0 copy rejection test; and the clearing loop of
  10.9 – exact counting stops adding subbands once `used` exceeds `left`, so the overshoot is at most
  the bits of the last coded subband, and the final fit can only cross the limit at that subband,
  above which nothing remains unprocessed.
- Not modelled: the tool's loop mode for WAVs with `smpl` chunks; other rows of the mode table
  (94/132/146/176 kbps) which the tool assigns to other encoders or does not offer.

## 15. Verification

- `tools/verify_lp4.py`: library-level frames (captured at the call return sites) and complete .at3
  files compared byte for byte with the emulated `psp_at3tool.exe`. Results at 66 and 105 kbps:
  100% identical frames on the six `gen_signals.py` signals, eleven `gen_lp4_signals.py` signals
  (one-sided stereo, harmonic combs, quiet tones, clipped squares, clicks, panned noise, DC, odd and
  exact lengths, tones near 16.5 kHz, broadband pressure), and the PSP SDK music `bgm01_1`, `bgm01_2`
  (592 library frames each), `overload_short` (2371) and `overload` (3:14; the complete 8380-frame
  .at3 file is identical at both rates).
- Stage-by-stage bit-exact comparisons (QMF, stereo conversion, peaks, gain points and modified
  samples, MDCT, budgets, allocation state) during development.
- Crafted-spectrum injection directly before the allocator (32 452 sound units with extreme and
  sparse spectra, tones at coefficients 1012–1023, scale factors 62/63): identical packed units.
- Line coverage of the model: every line is executed by the tests except the unreachable ones listed
  in section 14.
