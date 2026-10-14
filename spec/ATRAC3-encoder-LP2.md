# ATRAC3 encoder algorithm 1 (MiniDisc LP2, 132 kbps; 66 kbps mono)

**Chapter 7 of the [ATRAC3 specification](ATRAC3.md).** Derived from `psp_at3tool.exe` 3.0.0.0 (functions 0x43C520–0x441640).
Names in `monospace` are the corresponding functions of the reference model.

**Status.** Every stage below has been reimplemented in [`reference/atrac3_encoder_lp2.py`](reference/atrac3_encoder_lp2.py)
(+ `lp2_tables.py`, `lp2_analysis.py`, `lp2_alloc.py`, `lp2_bitstream.py`) and checked against
the emulated Sony encoder. The model's output files are **byte-identical** to at3tool's for every
test input (§16). Sections are therefore **[verified]** unless a statement is marked **[inferred]**,
which means the case never occurred in testing and is described from code reading.

Notation: `f32(x)` rounds to IEEE single precision, `f64` is double precision. `trunc` rounds
towards zero. Arrays are 0-based. "Line" means an MDCT coefficient index 0..1023 in the
bitstream's order (QMF band *b* covers lines 256·b … 256·b+255).

---

## 1. Overview

Algorithm 1 codes each channel independently into a 192-byte sound unit (1536 bits) per 1024-sample
frame. The LP2 frame is the two units side by side (384 bytes). The 66 kbps mono mode is the same
encoder fed with the mono signal in both channels, keeping only unit 0 (§2.4).

```
1024 PCM samples (one channel)
   │  int16 → float32
   ▼
QMF tree (3 × qmf_split) ── 4 bands × 256 samples (new block)
   │  band history: [old | current | next] (3 × 256)
   ├─► gain control analysis on [current | next] (band_gain_points, adjust_band0_gain) → gain points
   ▼
MDCT of [old | current] divided by the gain window
   (gain_window, mdct256)                      → spectrum A (1024 lines, gets quantised)
MDCT of [old | current] without gain window    → spectrum B (decisions only; = A if no gain points)
   ▼
encode_unit (one sound unit, 1536 bits):
   1. scale-factor index of every 4-line group of A and of B
   2. tonal components: either around one dominant isolated peak (single_tone_path)
      or wherever B is far above its mean (multitone_path); subtracted from A
   3. per-subband scale factors, initial word lengths from the scale-factor profile
      and a simple masking rule (word_lengths)
   4. A normalised per subband; global offset search to fit the bit budget
   5. dead-zone widening and word-length reduction if still too large
   6. refinement passes that spend the remaining bits
   7. scale-factor correction
   ▼
pack_unit → 192 bytes
```

The unit output by call *k* describes the input of call *k*−1 (the MDCT is one 256-sample block
behind the QMF output because gain control needs one block of look-ahead).

## 2. Configuration and framing

### 2.1 Fixed parameters

| parameter | value | origin |
|---|---|---|
| channels | 2, coded independently | encoder set-up: only accepts the codec-table code 0x84 with channel mode 1 |
| unit size | 192 bytes = 1536 bits per channel | |
| bandwidth | 15100 Hz | constant 0x3AFC |
| line limit | ⌊15100 · 2048 / 44100⌋ = 701 | `coded_layout` |
| default subbands coded | 29 (smallest *k* with start[*k*] ≥ 701 is 29; lines 0–703) | |
| default QMF bands coded | ⌊(701 + 255) / 256⌋ = 3 | |

The subband count and band count are reset to 29 and 3 at the start of every frame; the single-tone
path (§10.3) may raise them to 30/3 or 32/4 for one frame, and the final subband count is the last
subband with a non-zero word length (§10.10).

### 2.2 Library call

Each call receives up to 1024 samples per channel (`LP2Encoder.encode_call`). A short block is
zero-padded to 1024. The library keeps the per-channel sample counts of the last three calls
(initialised to 1024) and returns an **end flag** = (count of the call two calls earlier < 886).
Samples are converted exactly: `x = f32(int16)`.

The original library reports zero output bytes for call 0, although it advances all analysis
state. The executable model's `encode_call` deliberately returns the computed bytes even on
call 0; its file wrapper discards them. A library-compatible caller must suppress these bytes.
The end-flag call still produces a frame; stop after writing it. Do not call again after completion.

### 2.3 Tool framing (at3tool `-e`, no loop points)

For a WAV with N samples per channel (the tool requires N ≥ 6144):

1. The encoder input stream is 955 zero samples followed by the N source samples, cut into
   1024-sample calls; the last call carries the remainder (1..1024 samples).
2. After the data, calls with 0 samples are made until a call returns the end flag. That call's
   output is still written. (If the last data call had fewer than 886 samples, two flush calls are
   made, otherwise three.)
3. The output of call 0 is discarded by the library (it reports 0 bytes) and the output of call 1
   is discarded by the tool. Every later call writes one frame.

Hence the file holds (number of calls − 2) frames. The container (`fmt `, 12-byte `fact`
{N, 1024, 1024}, `data`) is described in ATRAC3.md §3. Input WAVs with a `smpl` chunk switch the
tool to its loop framing ([ATRAC3.md](ATRAC3.md) §9.2); that path is not part of this description.

### 2.4 Mono (66 kbps)

The tool duplicates the mono input into both channels, runs this encoder and writes only the first
192-byte unit of each frame (192-byte frames). Because the channels share nothing, this equals
encoding the mono signal as channel 0. The `fmt ` chunk has 1 channel, block align 192.

## 3. Arithmetic conventions

The reference runs on the x87 FPU with precision control set to 53 bits. Consequently:

* A value kept in an FPU register behaves as an IEEE **double**; a value stored to memory as
  `float` is rounded to **float32**. The description says which. A product of two float32 values
  is exact in double precision.
* "Truncate" means `ftol`: round towards zero.
* The quantisers compute `trunc(v + 31.5) − 31`, which equals round-half-up for v > −31.5.
* Comparisons between a float32 variable and a decimal constant (1.15, 0.9, …) are made in double
  precision with the constant as a double.

An implementation that wants Sony's exact output must follow these rules and the stated operation
orders; an implementation that only wants equivalent quality can use any precision.

## 4. QMF analysis (`qmf_split`)

The same 48-tap prototype *h*[0..47] as the decoder ([Appendix A](ATRAC3-tables.md); table 0x490FC0) splits a signal
x (even length *n*) into two half-rate signals. Each split keeps the last 46 input samples as state
(initially zero). With s = (46 state samples) ++ x:

```
for j in 0 .. n/2−1:
    P = 0; for k = 47, 45, …, 1:  P = P + h[k]·s[2j+k]          (double accumulation)
    Q = 0; for k = 46, 44, …, 0:  Q = h[k]·s[2j+k] + Q
    low[j]  = f32(Q + P)
    high[j] = f32(P − Q)
```

The band tree: split the 1024 input samples into L (512) and H (512); split L into band 0 (low)
and band 1 (high); split H into band 3 (low output) and band 2 (high output). Each split has its
own state (stage 1, stage 2 low, stage 2 high).

## 5. Band history and timing

Each band keeps the last three 256-sample blocks: **old** (from two calls ago), **current**
(previous call) and **next** (this call). After the QMF the history shifts by one block.

* The MDCT of this call transforms [old | current].
* Gain control describes transitions inside **current**, looking at [current | next] plus the
  last 32 samples of old.

## 6. Gain control (`band_gain_points`, `adjust_band0_gain`)

Per band and frame the result is a list of up to 7 gain points (level code 0..15, location 0..31
in 8-sample units), exactly as transmitted ([ATRAC3.md](ATRAC3.md) §4.3), plus a float32 **peak** used by the
next frame. Level code *L* means a log2 gain of *L* − 4 relative to the end of the block.

### 6.1 Measurements

* `blk[i]`, i = 0..63: maximum |x| over the 8-sample block i of [current | next] (float32).
* `grp[g]`, g = 0..15: max of blk[4g .. 4g+3].
* `pre`: max |x| over the last 32 samples of old.
* `prev.peak`: the previous frame's peak (= max of its blk[0..31], i.e. the peak of today's old block).
* `prev_lev`: max(0, largest *L* − 4 among the previous frame's points).
* This frame's **peak** = max(0, blk[0..31]) (stored float32).

### 6.2 Candidates

Attack candidates (transition from block i to i+1, i = 0..31): with a running maximum M,
initialised to prev.peak and updated M = max(M, blk[i]) before each test, block i+1 is a
candidate if `blk[i+1] > 10` and `blk[i+1] > 1.5·M` (double). Its ratio is
`r = blk[i+1]/M`, or `blk[i+1]/4` when M ≤ 4.

Release candidates (decays, g = 7 down to 0): with M2 initialised to max(blk[63], grp[8..15]) and
updated M2 = max(M2, grp[g]) before each test, let v = grp[g−1] (or `pre` for g = 0). It is a
candidate if `v > 10` and `v > 1.85f·M2` (1.85 rounded to float32, product in double). Ratio
`v/M2`, or `v/4` when M2 ≤ 4.

Candidate ratios are stored as float32 in a 40-entry array (attacks at index i, releases at
32 + g). The **7 largest** are selected (repeatedly take the first maximum, then clear it). Only
selected candidates are used below.

### 6.3 Levels

The level of a selected candidate is `a = trunc(log2(r') + 0.5)` with r' = max(r, 1e−20) and
r computed in double (not the float32 copy). (Sony evaluates log2 as
log10(r)/log10(2) with `fyl2x`.)

Attacks, scanned in increasing position: tot starts at prev_lev. For each selected attack with
a ≥ 0: if a + tot > 10 then a = 10 − tot; tot += a; record (position i, a) if a > 0. Stop after 7
recorded attacks.

Releases are processed only if at least one attack was recorded or the previous frame had gain
points. The cap is min(tot, 4). Scanned g = 7 down to 0, for each selected release with a ≥ 0:
clamp a so that the running release sum stays ≤ cap; record (position, a) if a > 0 and fewer than
7 attacks exist, where position = 4g (1 for g = 0). Stop when attacks + releases reach 7.

Accumulation: release values become suffix sums in recording order, each capped at 4 (the last
recorded, lowest position, keeps its value). With R = the largest of them (0 if none), attack
values become suffix sums in recording order capped at 11 + R.

Log-gain profile over the 33 positions 0..32, using the **already accumulated** values:

```
G[0..32] = 0
j = 0
for (pos, value) in attacks in recording order:       # increasing positions
    while j <= pos:
        G[j] += value
        j += 1
j = 32
for (pos, value) in releases in recording order:      # decreasing positions
    while j >= pos:
        G[j] += value
        j -= 1
```

Each list fills disjoint intervals; attack and release contributions are added where their
intervals overlap. Do not sum all accumulated values whose positions lie on one side of j:
that would count the suffix sums more than once. For example, accumulated attacks
`[(3, 5), (7, 2)]` and no releases give G[0..3] = 5, G[4..7] = 2, G[8..32] = 0.
Further profile and gain-point examples are in `vectors/gain-profiles.json`.

Gain points: for j = 0..31 in order, wherever G[j] ≠ G[j+1], emit (level G[j] − G[32] + 4,
location j). Then drop trailing points whose level is 4. More than 7 points is an error (never
occurs).

### 6.4 Band-0 adjustment (`adjust_band0_gain`)

For each band b with points, let *min* be the smallest level (initial value 4, first occurrence
wins) and *max* the largest level among the points up to and including that minimum (initial 0);
spread = max − min (0 for a band without points). If spread(band 1) ≥ 2, band 0 has no points in
this and the previous frame, the previous frame's band-0 peak ≤ 16384, and no sample among the
first 8·(loc + 1) samples of band 0's current block exceeds 16384 in magnitude (loc = location of
band 1's first point), then band 0 gets the single point (level 5, loc).

## 7. Gain window (`gain_window`)

Bands with gain points in this or the previous frame are divided by a 512-sample window before the
MDCT. Build 64 exponents E[0..63] (8-sample blocks of [old | current]), all 0 initially:

1. Current frame's points in order: point (L, loc) sets E[p .. loc+32] = L − 4, continuing from
   where the previous point stopped (the first point therefore also covers the whole old half).
2. Previous frame's points in order: point (L, loc) adds L − 4 to E[p .. loc], continuing likewise.

Then the window w[0..511] (float32), built from the end with `last = 0`:

```
for k = 63 down to 0:
    e = −E[k]
    if e == last:  w[8k .. 8k+7] = 2^e            (exact power of two)
    else:          w[8k+7−m] = f32(2^(((m+1)·e + (7−m)·last)/8)),  m = 0..7
    last = e
```

The MDCT input is `f32(f64(old[i]) / w[i])` and `f32(f64(current[i]) / w[256+i])`. Bands without
gain points in both frames use the samples unchanged.

## 8. MDCT (`mdct256`)

256 coefficients from 512 samples x, computed exactly as follows (all intermediate arrays float32
unless stated; `w[n] = MDCT_WIN[n]·x[n]` in double):

1. Folding into z[0..255]:
   z[n] = (−w[384+2n]) − w[383−2n] for n < 64;
   z[64+m] = w[2m] − w[255−2m] for m < 128;
   z[192+n] = w[256+2n] + w[511−2n] for n < 64.
2. Pre-rotation (k = 0..127): re[k] = c_k z[2k] − s_k z[2k+1], im[k] = s_k z[2k] + c_k z[2k+1]
   with c_k = cos(2πk/256), s_k = −sin(2πk/256).
3. 128-point complex FFT: 7-bit bit-reversal permutation, then 7 radix-2 stages
   (span 2, 4, …, 128; twiddle for butterfly b of a stage with *g* groups: cos/−sin(2π·b·g/128)).
   Each butterfly computes t = tw·bottom in double and stores bottom = f32(top − t),
   top = f32(t + top).
4. Post-rotation (k = 0..127, θ = π(2k+1)/1024, p1..p4 below, R/I the FFT output):
   y[k] = ((R[127−k]·p2 + I[k]·p3) + p1·R[k]) + p4·I[127−k],
   y[255−k] = ((p3·R[k] − p4·R[127−k]) − p1·I[k]) + p2·I[127−k],
   p1 = (cos θ − sin 5θ)/2, p2 = (cos θ + sin 5θ)/2, p3 = (sin θ + cos 5θ)/2, p4 = (cos 5θ − sin θ)/2.
5. out = y/128; for odd QMF bands the 256 outputs are reversed.

MDCT_WIN[n] = sin²(π(n+0.5)/512). All tables are float32 roundings of these double formulas,
with two quirks reproduced by the reference model: −sin(0) is stored as +0.0, and cos(π/2) is the
float32 with bit pattern 0x248D3000 instead of the correctly rounded value. Mathematically,
`out[k] = (1/128)·Σₙ x[n]·MDCT_WIN[n]·cos(π/256·(n + 0.5 + 128)·(k + 0.5))` (checked numerically
to float32 precision), i.e. the transform whose inverse is the decoder's IMDCT ([ATRAC3.md](ATRAC3.md) §5.1);
the FFT form above is only needed for bit-exactness.

Spectrum **A** = MDCT of the gain-compensated input. Spectrum **B** = MDCT of the same blocks
without gain compensation, computed with its own overlap block only when any band has gain points in
this or the previous frame; otherwise B = A.

## 9. Scale factors, quantiser, bit counts

* SF[i] = f32(2^((i−15)/3)), i = 0..63 ([Appendix A](ATRAC3-tables.md)).
* **Scale-factor index of a magnitude m** (float32): the smallest i with SF[i] > m, or 63 if none.
  Sony computes it with a 5-step binary search (`sf_index`: start 32, steps
  16, 8, 4, 2, 1, up when SF[e] ≤ m, then +1 if e < 63 and SF[e] ≤ m; 0 if m < SF[0]) and, for
  tonal components, with a linear search (`sf_index_tone`); both give exactly this result
  (checked exhaustively at every table boundary and on 20 000 random values).
* **Group index** (`group_sf`): the index of max |x| over each 4-line group g = 0..255.
* **Subband index** (`subband_sf`): maximum of the group indices in the subband.
  SA[i] is computed from A's groups, SB[i] from B's groups.
* **Word length** W (selector 0..7), maxq = {0,1,2,3,4,7,15,31}.
* **Dead zone** `dead_zone(e, W)`: 0 for e = 0, otherwise f32(T1[min(e,12)] / T2[W]) with
  T1[k] = 10^(k/40), T2 = 2·maxq+1.
* **Quantiser** `quantise_subband` on normalised lines x (float32):
  q = clip(trunc(f32(maxq+0.5)·x + 31.5) − 31, −maxq, maxq); q = 0 where |x| < dead zone.
  Cost = mantissa code bits (VLC table of W, [Appendix A](ATRAC3-tables.md); selector 1 codes pairs) + 6 for the SF.
  W = 0 costs nothing.

## 10. Sound-unit coding (`encode_unit`)

Division of labour between the two spectra: the **plain** spectrum B (no gain window) drives the
perceptual decisions (peak position, mean level, multitone detection, the SB profile behind the
offsets and the masking rule); the **gain-compensated** spectrum A is what gets coded (single-tone
check, transmitted scale factors SA, normalisation, quantisation). Without gain points they are
identical.

### 10.1 Peak and tonal class

`GA`, `GB` = group indices of A and B. `peak` = first group with the largest GB.
`single_tone_class(GA)` (class 0/1/2): let p be the first group with the largest GA, r = p mod 64,
its band b0 = p div 64 and the neighbour band b1 = b0−1 if r < 32 else b0+1 (where the QMF alias
image of a tone near the band edge appears, mirrored at position 63−r). mx = the largest GA over
all 256 groups except the groups r−2..r+2 of band b0 and (63−r)−2..(63−r)+2 of band b1 (positions
within the band; a band outside 0..3 excludes nothing). With v = GA[p]: if mx < max(v − 17, 9) the
class is 2 when v ≥ 36, 1 when v ≥ 21, else 0. (17 index steps are 34 dB.)

### 10.2 Bit budget

```
gain_bits(nb) = 3·(nb + 3·(number of gain points in bands 0..nb−1))
budget = 1536 − 3·nsub − 8 − gain_bits(nbands) − 13
```
(8 = unit id and band count; 13 = tonal group count 5 + tonal mode 2 + subband count 5 +
coding-mode bit 1; 3 per subband for the selectors.)

### 10.3 Single-tone path (`single_tone_path`)

Taken when class > 0 and nsub > 25 (always true at entry). If peak ≥ 160 the frame uses nsub = 32,
nbands = 4; else if peak ≥ 96, nsub = 30, nbands = 3 (else 29/3). The budget is computed with these
values. Then `class` tone groups are created (selector 7, length code 3 = 4 values, CLC), costing
(nbands + 6) bits each up front. Bands of interest: the peak's band b0 (centre r) and the
neighbour b1 (centre 63−r); width = 2 groups (1 if the unit had ≤ 1200 bits).

```
bits = (nbands+6)·class; go_on = true
for each of the `class` groups:
    if not go_on: leave the group empty
    for e in (b0, centre r), (b1, centre 63−r), skipping bands outside 0..nbands−1:
        mark the band in the group (bits += 12 for the four 3-bit cell counts)
        if bits ≥ budget: unmark (bits −= 12); continue with the next band
        for g in max(c−width, 0) .. min(c+width, 255), c = 64·band + centre:
            cell = g div 16 (64-line cell); skip if the cell already has 7 components
            quantise a component at lines 4g..4g+3 (§11) → cost
            if bits + cost > budget: go_on = false; continue
            accept: append to the cell, subtract it from A, GA[g] = group index of the residual,
                    bits += cost
```
The spectrum budget is reduced by `bits`. The masking slope used later is 60 (instead of 30), and
the multitone search is skipped.

### 10.4 Mean and attack state

`mean = f32(Σ GB[0 .. G−1] / G)` with G = start[nsub]/4 groups.
Attack state: 0 if any coded band has `is_attack` points now (levels including 4 span > 2),
else 1 if any had them in the previous frame, else −1.

### 10.5 Multitone path (`multitone_path`)

Only when the single-tone path produced no bits. Threshold t = f32(mean + 32). Up to two passes:

```
pass: count candidates = groups g < G with GB[g] > t, at most 7 per cell and 64 components in total
      if none: next pass;  if bits + nbands + 6 > budget: stop all passes
      bits += nbands + 6; new group (selector 7, length code 3, CLC)
      for g in 0..G−1 with GB[g] > t, cell not full (< 7), fewer than 64 components:
          band = g div 64; if the band is not yet marked: mark it if bits + 12 ≤ budget
                           (bits += 12), else set stop
          if marked: if bits + 36 > budget: set stop (no component)
                     else: add component at 4g, bits += its cost (§11), subtract it from A,
                           s = group index of the residual; GB[g] += s − GA[g]; GA[g] = s
```
(36 = 4 values × 6 CLC bits + 12.) "stop" only prevents a further pass. The budget is reduced by
`bits`. If the unit ends up with no tone groups at all, the budget gains 2 bits (no tonal mode
field).

### 10.6 Subband scale factors and initial word lengths

SA, SB = subband indices from the (updated) GA, GB. Offsets per subband:

`off[i] = f32((SB[i] − mean) · c_i)` with c_i = 1/3 (i<8), 1/3.3 (<12), 1/3.4 (<16), 1/3.5 (<18),
1/3.6 (<26), 1/3.8 (<28), 1/4 (<30), 1/4.2 (as double constants).

If the attack state is not −1: off[0..7] += 0.7, off[8..17] += 0.5 (double add, float32 store),
then off[0] = max(off[0], 6). Otherwise off[0] = max(off[0], 6) and off[1..3] = max(·, 3).

`word_lengths(slope, off)`:
1. W[i] = clamp(trunc(off[i] + 0.5), 1, 7).
2. quiet = 6 if the largest SB[1..nsub−1] ≥ 30, else max(1, that maximum div 6).
3. W[i] = 0 if SB[i] < quiet, or if SB[i] < SB[j] − slope for any of the 8 subbands j = i−1 … i−8
   (slope 30, or 60 on the single-tone path).

After the first call: **fixed[i]** = 1 when SA[i] < quiet (and then W[i] = 0); when the multitone
path produced no bits, also fixed[i] = 1 when W[i] = 0 or SA[i] < mean. W0 = copy of W (the upper
bound for later increases).

### 10.7 Normalisation and dead-zone limit

Every subband of A is multiplied by f32(1 / SF[SA[i]]) (product in double, stored float32).
Dead-zone limit from the density d = budget / (start[nsub] + 256) (double): 0 if d > 1.15; then
with d' = f32(d): 3 if d' > 1.1, 5 if > 1.05, 7 if > 1.0, 9 if > 0.95, 11 if > 0.9, else 12.

### 10.8 Offset search

E[i] = 0 (dead-zone index). Quantise all subbands; *total* = Σ cost.

```
offset = 2.0, step = 4.0 (float32)
for it = 0..14:
    if offset > 0:  off2[i] = f32(offset + off[i])
    else:           off2[i] = f32(off[i] + offset·w_i), w_i = 0.2 (i=0), 0.3 (i=1), 0.4 (i<8),
                    0.6 (i<18), else f32(offset + off[i])
    W = word_lengths(slope, off2), with W[i] = 0 wherever W0[i] = 0
    re-quantise the subbands whose W changed; update total
    if total > budget:          offset −= step
    elif it ≥ 6:                stop
    elif total < budget:        offset += step
    else:                       offset −= step          (exactly on budget)
    step ×= 0.5 if it < 7 else 1.5
```

### 10.9 Forced fitting

Four passes over i = nsub−1 down to 0 while total > budget: if W[i] > 0 and E[i] < limit,
E[i] += 1 and re-quantise (dead zone grows).

If still total > budget: order the subbands with `order_descending` (Shell sort, gaps …13, 4, 1,
descending; ties resolved by that algorithm) on the key `SB[i] − i div 2` (normal) or
`(SB[i]+1)·32 − i` (single-tone path). Starting from the **last** subband of that order, lower
W[i] one step at a time (re-quantising) until W[i] = 0 or total ≤ budget, then move to the
previous subband of the order.

The exact `order_descending(keys)` procedure is below. Integer division is floor division;
the returned permutation, including equal-key order, is normative. This is not a stable sort.

```
n = length(keys); K = copy(keys); order = [0, 1, ..., n-1]
h = 1
while h <= n: h = 3*h + 1
repeat:
    h = h div 3
    if h < 1: return order
    if h >= n: continue
    for t = 0 .. n-h-1:
        v = K[h+t]; j = t
        while j >= 0 and K[j] < v:                   # strict comparison
            K[j+h] = K[j]
            swap(order[j+h], order[j])
            j -= h
        K[j+h] = v
```

In particular, the index array is updated by swaps, not by a separate stable insertion.
Use this same procedure in §10.11(d). Equal margins in §10.10 are visited in descending
subband index within each threshold class; equal candidate ratios in §6.2 use the first
array index. Golden equal-key examples are in `vectors/ordering.json`.

### 10.10 Coded subband count and the remaining bits

nused = index of the last subband with W > 0, plus one (at least 1). The subband selectors are
W[0..nused−1], the scale factors SA[0..nused−1]. All unit bits so far:
`used = 1536 − (budget + 3·(nsub − nused) − total)`.

**Margin order.** margin[i] = f32(off2[i] − W[i]) for W[i] > 0, −9.99 otherwise (off2 from the
last search iteration). With lo/hi the minimum/maximum margin, for k = 1..10 the threshold is
f32(hi − k·(hi − lo)·0.1); scanning i = nused−1 down to 0, every subband not yet listed whose margin
is ≥ the threshold is appended. (Highest margins first; within a class, higher subbands first.)

### 10.11 Refinement

All steps re-quantise one subband and accept the change only if `used` stays ≤ 1536.

a. For each subband in margin order that is not fixed and has 0 < W < 7: if W < W0 raise W by one;
   else if E > 0 lower E by one; else (nothing to do) stop the whole pass if used > 1528. After an
   accepted change, stop the pass if used > 1528. A rejected change is undone.
b. For every subband with E > 0 and W > 0: re-quantise with E = 0; keep if it fits.
c. Up to 7 passes while used ≤ 1534: for each subband (ascending) with 0 < W < W0 and not fixed,
   raise W by one; if it does not fit, undo and mark the subband fixed.
d. If nused > 1: order the subbands with `order_descending` on `(SB[i]+1)·32 − i`; for p = 0..3
   (while used ≤ 1534): for each subband in that order with 0 < W < p+3, try W+1.
e. Up to 7 passes while used ≤ 1534: for each subband with 0 < W < 7, try W+1.

(The "used ≤ 1534" condition is re-checked before every subband.)

### 10.12 Scale-factor correction

For each coded subband: e_orig = Σ x² over the normalised lines, e_quant = Σ (q·s)² with
s = f32(1/(maxq + 0.5)), both accumulated in double in line order. If 1.25·e_orig < e_quant and
the SF index > 0, the transmitted SF index is decreased by one (the mantissas are kept). The bit
count is unaffected.

### 10.13 Result

nbands, gain points of bands 0..nbands−1, the tone groups, nused, coding mode 0 (VLC), the
selectors, SF indices and mantissas. The tonal coding-mode field (present when there are groups)
is 1 (all components CLC). The total size equals `used`, which never exceeds 1536.

## 11. Tonal components (`quantise_component`, `tone_mantissa`, `subtract_component`)

A component covers 4 lines starting at a multiple of 4 (selector 7, maxq 31, CLC).
* SF index: the scale-factor index (§9) of max|x| over the 4 lines of the current (partly cleaned)
  spectrum A.
* Mantissa: clip(trunc((maxq + 0.5)·(x / SF) + 31.5) − 31, ±31) with x/SF and the product in double.
* Cost: 12 (SF + position) + CLC bits (6 per value) = 36.
* Subtraction: x = f32(x − m·f32(SF/(maxq + 0.5))).

Components are stored per tone group and per 64-line cell (max 7 per cell, 64 per unit, 31 groups);
the bitstream order is group, band, cell ([ATRAC3.md](ATRAC3.md) §4.4). Position within the cell = line mod 64.

## 12. Serialisation (`pack_unit`)

Exactly the syntax of ATRAC3.md §4.2 (MSB first, zero padding): id 0x28 (6), nbands−1 (2); per
band: point count (3), then (level 4, location 5) per point; group count (5) and, if > 0, tonal
mode 1 (2); per group: band flags, length code 3 (3 bits), selector 7 (3), per marked band four
cell counts (3) each followed by its components (SF 6, position 6, 4 CLC values of 6 bits);
nsub−1 (5); coding mode 0 (1); selectors (3 each); SF indices (6) of subbands with W > 0; mantissa
codes (VLC; selector 1 codes pairs).

## 13. State carried between frames (per channel)

* QMF delay lines (3 × 46 samples).
* Band histories (4 × 3 × 256 samples).
* Two overlap blocks per band (gain-compensated and plain MDCT; both hold the raw "current" block).
* Previous frame's gain points and band peaks.
* Library: sample counts of the last three calls (end flag).

Nothing else survives a frame; in particular all bit-allocation state is rebuilt from scratch.
At the start everything is zero (the first call sees a zero history, zero gain data and peak 0).

## 14. Edge cases

* **Start:** the tool's 955 leading zeros, the zero history and the two discarded outputs (§2.3).
* **End:** flushing with zero calls until the end flag (§2.2, §2.3); the final partial block is
  zero-padded.
* **Silence:** all group indices 0; every subband falls below the quiet threshold (1), W = 0
  everywhere, and the unit codes one subband with selector 0: 31 bits (id, band count, three empty
  gain lists, no tone groups, subband count, coding mode, one selector), the rest zero.
* **Clipping / full scale:** handled by the normal path (tested with full-scale square waves and
  clipped noise). The band-0 gain adjustment is suppressed above |x| = 16384 (§6.4).
* **Never reached in this configuration [inferred]:** a frame with zero subbands (the code would
  emit one subband with W = 0); units of ≤ 20 bytes (the code would emit a minimal header); the
  state flag that suppresses packing; CLC spectrum coding; tonal mode 3; a component clipped at
  line 1024; internal error returns.

## 15. Tables (addresses in psp_at3tool.exe)

Exact values are exported in [`tables/lp2.json`](tables/lp2.json), with representation and
name mappings in [`tables/README.md`](tables/README.md). Use the exported IEEE bit patterns
for byte-exact work; this avoids host-math differences when regenerating the formulas below.

| table | address | definition |
|---|---|---|
| subband start / size | 0x48E980 / 0x48E900 | [Appendix A](ATRAC3-tables.md) (32 subbands) |
| SF | 0x48EAE8 | f32(2^((i−15)/3)) |
| maxq per selector | 0x48EA24 | 0,1,2,3,4,7,15,31 |
| CLC bits per value | 0x48EA04 | 0,2,3,3,4,4,5,6 |
| values per codeword | 0x48ED40 | selector 1: 2, others 1 |
| bit estimate (VLC, CLC) | 0x48EBE8 | {0,3,3,4,5,6,7,8}, {0,2,3,3,4,4,5,6} |
| dead zone T1 / T2 | 0x48EC68 / 0x48EC9C | f32(10^(k/40)), k=0..12 / 1,3,5,7,9,15,31,63 |
| masking slope | 0x48EC28 / 0x48EC48 | 8 × 30 / 8 × 60 |
| tone length code → values | 0x48EA88 | code+1 |
| cell → band | 0x48EA48 | cell div 4 |
| gain level → log2 | 0x48EAA8 | level − 4 |
| code tables VLC / CLC | 0x48EDC0 / 0x48F2C0 | (code, length) pairs, same codes as [Appendix A](ATRAC3-tables.md) |
| QMF prototype | 0x490FC0 | [Appendix A](ATRAC3-tables.md) |
| MDCT window | 0x4907C0 | f32(sin²(π(n+0.5)/512)) |
| pre-rotation cos / −sin | 0x48F7C0 / 0x48F9C0 | 2πk/256, k<128 |
| FFT twiddles cos / −sin | 0x4903C0 / 0x4904C0 | 2πk/128, k<64 |
| bit reversal | 0x4905C0 | 7-bit |
| post-rotation p1..p4 | 0x48FBC0, 0x48FDC0, 0x48FFC0, 0x4901C0 | §8 |
| offset factors | 0x491120 … 0x4910F0 | §10.6 |
| density thresholds | 0x4910D0 … 0x4910B0 | 1.15, 1.1, 1.05, 1.0, 0.95, 0.9 |
| negative-offset weights | 0x4910A8 … 0x491090 | 0.2, 0.3, 0.4, 0.6 |

## 16. Verification

Method: the encoder runs under the CPU emulator (`tools/at3emu.py`); each stage was compared with
the model at function boundaries (QMF outputs, gain points, MDCT outputs, spectra A/B, the internal
state of the allocator at ten points), then whole units and whole files. Script:
`tools/verify_lp2.py [--mono] [--extra] [--music]`.

| test material | frames | identical |
|---|---|---|
| `gen_signals.py` set (mix, transients, tones, corr, noise, silence_burst), 132 kbps | 786 | 786 (100 %), files byte-identical |
| same, 66 kbps mono (left channel) | 786 | 786 |
| edge cases (`--extra`): full-scale square/clipped noise, ±1 LSB, clicks, DC steps, sweep, weak tones, a length needing three flush calls; stereo and mono | 1100 | 1100 |
| PSP SDK music bgm01_1, bgm01_2, overload_short (PCM only), stereo and mono | 7108 | 7108 |

All float tables of the model are bit-identical to the binary. Over this material every branch of
§6–§11 listed above was exercised (both tonal classes, the 30/32-subband switches, multitone,
all dead-zone limits, attack states, word-length reduction, every refinement pass, the SF
correction, band-0 gain adjustment, attacks and releases), except those marked [inferred] in §14.

One numerical caveat: the level rounding in §6.3 uses `fyl2x`; the model (like the emulator)
evaluates log2 with the host's `log`. Only a ratio within an ulp of 2^(k+0.5) could round
differently on a real x87; no such case is known.
