# LP2 and LP4 golden vectors

Fixtures are synthetic and contain no SDK music or Sony executable. Expected encoded files
and checkpoints were produced by `spec/reference`. They are frozen regression targets, not
an independent proof of Sony equivalence. If `sony-verification.json` is present, it identifies
complete files subsequently compared byte for byte with the emulated Sony tool; checkpoints
remain model-derived.

## Encoding vectors

`cases.json` lists twelve stereo inputs and both expected modes for each: LP2 at 132 kbps with
384-byte independent-channel frames, and LP4 at 66 kbps with 192-byte joint-stereo frames.

Inputs are 16-bit little-endian, 44.1 kHz, two-channel WAVs. Encode the whole file with the
non-looping framing of ATRAC3.md §6.4, including the 955-sample lead-in and two discarded
analysis-call outputs. `NAME.MODE.at3` is the complete expected RIFF/WAVE file;
`NAME.MODE.frames` is just its `data` payload, concatenated fixed-size frames.
`cases.json` includes sample counts, frame counts, empty analysis-call counts and individual
frame hashes. `../assets-manifest.json` gives complete-file hashes and sizes, including inputs.

Cases cover silence, a correlated tone, weak tones and ±1 LSB input, transient noise and clicks,
left/right dominance and stereo transitions, clipping and broadband pressure, and final
source-plus-lead-in remainders 0, 885, 886 and 1023. Two cases separate readings that the
others do not: `tone_noise` (LP2) depends on Chapter 7 §10.6 fixing subbands on single-tone
frames, and `quiet_clicks` (LP4) on Chapter 8 §7.3 comparing the signed decision for gain
mode. Recipes use integer arithmetic, an explicit tone lookup table, and a specified 32-bit
LCG; their source is `generate_assets.py`.
These are compact regression vectors, not exhaustive branch coverage or a quality benchmark.
Mono, 105 kbps and PSP loop framing are outside this vector set.

## Intermediate checkpoints

`stages.json` uses **the first three 1024-sample blocks of `transients.wav` directly**, with
no lead-in or flush calls. State starts at zero as defined in each encoder chapter.
Floats use the bit-pattern format in `../tables/README.md`.

- `lp2`: channel 0; `qmf_bands` shape `[4,256]` in band order 0,1,2,3. Spectra A and B have
  shape `[1024]`, before tone extraction/allocation. Each gain entry has `(level, location)`
  pairs and a binary32 peak; all four bands and the previous gain state are included.
- `lp4`: two-unit lists `qmf_before_stereo` have shape `[256,4]` per unit (sample, band);
  `spectrum` shape `[4,256]` per unit, before allocation. `gain` has eight lists: unit 0
  bands 0–3 then unit 1 bands 0–3, with `(level, location)` pairs. Gain lists and budgets
  are recorded after packing, including empty-unit gain-count adjustments. `sent_matrix`
  and `sent_weight` are transmitted parameters; `budgets` are allocation budgets, not packed
  lengths. `frame_hex` is library output; call 0 is null because its full-call output is suppressed.

`ordering.json` holds four LP2 key arrays and expected permutations, including equal keys
at production sizes 29 and 32. Apply Chapter 7 §10.9's exact sort.

`gain-profiles.json` supplies hand-specified LP2 cases for Chapter 7 §6.3. Attack/release
pairs are `(position, already accumulated value)` in recording order, not raw step sizes.
Each case supplies the full 33-position profile and the resulting `(level, location)` points.

## Verification

From the repository root:

```sh
python spec/verify_assets.py                         # frozen hashes; standard library only
python spec/verify_assets.py --outputs /path/to/out  # compare twenty-four .at3 outputs
python spec/verify_assets.py --outputs /path/to/out --raw  # compare twenty-four .frames outputs
python spec/generate_assets.py --check              # recompute model outputs; needs NumPy
```

Candidate files must have the basenames in `cases.json`. The verifier reports the first
different byte and, for raw data, its frame and byte within the frame (zero-based). Missing,
shortened, lengthened or different files yield a nonzero exit status.

Maintainers with the private oracle can run `python spec/generate_assets.py --verify-sony`.
It checks stored asset hashes, compares complete files without changing expected outputs,
and writes `sony-verification.json` only after all comparisons pass. It needs the emulator's
additional dependencies and is not needed to implement or test an encoder against these fixtures.
