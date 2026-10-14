# Specification-driven ATRAC3 prototype

This is a separate Python/NumPy implementation developed from the Markdown specification,
`spec/tables/*.json`, and `spec/vectors/*`. It does not import or call the existing codec,
reference models, asset generator, emulator, or Sony executable. Existing implementation
source files are not consulted while developing this directory. NumPy supplies numeric arrays;
codec transforms and decisions are implemented here from the stated operation order.

Provenance limitation: the same chat inspected reference source during the preceding spec audit,
before the request for a clean-room implementation. This work follows a specification-only source
boundary from that request onward; it is not evidence of a separate, unexposed implementation team.

## Current implementation

- LP2 QMF, gain-point analysis, gain windows, both MDCT spectra, single-tone/multitone extraction, allocation,
  quantisation, refinement, packing and non-looping file framing.
- Shared sound-unit packing and independent/joint-stereo frame layout.
- LP4 QMF, joint stereo, gain detection/application, MDCT, bit split, tonal extraction,
  spectrum allocation, quantisation, packing and non-looping file framing.
- Exact JSON table reader, including binary32/binary64 patterns and indexed codebooks.

Both LP2 and LP4 match all twelve complete golden files per mode byte for byte: 24 files and
228 frames. All integrated stage fixtures also pass, including gain detection and frame packing
without substituting fixture gain points. `tone_noise` and `quiet_clicks` specifically check the
clarified fixed-subband and signed gain-mode rules. There are no expected failures or pending
specification questions in the conformance suite.

## Run

From the repository root, using a Python environment with NumPy:

```sh
python -m unittest cleanroom.test_conformance -v
python -m cleanroom.encoder spec/vectors/transients.wav /tmp/transients.lp2.at3
python -m cleanroom.encoder --mode lp4 spec/vectors/transients.wav /tmp/transients.lp4.at3
```

The CLI implements LP2 stereo, its independent-channel mono variant, and LP4 stereo at 44.1 kHz
with 16-bit PCM; short inputs and PSP loop metadata are rejected. The current golden set checks stereo.
This remains a standalone numerical prototype, not an integration into the main library or CLI.
Conformance is established for the published fixtures; exhaustive input coverage and production
performance have not been established. The code raises `SpecificationGap` if it encounters a
branch the specification marks unreachable but does not fully define for implementation.

## Clarified rules implemented

LP2 fixes subbands whenever the **multitone**
bit count is zero, including single-tone frames; LP4 compares absolute `sent` with signed `cur`
for gain mode; LP4 bandwidth limits are 27/12 and 28/28 as specified. The corrected side subband
limit gives the golden call-1 budget of 1149 using the unchanged bit-split formula.

LP4 all-tonal replacement reserves another `28 + 8·grp` bits without refunding its earlier cost.
It adds no tone to the count; a pass with replacements but no new tones ends iteration.
On `stereo.wav`, analysis call 10 (file frame 8), unit 1, there are 28 final tones and a reservation
of 1014 bits: 942 for the final records and headers plus two additional 36-bit reservations.
The reservation is subtracted from `left` without reclaiming the excess.
