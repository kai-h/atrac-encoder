# Exact encoder tables

`common.json`, `lp2.json`, and `lp4.json` are implementation-neutral, exact table data for
the specifications. They are part of the specification: for byte-exact work use these bit
patterns rather than re-evaluating formulas or abbreviated decimal constants. No reference
Python source or Sony binary is needed to read them.

Each file has `format_version: 1` and a `tables` object keyed by table name.

| Document name | Export |
|---|---|
| LP2 subband boundaries / widths | `lp2.json`: `SB_START`, `SB_SIZE` |
| LP2 SF / dead-zone tables | `SF`, `T1`, `T2` |
| LP2 QMF / MDCT | `QMF_WIN`, `MDCT_WIN`, `PRE_COS`, `PRE_SIN`, `FFT_COS`, `FFT_SIN`, `BITREV7`, `POST1`–`POST4` |
| LP2 allocation constants | `OFFSET_FACTORS`, `NEGATIVE_OFFSET_WEIGHTS`, `DENSITY_THRESHOLDS`, `DOUBLE_CONSTANTS`, `FLOAT_CONSTANTS` |
| LP2 codebooks | `ENC_CODES`, `CODE_LEN` |
| LP4 Q1 / Q2 | `lp4.json`: `QMF1`, `QMF2` |
| LP4 PERM / R0–R3 / C | `MDCT_PERM`, `MDCT_TW` (128 rows, columns R0–R3), `MDCT_FFT_TW` |
| LP4 gain / stereo / quantiser | `GAIN_POW2`, `GAIN_RAMP_BITS`, `CHSW_TABLE`, `CHWT_TABLE`, `SF_BITS`, `QSCALE_BITS`, `QMASK` |
| LP4 allocation / codes | `MIN_SFSUM`, `ZERO_LEN2`, `STEP_COST`, `EST_SAVE`, `EST_BITS`, `EST_THR`, `TONE_PENALTY`, `TONE_SF_OFFSET`, `SEL1_LEN`, `VLC` |
| Shared / decoder data | `common.json`: `CODES`, `MAX_QUANT`, `SUBBAND_START`, `SF_TABLE`, `QMF_HALF_SONY` |

Remaining named constants and mode tables are included too. Extra LP4 mode-table rows do not
extend the supported encoder scope. Unused entries are preserved, not recommendations to
encode invalid mantissas. LP2 `NEGATIVE_OFFSET_WEIGHTS` uses 1.0 for unweighted subbands 18–31.

## Representation

- `{"type":"binary32", "shape":[...], "bits":[...]}`: hexadecimal IEEE-754 single-precision
  bit patterns without a `0x` prefix. Reinterpret as floating point; do not numerically convert
  the integer to a float. `binary64` uses 16 hexadecimal digits for double precision.
- `shape: []` means a scalar with one bit pattern. Arrays are flattened in row-major order;
  the last dimension varies fastest. The product of dimensions equals the number of entries.
- `{"type":"integer", "shape":[...], "values":[...]}` holds mathematical integers with
  the same flattening convention. Plain JSON integers represent scalar integers.
- `{"type":"map", "entries":[{"key":..., "value":...}]}` preserves codebook keys.
  A key `[mode, selector]` represents a pair. Nested decoder maps use string keys containing
  the actual MSB-first codewords.
- Named groups such as `DOUBLE_CONSTANTS` are ordinary JSON objects containing typed values.

Encoder codebook rows are `(code, bit_length)`; emit the low `bit_length` bits of `code`, MSB
first. LP2 rows use Chapter 7's masked-mantissa indexing and `QBITS`/`VALS_PER_CODE`.
LP4 rows use Chapter 8's `QMASK` rule. Tables ending in `_BITS` are integer arrays of binary32
bit patterns stored as JSON integers; reinterpret them when the spec calls for a float value.
Byte order matters only when serialising a pattern to memory, not when reading its hex string.

## Provenance and regeneration

The files export existing reference-model constants. LP4 and shared source tables were
extracted from `psp_at3tool.exe`; LP2 source tables come from the documented formulas and quirks
and were previously compared with that binary (Chapter 7 §16). This export does not independently
re-establish those historical comparisons. Frozen JSON values are authoritative for table
values; formulas explain their construction.

`python spec/generate_assets.py --write` regenerates tables and vectors; `--check` recomputes
them and compares bytes, detecting model drift or platform-dependent formula changes. Both
need NumPy only. `python spec/verify_assets.py` checks frozen hashes using the standard library
alone. SHA-256 hashes and byte sizes are in `../assets-manifest.json`.
