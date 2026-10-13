"""Constants and tables of Sony's ATRAC3 encoder algorithm 1 (MiniDisc LP2 / 66 kbps mono).

Every table is given either as a closed formula or as a short integer list. Comments name the
address of the identical table in psp_at3tool.exe; tools/verify_lp2.py checks each one against
the binary (all float32 tables below are bit-identical to Sony's).

Float tables are float32, as in the reference.  Where the reference evaluates an expression in
x87 double precision, the companion modules convert to float64 explicitly.
"""
import numpy as np
from at3_tables import CODES as _DEC_CODES

f32, f64 = np.float32, np.float64

# ---------------------------------------------------------------- configuration (LP2)
SAMPLE_RATE = 44100
FRAME_SAMPLES = 1024           # samples per channel per frame
UNIT_BYTES = 192               # one sound unit (one channel)
UNIT_BITS = UNIT_BYTES * 8     # 1536-bit budget per channel per frame
BANDWIDTH_HZ = 15100           # 0x3afc, set by the encoder set-up
PRIME_ZEROS = 955              # zeros the tool feeds before the source (1024 - 69)

# ---------------------------------------------------------------- spectral layout
# 32 quantisation subbands ("qu"), edges in MDCT lines (0x48e980, 33 entries; sizes 0x48e900)
SB_START = [0, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160, 176, 192, 224, 256,
            288, 320, 352, 384, 416, 448, 480, 512, 576, 640, 704, 768, 896, 1024]
SB_SIZE = [SB_START[i + 1] - SB_START[i] for i in range(32)]


def coded_layout(bandwidth_hz=BANDWIDTH_HZ, rate=SAMPLE_RATE):
    """Number of subbands and QMF bands normally coded for a bandwidth."""
    lim = (bandwidth_hz << 11) // rate                    # 701 MDCT lines for 15100 Hz
    k = 0
    while k < 32:
        if lim <= SB_START[k + 1]:
            break
        k += 1
    nsub = k + 1                                          # 29 (lines 0..703)
    nbands = (lim + 255) >> 8                             # 3 QMF bands
    return nsub, nbands


NSUB_DEFAULT, NBANDS_DEFAULT = coded_layout()

# scale factors (0x48eae8): SF[i] = 2^((i-15)/3) rounded to float32
SF = np.array([2.0 ** ((i - 15) / 3.0) for i in range(64)], dtype=f32)

# word length (quantiser selector, "idwl") 0..7
MAXQ = [0, 1, 2, 3, 4, 7, 15, 31]             # largest |mantissa| (0x48ea24)
QBITS = [0, 2, 3, 3, 4, 4, 5, 6]              # fixed-length bits per value (0x48ea04)
VALS_PER_CODE = [0, 2, 1, 1, 1, 1, 1, 1]      # selector 1 codes value pairs (0x48ed40)
EST_BITS_CLC = [0, 2, 3, 3, 4, 4, 5, 6]       # per-value estimate, CLC column of 0x48ebe8
T2 = [1.0, 3.0, 5.0, 7.0, 9.0, 15.0, 31.0, 63.0]                # 2*MAXQ+1 (0x48ec9c)
T1 = np.array([10.0 ** (k / 40.0) for k in range(13)], f32)     # dead-zone steps (0x48ec68)

# tonal components
TONE_NVALS = [1, 2, 3, 4, 5, 6, 7, 8]          # length code -> coded values (0x48ea88)
TONE_SEL = 7                                   # quantiser selector used for every component
TONE_LEN_CODE = 3                              # -> 4 values per component
TONE_MODE = 1                                  # components are CLC coded
CELL_BAND = [c >> 2 for c in range(16)]        # 64-line cell -> QMF band (0x48ea48)

# masking slope used by word_lengths (0x48ec28 / 0x48ec48)
MASK_SLOPE_NORMAL = [30] * 8
MASK_SLOPE_TONAL = [60] * 8

# offset weights per subband range (doubles at 0x491120 ... 0x4910f0)
def offset_divisor(i):
    """(SB_i - mean) is multiplied by 1/d with d from this list (applied as a double constant)."""
    if i < 8: return 0.3333333333333333          # 1/3
    if i < 12: return 0.30303030303030304        # 1/3.3
    if i < 16: return 0.29411764705882354        # 1/3.4
    if i < 18: return 0.2857142857142857         # 1/3.5
    if i < 26: return 0.2777777777777778         # 1/3.6
    if i < 28: return 0.2631578947368421         # 1/3.8
    if i < 30: return 0.25                       # 1/4
    return 0.23809523809523808                   # 1/4.2


def offset_step_weight(i):
    """Weight of a negative global offset per subband (0x4910a8, 0x4910a0, 0x491098, 0x491090)."""
    if i < 1: return 0.2
    if i < 2: return 0.3
    if i < 8: return 0.4
    if i < 18: return 0.6
    return None                                  # added unweighted


def dead_zone_limit(density):
    """Maximum dead-zone index from bit density (thresholds 0x4910d0 ... 0x4910b0).
    `density` is the double budget/(lines+256); the first test uses the double, the rest the
    float32-rounded value, as in the reference."""
    if not (density <= 1.15):
        return 0
    d = f64(f32(density))
    if d > 1.1: return 3
    if d > 1.05: return 5
    if d > 1.0: return 7
    if d > 0.95: return 9
    if d > 0.9: return 11
    return 12

# ---------------------------------------------------------------- gain control
GAIN_LEVEL_EXP = [lev - 4 for lev in range(16)]   # level code -> log2 step (0x48eaa8)

# ---------------------------------------------------------------- analysis filter bank
# 48-tap QMF prototype (0x490fc0); identical to the decoder's (spec Appendix A)
_QMF_HALF = [-0.00001461907, -0.00009205479, -0.000056157569, 0.00030117269, 0.0002422519,
             -0.00085293897, -0.0005205574, 0.0020340169, 0.00078333891, -0.0042153862,
             -0.00075614988, 0.0078402944, -0.000061169922, -0.01344162, 0.0024626821,
             0.021736089, -0.007801671, -0.034090221, 0.01880949, 0.054326009,
             -0.043596379, -0.099384367, 0.13207909, 0.46424159]
QMF_WIN = np.array(_QMF_HALF + _QMF_HALF[::-1], f32)


def _twiddle(k, n):
    """float32 cos(2*pi*k/n) and -sin(2*pi*k/n), evaluated in double.  Two details of Sony's
    tables are reproduced: -sin(0) is stored as +0.0, and cos(pi/2) is not 0 but the float32
    with bit pattern 0x248d3000 (6.1230318e-17; the double formula would give 0x248d3132)."""
    a = 2.0 * np.pi * np.asarray(k, f64) / n
    c = np.cos(a).astype(f32)
    c[np.asarray(k) * 4 == n] = np.array([0x248d3000], np.uint32).view(f32)[0]
    return c, ((-np.sin(a)) + 0.0).astype(f32)


_n = np.arange(512)
MDCT_WIN = (np.sin(np.pi * (_n + 0.5) / 512.0) ** 2).astype(f32)          # 0x4907c0
PRE_COS, PRE_SIN = _twiddle(np.arange(128), 256)                         # 0x48f7c0 / 0x48f9c0
FFT_COS, FFT_SIN = _twiddle(np.arange(64), 128)                           # 0x4903c0 / 0x4904c0
BITREV7 = [int(format(i, '07b')[::-1], 2) for i in range(128)]           # 0x4905c0
_t = np.pi * (2 * np.arange(128) + 1) / 1024.0
POST1 = ((np.cos(_t) - np.sin(5 * _t)) / 2).astype(f32)                  # 0x48fbc0
POST2 = ((np.cos(_t) + np.sin(5 * _t)) / 2).astype(f32)                  # 0x48fdc0
POST3 = ((np.sin(_t) + np.cos(5 * _t)) / 2).astype(f32)                  # 0x48ffc0
POST4 = ((np.cos(5 * _t) - np.sin(_t)) / 2).astype(f32)                  # 0x4901c0

# ---------------------------------------------------------------- code tables
# Same codes as the decoder (spec Appendix A); here indexed by the two's-complement value masked
# to QBITS[sel] bits (selector 1: index = (m0 << 2) | m1).  ENC_CODES[(mode, sel)][idx] = (code, len)


def _build_codes():
    out = {}
    for (mode, sel), dec in _DEC_CODES.items():
        b = QBITS[sel]; m = (1 << b) - 1
        n = (1 << b) ** VALS_PER_CODE[sel]
        tab = [(0, 0)] * n
        for bits, val in dec.items():
            if VALS_PER_CODE[sel] == 2:
                idx = ((val[0] & m) << b) | (val[1] & m)
            else:
                idx = val & m
            tab[idx] = (int(bits, 2), len(bits))
        out[(mode, sel)] = tab
    return out


ENC_CODES = _build_codes()
CODE_LEN = {k: np.array([l for c, l in v], np.int64) for k, v in ENC_CODES.items()}
