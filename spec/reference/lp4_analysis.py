"""Signal analysis of the LP4 encoder: QMF band split, block envelope and the band MDCT.

All three stages are SSE single-precision code in Sony's encoder and are modelled here with
numpy float32 operations in the same order, which reproduces the binary bit for bit.

Data layout used throughout the model: a QMF frame is a float32 array X[256, 4]
(time index, band), i.e. the four band signals interleaved sample by sample, as in Sony's
working buffers.  Band 0 is 0-5.5 kHz, band 3 is 16.5-22 kHz.
"""
import numpy as np
from lp4_tables import QMF1, QMF2, MDCT_PERM, MDCT_TW, MDCT_FFT_TW, MDCT_SQRT_HALF, GAIN_POW2

F32 = np.float32
_Q1 = [np.array(v, F32) for v in QMF1]        # 0x48bfd0 + 16k, k = 0..11
_Q2 = [np.array(v, F32) for v in QMF2]        # 0x48c090 + 16k, k = 0..23
_NEG_HI = np.array([1, 1, -1, -1], F32)
_NEG_EVEN = np.array([-1, 1, -1, 1], F32)

QMF_STATE_LEN = 138       # 23 stage-1 output vectors (92 floats) + 46 past input samples


def qmf_analysis(x, state):
    """Split 1024 PCM samples (float32, int16 scale) into 4 bands of 256 samples.

    state: float32[138] carried between frames (initially zero).
    Returns (X[256, 4], new_state).

    Structure (0x435c50): a two-stage polyphase tree evaluated four input
    samples at a time.  Stage 1 forms, per group of 4 input samples, two 48-tap dot products
    SA and SB (12 vector multiplies each, summed in a fixed order) and combines their lanes
    into a vector u = [lo(t), lo(t+1), hi(t), hi(t+1)] style quantities; stage 2 runs a 24-tap
    vector filter over the last 24 stage-1 vectors and a final butterfly gives the 4 bands.
    """
    xh = np.concatenate([state[92:138], x]).astype(F32)     # xh[k] = input history + new samples
    s1_hist = state[0:92].reshape(23, 4)
    base = 40 + 4 * np.arange(256)

    def vec(o):                                              # 4 consecutive samples at offset o
        return np.stack([xh[base + o + j] for j in range(4)], 1)

    # stage 1: SA uses offsets +4, 0, -4, ..., -40 ; SB the same shifted by +2
    pa = [_Q1[10 - k] * vec(-4 * k) for k in range(11)] + [_Q1[11] * vec(4)]
    sa = pa[11] + pa[10]
    for k in range(9, -1, -1):
        sa = sa + pa[k]
    pb = [_Q1[10 - k] * vec(2 - 4 * k) for k in range(11)] + [_Q1[11] * vec(6)]
    sb = pb[11] + pb[10]
    for k in range(9, -1, -1):
        sb = sb + pb[k]
    u = np.stack([sa[:, 2] + sa[:, 0], sb[:, 2] + sb[:, 0], sa[:, 3] + sa[:, 1], sb[:, 3] + sb[:, 1]], 1)
    s1 = u[:, [2, 3, 0, 1]] * _NEG_HI + u                   # [u0+u2, u1+u3, u2-u0, u3-u1]
    # stage 2: 24-tap filter over stage-1 vectors (current one first, then 1..23 frames back)
    s1all = np.concatenate([s1_hist, s1])
    acc = _Q2[23] * s1
    taps = [_Q2[22 - k] * s1all[22 - k: 22 - k + 256] for k in range(22)] + [_Q2[0] * s1all[0:256]]
    for k in range(22, -1, -1):
        acc = acc + taps[k]
    r = acc[:, [1, 0, 3, 2]] + acc * _NEG_EVEN
    X = r[:, [1, 0, 2, 3]]
    new_state = np.concatenate([s1all[-23:].reshape(-1), xh[-46:]]).astype(F32)
    return np.ascontiguousarray(X), new_state


def envelope(X):
    """Peak magnitude of each block of 8 samples per band (0x435ab0): returns float32[4, 32]."""
    return np.abs(X).reshape(32, 8, 4).max(1).T.copy()


def mdct(time, state, first_levels):
    """Forward MDCT of all four bands (0x435380).

    time  : float32[512, 4] = [newest 256 samples | previous 256 samples] (Sony's buffer order);
            only the previous 256 samples (rows 256..511) enter this frame's transform.
    state : float32[128, 4] folded first-half data computed in the previous call.
    first_levels: gain level of the first gain point per band (4 = no gain point); the folded
            previous half is multiplied by 2**(level-4) so that it matches the gain-modified
            continuation of the signal.
    Returns (spec float32[4, 256] with odd bands reversed in frequency, new_state).

    The transform is computed as a pre-twiddle that folds 256 new samples into 128 complex
    values (half of the folding is kept as state for the next frame), a radix-2 complex FFT of
    length 128 computed on interleaved real data, and a post-processing step.  The result is
    numerically a 512-point MDCT (256 coefficients) per band with the windowing and scaling
    built into the twiddle tables.
    """
    g = np.array([GAIN_POW2[l] for l in first_levels], F32)
    k = np.array(MDCT_PERM)
    tw = np.array(MDCT_TW, F32)
    s2g = state * g
    xa = time[0x180 + k]
    xb = time[0x17f - k]
    w0, w1, w2, w3 = (tw[k, i][:, None] for i in range(4))
    new_state = w1 * (xa * w0 - xb)
    t = w2 * (xb * w0 + xa)
    d = w3 * (t - s2g)
    W = np.empty((256, 4), F32)
    W[0::2] = (t + s2g) + d
    W[1::2] = d
    t1 = np.array(MDCT_FFT_TW, F32)
    # first radix stage (groups of 4 vectors)
    V = W.reshape(64, 4, 4)
    T = t1[0:128:2][:, None]
    a, b, c, e = V[:, 0].copy(), V[:, 1].copy(), V[:, 2].copy(), V[:, 3].copy()
    x3 = T * (b - e)
    x0 = T * (a - c)
    V[:, 0] = (a + c) + x3
    V[:, 1] = (b + e) + x0
    V[:, 2] = x0
    V[:, 3] = x3
    # further stages, block size 2n for n = 4 .. 64
    n = 4
    while n != 128:
        nb = 256 // (2 * n)
        V = W.reshape(nb, 2 * n, 4)
        twd = t1[n // 2 - 1 + n * np.arange(nb)][:, None, None]
        lo = V[:, :n].copy()
        hi = V[:, n:].copy()
        hi2 = twd * (lo - hi)
        lo2 = lo + hi
        V[:, :n] = hi2[:, ::-1] + lo2
        V[:, n:] = hi2
        n *= 2
    lo = W[:128].copy()
    hi = W[128:].copy()
    W[128:] = F32(MDCT_SQRT_HALF) * (lo - hi)
    W[:128] = lo + hi
    i = np.arange(128)
    Wp = W[i] + W[255 - i]
    spec = np.empty((4, 256), F32)
    for band in range(4):
        if band % 2 == 0:
            spec[band, i] = Wp[:, band]
            spec[band, 255 - i] = W[255 - i, band]
        else:                                   # odd bands come out frequency-reversed
            spec[band, i] = W[255 - i, band]
            spec[band, 255 - i] = Wp[:, band]
    return spec, new_state
