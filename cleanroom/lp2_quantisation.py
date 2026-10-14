"""Chapter 7 §§9–11: scale factors, quantisers, bit costs, and priority ordering."""
import math

import numpy as np

from .bitstream import codewords, Tone
from .tables import load

T = load('lp2')


def scale_index(magnitude):
    return min(int(np.searchsorted(T['SF'], magnitude, side='right')), 63)


def group_scales(spectrum):
    peaks = np.max(np.abs(spectrum).reshape(256, 4), axis=1)
    return [scale_index(p) for p in peaks]


def subband_scales(groups, nsub):
    return [max(groups[int(T['SB_START'][i]) // 4:int(T['SB_START'][i + 1]) // 4])
            for i in range(nsub)]


def quantise(values, selector, dead_index=0):
    if not selector:
        return np.zeros(len(values), np.int32), 0
    maximum = int(T['MAXQ'][selector])
    q = np.clip(np.trunc(np.asarray(values, np.float64) * (maximum + 0.5) + 31.5) - 31,
                -maximum, maximum).astype(np.int32)
    zone = np.float32(float(T['T1'][min(dead_index, 12)]) / float(T['T2'][selector])) if dead_index else 0
    q[np.abs(values) < zone] = 0
    return q, 6 + sum(width for _, width in codewords(q, selector, 0, T))


def extract_component(spectrum, position):
    values = spectrum[position:position + 4].astype(np.float64)
    scale = scale_index(np.float32(np.max(np.abs(values))))
    q = np.clip(np.trunc(31.5 * (values / float(T['SF'][scale])) + 31.5) - 31,
                -31, 31).astype(np.int32)
    step = float(np.float32(float(T['SF'][scale]) / 31.5))
    spectrum[position:position + 4] = (values - q.astype(np.float64) * step).astype(np.float32)
    return Tone(position, scale, q.tolist())


def descending_order(keys):
    # Chapter 7 §10.9 explicitly specifies swaps of indices, including ties.
    values, indices = list(keys), list(range(len(keys)))
    gap = 1
    while gap <= len(values):
        gap = 3 * gap + 1
    while True:
        gap //= 3
        if gap == 0:
            return indices
        if gap >= len(values):
            continue
        for t in range(len(values) - gap):
            value = values[gap + t]
            j = t
            while j >= 0 and values[j] < value:
                values[j + gap] = values[j]
                indices[j + gap], indices[j] = indices[j], indices[j + gap]
                j -= gap
            values[j + gap] = value


def word_lengths(offsets, scales, slope):
    peak = max(scales[1:])
    quiet = 6 if peak >= 30 else max(1, peak // 6)
    selectors = []
    for i, offset in enumerate(offsets):
        value = min(7, max(1, math.trunc(float(offset) + 0.5)))
        if scales[i] < quiet or any(scales[i] < scales[j] - slope for j in range(max(0, i - 8), i)):
            value = 0
        selectors.append(value)
    return selectors, quiet
