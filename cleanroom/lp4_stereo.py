"""Chapter 8 §§6 and 9: joint-stereo analysis, conversion, and unit bit split."""
import math

import numpy as np

from .tables import load

T = load('lp4')


def f32(value):
    return float(np.float32(value))


class Stereo:
    def __init__(self, limit=1):
        self.limit = limit
        self.current = [-3] * 4
        self.sent = [0] * 4
        self.weight = 15
        self.sent_weight = 0
        self.side = [1.0] * 4
        self.energies = np.zeros((4, 2), np.float32)

    def process(self, left, right):
        # Inputs are already binary32 QMF samples; conversion mutates copies only.
        left, right = left.copy(), right.copy()
        accum = [f32(0.1), f32(0.1)]
        for band in range(4):
            sums = [sum(abs(float(v)) for v in channel[:, band]) for channel in (left, right)]
            e0, e1 = map(f32, sums)
            self.energies[band] = [e0, e1]
            k0 = k1 = float(T['JS_K_SIDE'])
            if band < self.limit:
                k0 = k1 = float(T['JS_K_MATRIX'])
                if self.current[band] == 0:
                    k0 = float(T['JS_K_HYST'])
                elif self.current[band] == 1:
                    k1 = float(T['JS_K_HYST'])
            if k0 * sums[0] < sums[1]:
                decision = 0
            elif k1 * sums[1] < sums[0]:
                decision = 1
            else:
                decision = -3 if sums[0] + sums[1] < float(T['JS_QUIET']) else 3
            self.sent[band], self.current[band] = abs(self.current[band]), decision
            if band >= self.limit and abs(decision) == 3:
                accum[0], accum[1] = 16 * accum[0] + e0, 16 * accum[1] + e1
        self.sent_weight = self.weight
        total = f32(accum[1] + accum[0])
        minimum, side = (accum[0], 0) if accum[0] <= accum[1] else (f32(accum[1]), 8)
        rounded = np.float32(minimum / total * 14 + 12582912.0)
        self.weight = (int(rounded.view(np.uint32)) & 7) + side
        if self.weight in (0, 8):
            self.weight += 1
        for band in range(4):
            d, p = abs(self.current[band]), self.sent[band]
            x0, x1 = left[:, band].astype(np.float64), right[:, band].astype(np.float64)
            if band < self.limit:
                e0, e1 = (float(v) for v in self.energies[band])
                r = 1.0
                factor = float(T['JS_RATIO_GAIN'])
                if d == 3:
                    if e1 > 3.94 * e0 and e0 * 256 > e1:
                        r = f32((e1 - e0) * factor / (e1 + e0))
                    elif e0 > 3.94 * e1 and e1 * 256 > e0:
                        r = f32((e0 - e1) * factor / (e1 + e0))
                    elif e1 == 0 or e0 == 0:
                        r = 1.0
                    elif e1 * 256 <= e0 or e0 * 256 <= e1:
                        r = f32(5 / 3)
                q = self.side[band]
                if r - q >= 0.05:
                    r = f32(q + 0.05)
                if q - r >= 0.05:
                    r = f32(q - 0.05)
                self.side[band] = r
                c0, c1 = (float(v) for v in T['CHSW_TABLE'][d:d + 2])
                p0, p1 = (float(v) for v in T['CHSW_TABLE'][p:p + 2])
                c1r = c1 * r
                dc, dw = f32((p0 - c0) / 8), (p1 * q - c1r) / 8
                mid = (x0 + x1) * 0.5
                difference = (x0 - mid * c0) * (1 / c1r)
                for t in range(8):
                    k = t - 8
                    difference[t] = (x0[t] - (c0 - dc * k) * mid[t]) / (c1r - k * dw)
                left[:, band], right[:, band] = mid, difference
            else:
                def weight_index(w):
                    return (w ^ 8) + 8 if d == 0 else w + 8 if d == 1 else w & 7

                tc = float(T['CHWT_TABLE'][weight_index(self.weight)])
                tp = float(T['CHWT_TABLE'][weight_index(self.sent_weight)])
                delta = (tp - tc) / 8
                intensity = (x1 + x0) * (1 / tc)
                for t in range(8):
                    intensity[t] = (x1[t] + x0[t]) / (tc - (t - 8) * delta)
                left[:, band] = intensity
        return left, right


def apply_gain(samples, points):
    """§7.4 gain application with specified points; detection is a separate stage."""
    output = samples.copy()
    start = 0
    for i, (level, location) in enumerate(points):
        next_level = points[i + 1][0] if i + 1 < len(points) else 4
        end = 8 * location
        output[start:end + 1] *= np.float32(2 ** (level - 4))
        for j in range(1, 8):
            factor = np.float32(2 ** ((8 * level + j * (next_level - level)) / 8 - 4))
            output[end + j] *= factor
        start = end + 8
    return output


def split_budget(spectrum0, spectrum1, matrix0, side_gain_count):
    a, b = spectrum0.reshape(1024), spectrum1.reshape(1024)
    edges, widths = T['SB_START'], T['SB_WIDTH']
    s0 = 0.0
    for i in range(18, -1, -1):
        s0 += int(widths[i]) * float(np.max(np.abs(a[int(edges[i]):int(edges[i + 1])])))
    s0 = f32(s0)
    s1 = 0.0
    for i in range(12):
        s1 += int(widths[i]) * float(np.max(np.abs(b[int(edges[i]):int(edges[i + 1])])))
    total_side = f32(2 * s1)
    if matrix0 != 3:
        total_side = f32(2 * total_side)
    total, base = 1477, 1133
    if total_side < 1 or (matrix0 != 3 and total_side < 4):
        return total
    if side_gain_count == 0 and ((3 * total_side < s0 and matrix0 == 3) or 6 * total_side < s0):
        spare = total - base
        if not 12 * total_side > s0:
            extra = spare
        else:
            if not 10 * total_side >= s0:
                spare = math.trunc(1.5 * spare)
            elif not 5 * total_side >= s0:
                spare = math.trunc(1.25 * spare)
            extra = math.trunc((s0 - total_side) / (2 * s0) * spare)
        return base + (extra & ~7)
    if total_side >= s0 or matrix0 != 3:
        return base
    return base + (math.trunc((s0 - total_side) / (10 * s0) * (total - base)) & ~7)
