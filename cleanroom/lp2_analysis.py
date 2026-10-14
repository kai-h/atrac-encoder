"""Chapter 7 §§4–8: LP2 QMF, gain analysis, and FFT-based MDCT."""
from dataclasses import dataclass, field
import math

import numpy as np

from .tables import load

T = load('lp2')


@dataclass
class Gain:
    points: list = field(default_factory=list)  # (level, location)
    peak: float = 0.0


def profile(attacks, releases):
    result = [0] * 33
    position = 0
    for end, value in attacks:
        while position <= end:
            result[position] += value
            position += 1
    position = 32
    for start, value in releases:
        while position >= start:
            result[position] += value
            position -= 1
    return result


def profile_points(values):
    points = [(values[i] - values[32] + 4, i)
              for i in range(32) if values[i] != values[i + 1]]
    while points and points[-1][0] == 4:
        points.pop()
    if len(points) > 7:
        raise ValueError('Gain profile exceeds seven points')
    return points


def detect_gain(history, previous):
    block = np.max(np.abs(history[256:]).reshape(64, 8), axis=1)
    group = np.max(block.reshape(16, 4), axis=1)
    pre = float(np.max(np.abs(history[224:256])))
    prior_level = max([0] + [level - 4 for level, _ in previous.points])
    ratios = np.zeros(40, np.float32)
    exact = [0.0] * 40
    maximum = float(previous.peak)
    for i in range(32):
        maximum = max(maximum, float(block[i]))
        value = float(block[i + 1])
        if value > 10 and value > 1.5 * maximum:
            exact[i] = value / (maximum if maximum > 4 else 4)
            ratios[i] = exact[i]
    maximum = max(float(block[63]), *(float(x) for x in group[8:]))
    for g in range(7, -1, -1):
        maximum = max(maximum, float(group[g]))
        value = float(group[g - 1]) if g else pre
        if value > 10 and value > float(np.float32(1.85)) * maximum:
            exact[32 + g] = value / (maximum if maximum > 4 else 4)
            ratios[32 + g] = exact[32 + g]
    chosen = set()
    for _ in range(7):
        i = int(np.argmax(ratios))
        if ratios[i] == 0:
            break
        chosen.add(i)
        ratios[i] = 0

    def level(i):
        return math.trunc(math.log10(max(exact[i], 1e-20)) / math.log10(2) + 0.5)

    attacks = []
    total = prior_level
    for i in range(32):
        if i not in chosen:
            continue
        step = level(i)
        if step >= 0:
            step = min(step, 10 - total)
            total += step
            if step > 0:
                attacks.append([i, step])
            if len(attacks) == 7:
                break
    releases = []
    if attacks or previous.points:
        cap = min(total, 4)
        total_release = 0
        for g in range(7, -1, -1):
            if 32 + g not in chosen:
                continue
            step = level(32 + g)
            if step >= 0:
                step = min(step, cap - total_release)
                total_release += step
                if step > 0 and len(attacks) < 7:
                    releases.append([4 * g if g else 1, step])
            if len(attacks) + len(releases) == 7:
                break
    for events, cap in ((releases, 4), (attacks, 11 + min(4, sum(v for _, v in releases)))):
        running = 0
        for event in reversed(events):
            running = min(cap, running + event[1])
            event[1] = running
    return Gain(profile_points(profile(attacks, releases)), float(np.max(block[:32])))


def gain_window(previous, current):
    exponents = [0] * 64
    for points, shift, add in ((current.points, 32, False), (previous.points, 0, True)):
        start = 0
        for level, location in points:
            end = location + shift
            for k in range(start, end + 1):
                exponents[k] = exponents[k] + level - 4 if add else level - 4
            start = end + 1
    window = np.empty(512, np.float32)
    last = 0
    for k in range(63, -1, -1):
        exponent = -exponents[k]
        for m in range(8):
            window[8 * k + 7 - m] = 2 ** (((m + 1) * exponent + (7 - m) * last) / 8)
        last = exponent
    return window


def mdct(samples, odd=False):
    windowed = np.asarray(samples, np.float64) * T['MDCT_WIN'].astype(np.float64)
    z = np.empty(256, np.float32)
    n = np.arange(64)
    z[:64] = -windowed[384 + 2 * n] - windowed[383 - 2 * n]
    m = np.arange(128)
    z[64:192] = windowed[2 * m] - windowed[255 - 2 * m]
    z[192:] = windowed[256 + 2 * n] + windowed[511 - 2 * n]
    c, s = T['PRE_COS'].astype(np.float64), T['PRE_SIN'].astype(np.float64)
    even, odd_values = z[::2].astype(np.float64), z[1::2].astype(np.float64)
    re = (c * even - s * odd_values).astype(np.float32)[T['BITREV7']]
    im = (s * even + c * odd_values).astype(np.float32)[T['BITREV7']]
    span = 2
    while span <= 128:
        half = span // 2
        for start in range(0, 128, span):
            for b in range(half):
                top, bottom = start + b, start + b + half
                idx = b * (128 // span)
                c, s = float(T['FFT_COS'][idx]), float(T['FFT_SIN'][idx])
                tr = c * float(re[bottom]) - s * float(im[bottom])
                ti = s * float(re[bottom]) + c * float(im[bottom])
                r, i = float(re[top]), float(im[top])
                re[bottom], im[bottom] = r - tr, i - ti
                re[top], im[top] = tr + r, ti + i
        span *= 2
    re, im = re.astype(np.float64), im.astype(np.float64)
    p1, p2, p3, p4 = (T[f'POST{k}'].astype(np.float64) for k in range(1, 5))
    y = np.empty(256, np.float32)
    y[:128] = ((re[::-1] * p2 + im * p3) + p1 * re) + p4 * im[::-1]
    y[255:127:-1] = ((p3 * re - p4 * re[::-1]) - p1 * im) + p2 * im[::-1]
    y /= np.float32(128)
    return y[::-1].copy() if odd else y


class Analysis:
    def __init__(self):
        self.delays = [np.zeros(46, np.float32) for _ in range(3)]
        self.history = np.zeros((4, 768), np.float32)
        self.previous_gain = [Gain() for _ in range(4)]

    def split(self, samples, stage):
        samples = np.asarray(samples, np.float32)
        combined = np.concatenate([self.delays[stage], samples]).astype(np.float64)
        self.delays[stage] = samples[-46:].copy()
        count = len(samples) // 2
        p, q = np.zeros(count), np.zeros(count)
        for k in range(47, 0, -2):
            p = p + float(T['QMF_WIN'][k]) * combined[k:k + 2 * count:2]
        for k in range(46, -1, -2):
            q = float(T['QMF_WIN'][k]) * combined[k:k + 2 * count:2] + q
        return (q + p).astype(np.float32), (p - q).astype(np.float32)

    def process(self, pcm):
        if np.shape(pcm) != (1024,):
            raise ValueError('Analysis requires exactly 1024 samples')
        low, high = self.split(pcm, 0)
        b0, b1 = self.split(low, 1)
        b3, b2 = self.split(high, 2)
        bands = np.stack([b0, b1, b2, b3])
        self.history[:, :512] = self.history[:, 256:].copy()
        self.history[:, 512:] = bands
        previous = self.previous_gain
        current = [detect_gain(self.history[b], previous[b]) for b in range(4)]
        # §6.4: minimum starts at unity; track the maximum before its first occurrence.
        if current[1].points:
            levels = [v for v, _ in current[1].points]
            minimum = min([4] + levels)
            end = levels.index(minimum) + 1 if minimum in levels else len(levels)
            spread = max([0] + levels[:end]) - minimum
            location = current[1].points[0][1]
            if (spread >= 2 and not current[0].points and not previous[0].points
                    and previous[0].peak <= 16384
                    and np.max(np.abs(self.history[0, 256:256 + 8 * (location + 1)])) <= 16384):
                current[0].points = [(5, location)]
        a, b = np.empty((4, 256), np.float32), np.empty((4, 256), np.float32)
        for band in range(4):
            raw = self.history[band, :512]
            b[band] = mdct(raw, band % 2 == 1)
            if previous[band].points or current[band].points:
                compensated = (raw.astype(np.float64) / gain_window(previous[band], current[band])).astype(np.float32)
                a[band] = mdct(compensated, band % 2 == 1)
            else:
                a[band] = b[band]
        self.previous_gain = current
        return a.reshape(1024), b.reshape(1024), current, previous
