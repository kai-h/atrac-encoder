"""Chapter 8 §7: stateful LP4 gain detection and sample compensation."""
from dataclasses import dataclass, field

import numpy as np

from .lp4_stereo import T, f32, apply_gain


def float_bits(value):
    with np.errstate(over='ignore'):
        return int(np.float32(value).view(np.uint32))


def gain_mode(band, joint, sent, current):
    if not joint:
        return band
    return 5 if sent != current else 0 if sent != 3 else -1


@dataclass
class Band:
    points: list = field(default_factory=list)
    envelope: np.ndarray = field(default_factory=lambda: np.zeros(32, np.float32))
    peak: float = 0.0
    last: float = 0.0

    def detect(self, newest, mode):
        old_envelope = self.envelope.copy()
        a = np.concatenate([old_envelope, newest]).astype(np.float64)
        blocks = [float(np.max(a[4 * j:4 * j + 4])) for j in range(8)]
        before = self.last
        self.last = blocks[7]
        end = 64 if mode <= 0 else 8 * (8 - mode)
        floor, minimum, ratio = (float(T[name]) for name in ('GC_FLOOR', 'GC_MIN', 'GC_RATIO'))
        ref = max([floor, blocks[7]] + a[32:end].tolist())
        threshold = f32(ratio * ref)
        budget = 4
        decays = []

        def level(z):
            return (float_bits(z * float(T['GC_SQRT2'])) >> 23) - 127

        for j in range(7, -2, -1):
            value = blocks[j] if j >= 0 else before
            if value < ref:
                continue
            if value > minimum and value > threshold:
                e = 4 * (j + 1)
                if j == -1 or a[e] >= threshold or a[e - 1] >= threshold:
                    location = e
                else:
                    location = e - 1 if a[e - 2] >= threshold else e - 2
                step = min(level(value / ref), budget)
                budget -= step
                decays.insert(0, (location, -step))
                if budget <= 0 or len(decays) == 2:
                    break
            ref = value
            threshold = f32(value * ratio)
        stop = decays[0][0] if decays else 32
        room = min(15, 131 - (float_bits(self.peak * float(T['GC_EXPSCALE'])) >> 23)) - budget
        attacks = []
        if room > 0:
            factor = 2.0 if mode == -1 else f32(1.6)
            ref = max(self.peak, float(a[0]))
            threshold = factor * ref
            for i in range(stop):
                value = float(a[i + 1])
                if value < ref:
                    continue
                if value > minimum and value > threshold:
                    step = level(value / ref if ref else float('inf'))
                    if attacks and attacks[-1][0] == i - 1 and step >= attacks[-1][1]:
                        _, merged = attacks.pop()
                        room += merged
                        step += merged
                    step = min(step, room)
                    room -= step
                    attacks.append((i, step))
                    if len(attacks) == 7 - len(decays) or room <= 0:
                        break
                ref = value
                threshold = value * factor
        events = attacks + decays
        points = []
        level_sum = 4
        for location, step in reversed(events):
            level_sum += step
            points.append((level_sum, location))
        self.points = list(reversed(points))
        self.peak = max(floor, float(np.max(old_envelope)))
        self.envelope = np.asarray(newest, np.float32).copy()
        return old_envelope


class GainControl:
    def __init__(self):
        self.bands = [Band() for _ in range(4)]

    def process(self, newest, previous, joint=False, sent=None, current=None):
        peaks = np.max(np.abs(newest).reshape(32, 8, 4), axis=1).T
        old_count0 = len(self.bands[0].points)
        old_band0 = None
        for band in (3, 2, 1, 0):
            mode = gain_mode(band, joint, sent[band] if joint else 0, current[band] if joint else 0)
            old = self.bands[band].detect(peaks[band], mode)
            if band == 0:
                old_band0 = old
        zero, one = self.bands[:2]
        bound = float(T['GC_BAND0_MAX'])
        if not joint and not zero.points and not old_count0 and one.points and zero.peak <= bound:
            levels = [4] + [v for v, _ in one.points]
            location = one.points[0][1]
            if max(levels) - min(levels) >= 2 and np.max(old_band0[:location + 1]) <= bound:
                zero.points = [(5, location)]
        modified = np.stack([apply_gain(previous[:, b], self.bands[b].points) for b in range(4)], axis=1)
        return modified, old_count0
