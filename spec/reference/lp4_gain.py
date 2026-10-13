"""Gain control (transient pre-processing) of the LP4 encoder (x87 code, 0x439820/0x439bb0).

Per channel and QMF band the encoder keeps the 32 block peaks (blocks of 8 samples) of the
previous frame's new samples.  It looks at those 32 peaks plus a look-ahead into the peaks of
the samples that just arrived, finds level drops (decays after a loud part, scanning
backwards) and level rises (attacks, scanning forwards), turns them into at most 7 gain points
(location in units of 8 samples, level 0..15 where 4 means unity), and multiplies the
previous frame's 256 band samples by the corresponding gain curve before they are transformed.
The decoder undoes the curve with the transmitted points.
"""
import numpy as np
from lp4_fp import f32, bits, from_bits
from lp4_tables import (GC_FLOOR, GC_MIN, GC_RATIO, GC_SQRT2, GC_EXPSCALE, GC_BAND0_MAX,
                        GAIN_POW2, GAIN_RAMP_BITS)

F32 = np.float32


class BandGain:
    """Per-band gain state (0xcc bytes at unit+0x2a28+band*0xcc in Sony's layout)."""
    def __init__(self):
        self.loc = [0] * 8          # +0x00 gain point locations (units of 8 samples)
        self.lev = [0] * 9          # +0x20 gain levels, lev[count] = 4 sentinel
        self.env = [0.0] * 32       # +0x40 block peaks of the newest 256 samples
        self.peak = 0.0             # +0xc0 max(floor, peaks of the previous-but-one frame)
        self.last_block = 0.0       # +0xc4 peak of the last 32-sample block of the analysed frame
        self.count = 0              # +0xc8 number of gain points


def _log2_level(ratio):
    """floor(log2(ratio * sqrt2)) from the float32 exponent field, i.e. round(log2(ratio))."""
    return (bits(f32(ratio * GC_SQRT2)) >> 23) - 127


def detect(G, A, mode):
    """Find gain points for one band.  A = 64 block peaks [previous frame 32 | new frame 32].

    mode: band index for ordinary units (bands 1..3 shorten the look-ahead), or for the second
    (joint-stereo) unit -1 / 0 / 5 depending on whether the stereo decision of the band changed.
    Fills G.loc/G.lev and returns the number of points.
    """
    bm = [max(A[4 * j:4 * j + 4]) for j in range(8)]        # 32-sample block peaks
    prev_last = G.last_block
    G.last_block = bm[7]
    end = 64 if mode <= 0 else (8 - mode) * 8
    ref = max([bm[7]] + list(A[32:end]))                     # look-ahead reference level
    if ref < GC_FLOOR:
        ref = GC_FLOOR
    thr = f32(GC_RATIO * ref)
    budget = 4                                               # total attenuation available (levels)
    slot = 7                                                 # decay points are stored from slot 6 down
    G.loc[7] = 32
    blk = [prev_last] + bm
    # part 1: scan backwards for decays (blocks louder than 1.85x what follows them)
    for e in range(8, -1, -1):
        v = blk[e]
        if v < ref:
            continue
        if v > GC_MIN and v > thr:
            loc = 4 * e
            if e != 0:
                if A[4 * e] >= thr or A[4 * e - 1] >= thr:
                    loc = 4 * e
                elif A[4 * e - 2] >= thr:
                    loc = 4 * e - 1
                else:
                    loc = 4 * e - 2
            slot -= 1
            G.loc[slot] = loc
            ex = min(_log2_level(v / ref), budget)
            budget -= ex
            G.lev[slot] = -ex                                # relative step, made absolute below
            if budget <= 0 or slot == 5:
                break
        ref = v
        thr = f32(v * GC_RATIO)
    # part 2: scan forwards for attacks before the first decay point
    n = 0
    room = min(15, 0x83 - (bits(f32(G.peak * GC_EXPSCALE)) >> 23)) - budget
    if room > 0:
        fac = 2.0 if mode == -1 else f32(1.6)
        ref = G.peak
        if ref < A[0]:
            ref = A[0]
        thr = fac * ref
        i = 0
        while i < G.loc[slot]:
            a = A[i + 1]
            if not (a < ref):
                if a > GC_MIN and a > thr:
                    G.loc[n] = i
                    ex = _log2_level(a / ref)
                    if n >= 1 and G.loc[n - 1] == i - 1 and ex >= G.lev[n - 1]:
                        n -= 1                               # merge with an adjacent attack
                        G.loc[n] = i
                        room += G.lev[n]
                        ex += G.lev[n]
                    ex = min(ex, room)
                    room -= ex
                    G.lev[n] = ex
                    n += 1
                    if n == slot or room <= 0:
                        break
                ref = a
                thr = a * fac
            i += 1
    while slot < 7:                                          # append the decay points
        G.lev[n] = G.lev[slot]
        G.loc[n] = G.loc[slot]
        slot += 1
        n += 1
    acc = 4                                                  # cumulative levels, last point -> 4
    for k in range(n - 1, -1, -1):
        acc += G.lev[k]
        G.lev[k] = acc
    return n


def gain_band(bands, b, env_new, Xprev, js_unit, dec_prev, dec_cur, unit):
    """Gain control of band b of one unit (0x439820).  Modifies Xprev[:, b] (the previous
    frame's samples, about to be transformed) and bands[b]; returns the point count."""
    G = bands[b]
    A = list(G.env) + [float(v) for v in env_new]
    peak = max([GC_FLOOR] + list(G.env))
    G.env = [float(v) for v in env_new]
    mode = b
    if js_unit:
        if dec_prev[b] != dec_cur[b]:
            mode = 5
        elif dec_prev[b] != 3:
            mode = 0
        else:
            mode = -1
    n = detect(G, A, mode)
    G.peak = peak
    if b == 0:
        unit.prev_count0 = G.count
    if n == 0:
        # band 0 may copy a strong attack found in band 1 (bass smearing protection)
        g1 = bands[1]
        if b == 0 and not js_unit and G.count == 0 and g1.count != 0 and peak <= GC_BAND0_MAX:
            levs = [4] + g1.lev[:g1.count]
            if max(levs) - min(levs) < 2:
                return 0
            l0 = g1.loc[0]
            for k in range(l0 + 1):
                if A[k] > GC_BAND0_MAX:
                    return 0
            n = 1
            G.loc[0] = l0
            G.lev[0] = 5
        else:
            return 0
    apply_gain(G, n, Xprev[:, b])
    return n


def apply_gain(G, n, col):
    """Multiply a band signal by the gain curve: constant 2**(lev-4) up to and including sample
    8*loc, then a 7-sample geometric ramp towards the next level (finally towards unity)."""
    lev = G.lev[:n] + [4]
    pos = 0
    for i in range(n + 1):
        c = lev[i]
        if c == 4:
            if i == n:
                return
            pos = G.loc[i] * 8 + 1
        else:
            end = G.loc[i] * 8
            p = F32(GAIN_POW2[c])
            while pos < end:
                col[pos:pos + 8] = col[pos:pos + 8] * p
                pos += 8
            col[pos] = col[pos] * p
            pos += 1
        d = lev[i + 1] - c
        v = c * 8 + d
        for j in range(7):
            g = from_bits((v << 20) + GAIN_RAMP_BITS[v & 7])     # = 2**(v/8 - 4)
            col[pos] = col[pos] * F32(g)
            pos += 1
            v += d
