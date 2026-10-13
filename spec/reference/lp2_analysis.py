"""LP2 encoder analysis front end: QMF band split, gain control, MDCT.

Arithmetic note: Sony's code runs on the x87 FPU with 53-bit precision.  Intermediate values
that stay in FPU registers are therefore IEEE doubles; values stored to memory are float32.
This module reproduces that exactly: float64 for register values, float32 for stored ones,
and the same operation order wherever the result is rounded.
"""
import collections, math
import numpy as np
from lp2_tables import (QMF_WIN, MDCT_WIN, PRE_COS, PRE_SIN, FFT_COS, FFT_SIN, BITREV7,
                        POST1, POST2, POST3, POST4, GAIN_LEVEL_EXP)

f32, f64 = np.float32, np.float64
STATS = collections.Counter()      # branch statistics (informational only)
_QW = QMF_WIN.astype(f64)


# ================================================================ QMF analysis (0x440bc0)
def qmf_split(x, delay):
    """Split x (float32, even length n) into a low and a high half-rate signal.
    `delay` (46 float32) holds the last 46 input samples of the previous call; updated in place.
    out_lo[j] = Q + P, out_hi[j] = P - Q with
      P = sum over odd k (47, 45, ..., 1) of w[k] * s[2j + k]   (accumulated in that order)
      Q = sum over even k (46, 44, ..., 0) of w[k] * s[2j + k]
    where s = delay ++ x."""
    n = len(x)
    s = np.concatenate([delay, x]).astype(f64)
    delay[:] = x[n - 46:]
    m = n // 2
    P = np.zeros(m); Q = np.zeros(m)
    for k in range(47, 0, -2):
        P = P + _QW[k] * s[k:k + 2 * m:2]
    for k in range(46, -1, -2):
        Q = _QW[k] * s[k:k + 2 * m:2] + Q
    return (Q + P).astype(f32), (P - Q).astype(f32)


# ================================================================ gain control
class GainInfo:
    """Gain-control data of one QMF band for one frame (the bitstream's gain points)."""
    __slots__ = ('pts', 'peak')

    def __init__(self):
        self.pts = []          # list of (level_code 0..15, location 0..31)
        self.peak = f32(0.0)   # peak |x| of the band's current 256-sample block (float32)


def _log2_round(r):
    """Level from a ratio: int(log2(max(r, 1e-20)) + 0.5) as computed by the reference
    (log10 via fyl2x divided by log10(2), +0.5, truncation)."""
    if r <= 1e-20:
        r = 1e-20
    lg2 = math.log10(2.0)
    a = (math.log(r) / math.log(2.0)) * lg2
    return int(a / lg2 + 0.5)


def band_gain_points(x, prev):
    """Gain-control analysis of one band (0x440d20).
    x: 768-sample band history [previous block, current block, next block] (float32).
    prev: GainInfo of the previous frame.  Returns the GainInfo of the current frame, which
    describes the transitions inside the current block (and its boundary to the next one)."""
    blk = np.abs(x[256:768]).reshape(64, 8).max(axis=1).astype(f32)  # 8-sample peaks
    pre = f32(np.abs(x[224:256]).max())                              # tail of previous block
    prev_lev = 0
    for lev, loc in prev.pts:
        prev_lev = max(prev_lev, GAIN_LEVEL_EXP[lev])
    grp = blk.reshape(16, 4).max(axis=1)                             # 32-sample peaks
    k185 = f64(f32(1.85))

    # --- candidate list: attack ratios R[0..31] and release ratios R[32..39]
    R = np.zeros(40, f32)
    cnt = 0
    M = f64(prev.peak)
    for i in range(32):
        if M < blk[i]: M = f64(blk[i])
        v = f64(blk[i + 1])
        if v > 10.0 and v > 1.5 * M:
            R[i] = f32(v * 0.25) if M <= 4.0 else f32(v / M)
            cnt += 1
    M2 = f64(blk[63])
    for g in range(15, 7, -1):
        if M2 < grp[g]: M2 = f64(grp[g])
    for g in range(7, -1, -1):
        if M2 < grp[g]: M2 = f64(grp[g])
        v = f64(pre) if g == 0 else f64(grp[g - 1])
        if v > 10.0 and v > M2 * k185:
            R[32 + g] = f32(0.25 * v) if M2 <= 4.0 else f32(v / M2)
            cnt += 1
    sel = np.zeros(40, np.int32)            # keep the 7 strongest candidates
    for _ in range(min(cnt, 7)):
        best = f64(0.0); idx = 0
        for i in range(40):
            if best < R[i]: best = f64(R[i]); idx = i
        R[idx] = 0; sel[idx] = 1

    # --- attacks: level steps (log2 of the jump, rounded), at most 7, total level <= 10
    att = []
    tot = prev_lev
    M = f64(prev.peak)
    for i in range(32):
        if M < blk[i]: M = f64(blk[i])
        v = f64(blk[i + 1])
        if not (v > 10.0 and v > 1.5 * M and sel[i]):
            continue
        a = _log2_round(v * 0.25 if M <= 4.0 else v / M)
        if a >= 0:
            if a + tot > 10: a = 10 - tot
            tot += a
            if a > 0:
                att.append([i, a])
        if len(att) >= 7:
            break

    # --- releases (decays), scanned backwards in 32-sample groups, total <= min(tot, 4)
    cap = min(tot, 4)
    rel = []
    if att or prev.pts:
        M2 = f64(blk[63]); rs = 0; n_all = len(att)
        for g in range(15, 7, -1):
            if M2 < grp[g]: M2 = f64(grp[g])
        for g in range(7, -1, -1):
            if M2 < grp[g]: M2 = f64(grp[g])
            v = f64(pre) if g == 0 else f64(grp[g - 1])
            if v > 10.0 and v > M2 * k185 and sel[32 + g]:
                a = _log2_round(v * 0.25 if M2 <= 4.0 else v / M2)
                if a >= 0:
                    if a + rs > cap: a = cap - rs
                    rs += a
                    if a > 0 and len(att) < 7:
                        rel.append([1 if g == 0 else 4 * g, a])
                        n_all += 1
            if n_all >= 7:
                break

    out = GainInfo()
    out.peak = f32(max(0.0, float(blk[:32].max())))
    # cumulative levels (from the end), releases capped at 4, attacks at 11 + total release
    acc = 0
    for k in range(len(rel) - 1, -1, -1):
        rel[k][1] = min(rel[k][1] + acc, 4); acc = rel[k][1]
    cap_a = 11 + acc
    acc = 0
    for k in range(len(att) - 1, -1, -1):
        att[k][1] = min(att[k][1] + acc, cap_a); acc = att[k][1]
    lv = [0] * 33                      # log2 gain per 8-sample position
    j = 0
    for pos, val in att:
        while j <= pos:
            lv[j] += val; j += 1
    j = 32
    for pos, val in rel:
        while j >= pos:
            lv[j] += val; j -= 1
    pts = [(lv[j] - lv[32] + 4, j) for j in range(32) if lv[j] != lv[j + 1]]
    while pts and pts[-1][0] == 4:
        pts.pop()
    STATS['gain_attacks'] += len(att); STATS['gain_releases'] += len(rel)
    STATS['gain_points'] += len(pts)
    if len(pts) > 7:
        raise ValueError('more than 7 gain points')
    out.pts = pts
    return out


def adjust_band0_gain(hist, prev, cur):
    """adjust_band0_gain (0x4414f0): give band 0 an attack point when band 1 has a strong attack
    (level spread >= 2) and band 0 has none and is not loud."""
    spread = []
    for b in range(4):
        p = cur[b].pts
        if not p:
            spread.append(0); continue
        mn, mi = 4, 0
        for i, (lev, loc) in enumerate(p):
            if lev < mn: mn, mi = lev, i
        mx = max([0] + [p[i][0] for i in range(mi + 1)])
        spread.append(mx - mn)
    if spread[1] < 2 or prev[0].pts or cur[0].pts:
        return
    loc = cur[1].pts[0][1]
    ok = not (prev[0].peak > 16384.0)
    if np.any(np.abs(hist[0][256:256 + 8 * (loc + 1)].astype(f64)) > 16384.0):
        ok = False
    if ok:
        cur[0].pts = [(5, loc)]
        STATS['adjust_band0_copies'] += 1


def gain_window(prev, cur):
    """gain_window (0x440470): 512-sample divisor applied before the MDCT.
    First half from the previous frame's points, second half from the current frame's."""
    G = [0] * 64
    p = 0
    for lev, loc in cur.pts:
        e = GAIN_LEVEL_EXP[lev]
        while p <= loc + 32:
            G[p] = e; p += 1
    p = 0
    for lev, loc in prev.pts:
        e = GAIN_LEVEL_EXP[lev]
        while p <= loc:
            G[p] += e; p += 1
    out = np.empty(512, f32)
    o = 511; last = 0
    for k in range(63, -1, -1):            # built backwards in 8-sample blocks
        e = -G[k]
        if e == last:
            v = f32(1.0 / (1 << -e)) if e < 0 else f32(1 << e)
            out[o - 7:o + 1] = v; o -= 8
        else:                              # 8-step geometric interpolation
            for m in range(8):
                out[o] = f32(math.pow(2.0, ((m + 1) * e + (7 - m) * last) * 0.125)); o -= 1
        last = e
    return out


# ================================================================ MDCT (mdct256, 0x440620)
def mdct256(x, odd):
    """512 windowed samples -> 256 coefficients, via a 128-point complex FFT.
    Odd QMF bands are returned in reversed order (spectral inversion)."""
    w = MDCT_WIN.astype(f64) * np.asarray(x, f32).astype(f64)
    z = np.empty(256, f32)
    n = np.arange(64)
    z[0:64] = (-w[384 + 2 * n]) - w[383 - 2 * n]
    m = np.arange(128)
    z[64:192] = w[2 * m] - w[255 - 2 * m]
    z[192:256] = w[256 + 2 * n] + w[511 - 2 * n]
    zd = z.astype(f64)
    c, s = PRE_COS.astype(f64), PRE_SIN.astype(f64)
    re = (c * zd[0::2] - s * zd[1::2]).astype(f32)
    im = (s * zd[0::2] + c * zd[1::2]).astype(f32)
    re = re[BITREV7]; im = im[BITREV7]
    groups, span = 64, 2
    for _ in range(7):                     # radix-2 decimation-in-time butterflies
        half = span // 2
        tw = np.arange(half) * groups
        tc, ts = FFT_COS[tw].astype(f64), FFT_SIN[tw].astype(f64)
        for g in range(groups):
            top = g * span + np.arange(half); bot = top + half
            rb, ib = re[bot].astype(f64), im[bot].astype(f64)
            tr = tc * rb - ts * ib            # kept in double (FPU register)
            ti = ts * rb + tc * ib
            rt, it = re[top].astype(f64), im[top].astype(f64)
            re[bot] = (rt - tr).astype(f32); im[bot] = (it - ti).astype(f32)
            re[top] = (tr + rt).astype(f32); im[top] = (ti + it).astype(f32)
        groups //= 2; span *= 2
    R, I = re.astype(f64), im.astype(f64)
    k = np.arange(128)
    p1, p2, p3, p4 = (t.astype(f64) for t in (POST1, POST2, POST3, POST4))
    y = np.empty(256, f32)
    y[k] = ((R[127 - k] * p2 + I[k] * p3) + p1 * R[k]) + p4 * I[127 - k]
    y[255 - k] = ((p3 * R[k] - p4 * R[127 - k]) - p1 * I[k]) + p2 * I[127 - k]
    out = (y.astype(f64) * 0.0078125).astype(f32)
    return out[::-1].copy() if odd else out


# ================================================================ per-channel front end
class ChannelAnalysis:
    """Analysis state for one channel."""

    def __init__(self):
        self.qmf_delay = [np.zeros(46, f32) for _ in range(3)]
        self.hist = [np.zeros(768, f32) for _ in range(4)]   # per band: 3 blocks of 256
        self.overlap_gc = [np.zeros(256, f32) for _ in range(4)]
        self.overlap_plain = [np.zeros(256, f32) for _ in range(4)]
        self.prev_gain = [GainInfo() for _ in range(4)]

    def process(self, pcm):
        """pcm: 1024 float32 samples.  Returns (spec_gc, spec_plain, gain, prev_gain):
        spec_gc    MDCT of the gain-compensated block (what gets quantised),
        spec_plain MDCT without gain compensation (used for tonal / masking decisions),
        gain       GainInfo per band for this frame, prev_gain the previous frame's."""
        lo, hi = qmf_split(np.asarray(pcm, f32), self.qmf_delay[0])
        b0, b1 = qmf_split(lo, self.qmf_delay[1])
        b3, b2 = qmf_split(hi, self.qmf_delay[2])
        for b, s in enumerate((b0, b1, b2, b3)):
            self.hist[b] = np.concatenate([self.hist[b][256:], s])
        prev = self.prev_gain
        cur = [band_gain_points(self.hist[b], prev[b]) for b in range(4)]          # gain control
        adjust_band0_gain(self.hist, prev, cur)
        spec_gc = np.empty(1024, f32)
        for b in range(4):                                                  # MDCT
            mid = self.hist[b][256:512]
            if not cur[b].pts and not prev[b].pts:
                buf = np.concatenate([self.overlap_gc[b], mid])
            else:
                w = gain_window(prev[b], cur[b])
                buf = np.concatenate([(self.overlap_gc[b].astype(f64) / w[:256]).astype(f32),
                                      (mid.astype(f64) / w[256:]).astype(f32)])
            spec_gc[b * 256:(b + 1) * 256] = mdct256(buf, b & 1)
            self.overlap_gc[b] = mid.copy()
        if any(cur[b].pts or prev[b].pts for b in range(4)):
            spec_plain = np.empty(1024, f32)
            for b in range(4):
                mid = self.hist[b][256:512]
                spec_plain[b * 256:(b + 1) * 256] = mdct256(
                    np.concatenate([self.overlap_plain[b], mid]), b & 1)
                self.overlap_plain[b] = mid.copy()
        else:
            spec_plain = spec_gc.copy()
            for b in range(4):
                self.overlap_plain[b] = self.hist[b][256:512].copy()
        self.prev_gain = cur
        return spec_gc, spec_plain, cur, prev
