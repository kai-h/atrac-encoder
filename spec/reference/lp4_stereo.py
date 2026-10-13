"""Joint-stereo stage of the LP4 encoder (x87 code, modelled with doubles + float32 stores).

Addresses: JointStereo.analyse (decision, 0x439630), JointStereo.convert (per-band conversion,
inlined at 0x437150, using CHSW_TABLE/CHWT_TABLE) and split_bits (bit split, 0x4378c0).

Per QMF band b the encoder decides a 'conversion' d:
    0  : channel 1 (right) dominates          -> matrix selector 0
    1  : channel 0 (left) dominates           -> matrix selector 1
    3  : comparable, loud                      -> matrix selector 3 (sum/difference)
   -3  : comparable, quiet (|d| = 3 is used)   -> matrix selector 3
For bands below `limit` (the bands in which the second sound unit carries data) the two
channels are matrixed into A = (L+R)/2 and a weighted difference B.  For the remaining bands
only A is coded and the decoder rebuilds L and R from A with a pair of weights (intensity
stereo); the weight index is one 4-bit value per frame.
Decisions and weights take effect with a 8-sample linear crossfade from the previous frame's
values, and are transmitted one frame later (the bitstream of frame n carries the values
decided in frame n-1, which describe the QMF samples that frame n's MDCT finishes).
"""
import struct
from lp4_fp import f32, bits
from lp4_tables import (CHSW_TABLE, CHWT_TABLE, JS_K_SIDE, JS_K_MATRIX, JS_K_HYST, JS_QUIET,
                        JS_RATIO_GAIN, SB_START, SB_WIDTH)


class JointStereo:
    def __init__(self, limit):
        self.limit = limit                 # number of matrixed bands (state[0xb64]); 1 for LP4
        self.dec_prev = [0] * 4            # |decision| of previous frame (state 0xb68), sent in the stream
        self.dec = [-3] * 4                # decision of this frame (state 0xb78)
        self.w_prev = 0                    # weight index of previous frame (0xb88), sent in the stream
        self.w = 15                        # weight index of this frame (0xb8c)
        self.r_prev = [1.0] * 4            # side scaling of previous frame (0xb90)
        self.r = [1.0] * 4                 # side scaling of this frame (0xba0)
        self.energy = [[(0.0, 0.0)] * 4 for _ in range(3)]  # 0xbb0: [cur, prev, prev2][band] = (L, R)

    # ------------------------------------------------------------------ decision
    def analyse(self, X0, X1):
        """Band energies (sum of |x|) and per-band decisions; weight index for the intensity bands."""
        acc0 = acc1 = f32(0.1)
        cur = []
        for b in range(4):
            s0 = 0.0
            s1 = 0.0
            col0 = X0[:, b].tolist()
            col1 = X1[:, b].tolist()
            for t in range(256):
                s0 += abs(col0[t])
                s1 += abs(col1[t])
            e0, e1 = f32(s0), f32(s1)
            cur.append((e0, e1))
            prev = self.dec[b]
            k0 = k1 = JS_K_SIDE
            if b < self.limit:                 # hysteresis: easier to stay in a one-sided mode
                k0 = k1 = JS_K_MATRIX
                if prev == 0:
                    k0 = JS_K_HYST
                elif prev == 1:
                    k1 = JS_K_HYST
            if k0 * s0 < s1:
                d = 0
            elif k1 * s1 < s0:
                d = 1
            elif s1 + s0 < JS_QUIET:
                d = -3
            else:
                d = 3
            self.dec_prev[b] = abs(prev)
            self.dec[b] = d
            if b >= self.limit and d in (3, -3):
                acc0 = acc0 * 16.0 + e0
                acc1 = acc1 * 16.0 + e1
        self.energy = [cur] + self.energy[:2]
        self.w_prev = self.w
        tot = f32(acc1 + acc0)
        if acc0 <= acc1:
            mn, side = acc0, 0
        else:
            mn, side = f32(acc1), 8
        w = (bits(f32(mn / tot * 14.0 + 12582912.0)) & 7) + side
        if w in (0, 8):
            w += 1
        self.w = w

    # ------------------------------------------------------------------ conversion
    def convert(self, b, X0, X1):
        """Apply the matrix (b < limit) or the intensity weighting (b >= limit) to band b in place."""
        cur = abs(self.dec[b])
        if b >= self.limit:
            wc, wp = self.w, self.w_prev
            if cur == 0:
                wc, wp = (wc ^ 8) + 8, (wp ^ 8) + 8
            elif cur == 1:
                wc, wp = wc + 8, wp + 8
            elif cur == 3:
                wc, wp = wc & 7, wp & 7
            tc, tp = CHWT_TABLE[wc], CHWT_TABLE[wp]
            step = (tp - tc) * 0.125
            for i, k in enumerate(range(-8, 0)):
                X0[i, b] = (float(X1[i, b]) + float(X0[i, b])) / (tc - k * step)
            inv = 1.0 / tc
            col0 = X0[8:, b].tolist()
            col1 = X1[8:, b].tolist()
            X0[8:, b] = [(c1 + c0) * inv for c0, c1 in zip(col0, col1)]
            return
        e0, e1 = self.energy[0][b]
        prv = self.dec_prev[b]
        c0, c1 = CHSW_TABLE[cur], CHSW_TABLE[cur + 1]
        p0, p1 = CHSW_TABLE[prv], CHSW_TABLE[prv + 1]
        r = 1.0
        if cur == 3:                           # side gain from the energy imbalance
            if not (e1 <= e0 * 3.94) and not (e0 * 256.0 <= e1):
                r = f32((e1 - e0) * JS_RATIO_GAIN / (e1 + e0))
            elif not (e0 <= e1 * 3.94) and not (e1 * 256.0 <= e0):
                r = f32((e0 - e1) * JS_RATIO_GAIN / (e1 + e0))
            elif e1 == 0.0 or e0 == 0.0:
                r = 1.0
            elif e1 * 256.0 <= e0 or e0 * 256.0 <= e1:
                r = f32(5.0 / 3.0)            # 0x3fd55555
        q = self.r_prev[b]
        if not (r - q < 0.05):                 # limit the change to 0.05 per frame
            r = f32(q + 0.05)
        if not (q - r < 0.05):
            r = f32(q - 0.05)
        self.r[b] = r
        c1r = c1 * r
        dc = f32((p0 - c0) * 0.125)
        dw = (p1 * q - c1r) * 0.125
        self.r_prev[b] = r
        for i, k in enumerate(range(-8, 0)):
            x0, x1 = float(X0[i, b]), float(X1[i, b])
            m = (x1 + x0) * 0.5
            X1[i, b] = (x0 - (c0 - dc * k) * m) / (c1r - k * dw)
            X0[i, b] = m
        inv = 1.0 / c1r
        col0 = X0[8:, b].tolist()
        col1 = X1[8:, b].tolist()
        ms = [(c1_ + c0_) * 0.5 for c0_, c1_ in zip(col0, col1)]
        X1[8:, b] = [(c0_ - m * c0) * inv for c0_, m in zip(col0, ms)]
        X0[8:, b] = ms


# ---------------------------------------------------------------------- split_bits
def _ftol(x):
    return int(x)              # x87 _ftol: truncation toward zero


def split_bits(spec0, spec1, nsub0, nsub1, bytes_per_ch, base0, dec_prev0, gaincount1_b0):
    """Bits granted to the first sound unit in joint-stereo mode (0x4378c0).

    spec0/spec1: float32[1024] spectra of the two units (A and B); nsub0/nsub1: bandwidth limits
    in subbands; base0: the static budget of unit 0 (bytes0*8 - 19); dec_prev0: |decision| of
    band 0 that is transmitted in this frame; gaincount1_b0: gain points of unit 1 in band 0.
    Returns the bit budget of unit 0; unit 1 later receives whatever bytes unit 0 leaves.
    """
    total = bytes_per_ch * 16 - 0x3b
    # weighted peak sums: unit 0 over subbands min(nsub0,19)-1 .. 0, unit 1 over 0 .. nsub1-1
    sum0 = 0.0
    for i in range(min(nsub0 - 1, 18), -1, -1):
        mx = 0.0
        for v in spec0[SB_START[i]:SB_START[i + 1]].tolist():
            if not (mx >= abs(v)):
                mx = abs(v)
        sum0 += SB_WIDTH[i] * mx
    sum0 = f32(sum0)
    sum1 = 0.0
    for i in range(nsub1):
        mx = 0.0
        for v in spec1[SB_START[i]:SB_START[i + 1]].tolist():
            if not (mx >= abs(v)):
                mx = abs(v)
        sum1 += SB_WIDTH[i] * mx
    t = f32(sum1 + sum1)
    if dec_prev0 != 3:
        t = f32(t + t)
    if t < 1.0 or (dec_prev0 != 3 and t < 4.0):
        return total
    spare = total - base0
    if gaincount1_b0 == 0 and ((t * 3.0 < sum0 and dec_prev0 == 3) or t * 6.0 < sum0):
        if not (t * 12.0 > sum0):
            x = spare
        else:
            if not (t * 10.0 >= sum0):
                spare = _ftol(spare * 1.5)
            elif not (t * 5.0 >= sum0):
                spare = _ftol(spare * 1.25)
            x = _ftol((sum0 - t) / (sum0 + sum0) * spare)
        return base0 + (x & ~7)
    if t >= sum0 or dec_prev0 != 3:
        return base0
    x = _ftol((sum0 - t) / (sum0 * 10.0) * spare)
    return base0 + (x & ~7)
