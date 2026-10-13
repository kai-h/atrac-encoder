"""Bit allocation, tonal component extraction and quantisation of the LP4 encoder.

Addresses: allocate 0x437b40 and its helpers 0x4385c0..0x4395d0, best_sf_tone 0x438bb0, bitcount
(an SSE kernel) 0x436740.  Integer logic is reproduced literally; the tone quantiser
is x87 code (doubles, float32 stores) and the spectral quantiser/bit counter is SSE code
(float32 multiply, then float32 add of 1.5*2**23 whose bit pattern holds the rounded value).

Terminology: 'subband' = one of the 32 quantisation units of the bitstream (SB_START), 'sel' =
quantiser selector (word length) 0..7, 'sf' = scale factor index 0..63 (2**((sf-15)/3)),
'group' = 4 consecutive spectral coefficients.
"""
import numpy as np
from lp4_fp import F32, MAGIC, f32, bits, from_bits, low16
from lp4_tables import (SB_START, SB_WIDTH, SF_BITS, QSCALE_BITS, QMASK, MIN_SFSUM, ZERO_LEN2,
                        STEP_COST, EST_SAVE, EST_BITS, EST_THR, TONE_PENALTY, TONE_SF_OFFSET,
                        SEL1_LEN, VLC, TONE_BEST_INIT)

_LEN = {key: np.array([ln for _, ln in tab]) for key, tab in VLC.items()}
_SEL1_LEN = np.array(SEL1_LEN)


# ---------------------------------------------------------------------- basic helpers
def sf_index4(v):
    """Scale factor index of a group of 4 float32 values (0x438760 / 0x4365e0): the smallest
    index whose value 2**((i-15)/3) is not below the peak magnitude, computed from the bit
    pattern of the peak (exponent * 3 plus a third-octave from the mantissa)."""
    u = max((bits(float(x)) << 1) & 0xffffffff for x in v)
    e, m = u >> 24, u & 0xffffff
    i = 3 * e - 0x16c
    if m > 0x965fe9:
        i += 1
    if m < 0x428a30:
        i -= 1
    return i if 0 <= i <= 63 else 0


def sf_index_all(spec):
    """sf_index4 for all 256 groups of a 1024-coefficient spectrum (vectorised)."""
    u = (np.abs(spec).astype(F32).view(np.uint32).reshape(256, 4).max(1).astype(np.int64) << 1)
    e, m = u >> 24, u & 0xffffff
    i = 3 * e - 0x16c + (m > 0x965fe9) - (m < 0x428a30)
    return [int(x) if 0 <= x <= 63 else 0 for x in i]


def qscale(sf, sel):
    """Quantiser scale MAX_QUANT[sel] / 2**((sf-15)/3) as float32, built like the binary does:
    a mantissa table indexed by 3*sel - sf%3 and an exponent decrement of sf//3."""
    u = (sf * 0x2b0000) & 0xffffffff
    return from_bits(QSCALE_BITS[((u >> 23) + sel) * 3 - sf] - (u & 0x7f800000))


def quantise(x, sel, sf):
    """SSE quantiser: rounded integers (as int32 bit patterns - 0x4b400000) of x * scale."""
    v = (np.asarray(x, F32) * F32(qscale(sf, sel))).astype(F32) + F32(MAGIC)
    return v.view(np.int32).astype(np.int64) - 0x4b400000, v


def bitcount(sel, sf, n, x):
    """Bits of one subband: 6 (scale factor) + VLC lengths of the quantised values (0x436740).
    Values outside the selector's range wrap through the index mask exactly as in Sony's code."""
    q, v = quantise(x[:n], sel, sf)
    if sel == 1:
        nz = (v != F32(MAGIC)).reshape(-1, 4)
        idx = nz[:, 0] + 2 * nz[:, 1] + 4 * nz[:, 2] + 8 * nz[:, 3]
        return 6 + int(_SEL1_LEN[idx].sum())
    return 6 + int(_LEN[(0, sel)][q & QMASK[sel]].sum())


# ---------------------------------------------------------------------- work state
class Comp:
    """Tonal component record (6 ints at S+0x540+24n): q[0..3] (masked), pos, sf."""
    __slots__ = ('r',)

    def __init__(self):
        self.r = [0] * 6

    q = property(lambda s: s.r)
    pos = property(lambda s: s.r[4], lambda s, v: s.r.__setitem__(4, v))
    sf = property(lambda s: s.r[5], lambda s, v: s.r.__setitem__(5, v))


class ToneGroup:
    """One tonal coding group: band flags, quantiser step, values per component - 1 and, per
    block of 64 coefficients, a list of at most 7 components."""
    def __init__(self):
        self.flag = [0] * 4
        self.qstep = 0
        self.nvm1 = 0
        self.lst = [[] for _ in range(16)]


class AllocState:
    """Sony's shared allocation work area (state+0x000..0x10ff and state+0xc80 spectrum)."""
    def __init__(self, spec):
        self.sel = [0] * 32
        self.sf = [0] * 32
        self.nsub = 0
        self.nbands = 0
        self.sel_mode = 0          # 0 = tonal values VLC coded, 1 = CLC
        self.ngroups = 0
        self.groups = [ToneGroup(), ToneGroup()]
        self.comps = [Comp() for _ in range(80)]
        self.spec = np.array(spec, F32)      # residual after tone extraction


class UnitInfo:
    """Per sound unit inputs to the allocator (fields of Sony's per-channel block)."""
    def __init__(self, budget, js, nsub_max, gain_counts, prev_count0, weight_sent=7):
        self.budget = budget              # +0x2d60 bits available for this unit
        self.js = js                      # +0x2d64 nonzero for the joint-stereo second unit
        self.nsub_max = nsub_max          # +0x2d5c bandwidth limit (subbands)
        self.gain_counts = gain_counts    # gain points per QMF band (this frame)
        self.prev_count0 = prev_count0    # +0x2d58 band-0 gain points of the previous frame
        self.weight_sent = weight_sent    # stereo weight index transmitted in this frame


def _clamp_sel(v, shift):
    v >>= shift
    return 1 if v < 1 else (7 if v > 7 else v)


# ---------------------------------------------------------------------- tones
def best_sf_tone(spec, pos, rec, qstep, n=4):
    """best_sf_tone (0x438bb0): try the 3 scale factors above the peak's index and keep the
    one with the smallest maximum rounding error (relative to the step), then move up in steps
    of 3 (one octave) while all quantised values are even."""
    best = TONE_BEST_INIT
    x = [float(v) for v in spec[pos:pos + 4]]
    lo = sf_index4(x)
    hi = min(lo + 3, 64) - 1
    for sf in range(hi, lo - 1, -1):
        scale = qscale(sf, qstep)
        err = 0.0
        acc = 0
        for k in range(n):
            p = scale * x[k]
            q = low16(f32(p + MAGIC))
            acc |= q
            d = abs(p - q)
            if not (err >= d):
                err = f32(d)
        if not (scale * best <= err):
            best = err / scale
            rec.sf = sf
            if (acc & 1) == 0:
                while rec.sf < 61:
                    acc >>= 1
                    rec.sf += 3
                    if acc & 1:
                        break


def quant_tone(spec, rec, qstep, n=4):
    """Quantise a tonal component at rec.pos (0x438af0); subtract it from the spectrum.
    Returns its VLC cost in bits (6 sf + 6 position + codes)."""
    best_sf_tone(spec, rec.pos, rec, qstep, n)
    scale = qscale(rec.sf, qstep)
    inv = f32(1.0 / scale)
    nb = 12
    lens = _LEN[(0, qstep)]
    for k in range(n):
        p = rec.pos + k
        x = float(spec[p])
        q = low16(f32(scale * x + MAGIC))
        idx = q & QMASK[qstep]
        spec[p] = x - q * inv
        rec.q[k] = idx
        nb += int(lens[idx])
    return nb


def restore_tone(spec, rec, nvm1, qstep):
    """Add a component back into the spectrum (0x438dc0); uses the unrounded reciprocal."""
    inv = 1.0 / qscale(rec.sf, qstep)
    nbits = VLC[(1, qstep)][0][1]
    for k in range(nvm1, -1, -1):
        v = rec.q[k]
        if v >= 1 << (nbits - 1):
            v -= 1 << nbits
        p = rec.pos + k
        spec[p] = v * inv + float(spec[p])


def find_peaks(spec, sfthr, start, end, maxn):
    """Positions with |x| >= SF value of sfthr, skipping 3 after each hit (0x4395d0).
    Returns the list, or None if more than maxn peaks exist."""
    thr2 = (SF_BITS[sfthr] * 2) & 0xffffffff
    out = []
    k = start
    while True:
        if ((bits(float(spec[k])) << 1) & 0xffffffff) >= thr2:
            if len(out) == maxn:
                return None
            out.append(k)
            k += 3
        k += 1
        if k >= end:
            return out


def shift_comp(rec, i):
    """Drop leading zero values of a component by moving it right (0x439540)."""
    e = 1
    while e < 4 and rec.q[e] == 0:
        e += 1
    np_ = rec.pos + e
    if np_ >= SB_START[i + 1]:
        return -1
    if np_ > 0x3fc:
        return 0
    r = rec.r
    ncopy = max(4 - e, 1)
    for c in range(ncopy):
        r[c] = r[c + e]
    for c in range(ncopy, 4):
        r[c] = 0
    rec.pos = np_
    return e


def shrink_values(S, n):
    """Reduce the number of values per component while the trailing value is zero for all
    components (0x439360).  Returns how many values were removed."""
    g = S.groups[0]
    while True:
        flag = False
        nv = g.nvm1
        for j in range(n - 1, -1, -1):
            rec = S.comps[j]
            if rec.q[nv] != 0:
                rec.pos += 1
                for c in range(nv):
                    rec.q[c] = rec.q[c + 1]
            if rec.q[nv - 1] != 0 and (rec.q[0] != 0 or (rec.pos & 0x3f) == 0x3f):
                flag = True
        g.nvm1 = nv - 1
        if flag:
            return 3 - g.nvm1


def split_comps(S, n):
    """Split every 4-value component whose last value is nonzero into two 2-value components
    (0x439400).  Returns the new component count, or 0 if not possible."""
    g = S.groups[0]
    cnt = [len(l) for l in g.lst]
    total = n
    for j in range(n):
        rec = S.comps[j]
        if rec.q[3] != 0:
            p = rec.pos + 2
            blk = p >> 6
            if g.flag[blk >> 2] == 0 or cnt[blk] >= 7 or total == 64 or p > 0x3fd:
                return 0
            cnt[blk] += 1
            total += 1
    total = n
    for j in range(n):
        rec = S.comps[j]
        if rec.q[3] != 0:
            new = S.comps[total]
            new.r = [rec.q[2], rec.q[3], 0, 0, rec.pos + 2, rec.sf]
            rec.q[2] = rec.q[3] = 0
            g.lst[new.pos >> 6].append(new)
            total += 1
    return total


def extract_tones(bits_left, P, S):
    """Tonal component extraction for the normal allocation path (0x438e60).  Peaks well above
    their subband's level are coded as 4-value components (quantiser step 3) and removed from
    the spectrum.  Marks affected subbands with a negative sf and reuses sel[] as the number of
    coefficients covered by tones.  Returns the bits spent on the tonal block."""
    S.groups[0] = g = ToneGroup()
    g.qstep, g.nvm1 = 3, 3
    hdr = 8
    ncomp = 0
    bits_tone = 0
    nzero = 0
    i = 8 if bits_left < 1100 else 0
    while i < S.nsub:
        if ncomp > 64 - (SB_WIDTH[i] >> 4) or 24 * ncomp + hdr > bits_left - 600:
            break
        sel = S.sel[i]
        if sel != 0:
            if P[i] < 0x500:
                sel = 0
            start = SB_START[i]
            band = 1 if SB_WIDTH[i] == 32 else start >> 8
            thr_sf = max(S.sf[i] - TONE_SF_OFFSET[sel + band * 8], 0)
            st = start
            if ncomp > 0 and S.comps[ncomp - 1].pos + 4 > start:
                S.sf[i] = -S.sf[i]
                ov = S.comps[ncomp - 1].pos - start
                S.sel[i] = ov + 4
                S.sel[i - 1] += ov
                st = S.comps[ncomp - 1].pos + 4
            peaks = find_peaks(S.spec, thr_sf, st, SB_START[i + 1], (SB_WIDTH[i] + 8) >> 4)
            if peaks:
                if S.sf[i] >= 0:
                    S.sf[i] = -S.sf[i]
                    S.sel[i] = 0
                for pk in peaks:
                    c = min(pk, 0x3fc)
                    nv = 4
                    if ncomp > 0 and S.comps[ncomp - 1].pos + 4 == c:
                        while True:
                            if abs(float(S.spec[c + 3])) > abs(float(S.spec[c - 1])):
                                break
                            c -= 1
                            nv -= 1
                            if nv <= 1:
                                break
                    blk = c >> 6
                    if len(g.lst[blk]) >= 7:
                        continue
                    rec = S.comps[ncomp]
                    ncomp += 1
                    rec.pos = c
                    nb = quant_tone(S.spec, rec, 3, 4)
                    if rec.q[0] == 0 and rec.pos <= 0x3fc:
                        if shift_comp(rec, i) < 0:
                            restore_tone(S.spec, rec, 3, 3)
                            ncomp -= 1
                            break
                        blk = rec.pos >> 6
                    bits_tone += nb
                    S.sel[i] += nv
                    if g.flag[blk >> 2] == 0:
                        g.flag[blk >> 2] = 1
                        hdr += 12
                    g.lst[blk].append(rec)
                    if rec.q[3] == 0:
                        nzero += 1
        i += 1
    if ncomp == 0:
        return 0
    if ncomp == 1 and (S.comps[0].pos < 128 or bits_left < 1100):
        restore_tone(S.spec, S.comps[0], 3, 3)
        return 0
    hdr += S.nbands
    S.ngroups = 1
    len0v = VLC[(0, g.qstep)][0][1]
    len0c = VLC[(1, g.qstep)][0][1]
    vlc, clc = bits_tone, 24 * ncomp
    do_shrink = nzero == ncomp
    if not do_shrink:
        nz = 2 * ncomp - 3 * nzero
        x12 = 12 * (ncomp - nzero)
        est_v = len0v * nz + x12 + bits_tone
        est_c = len0c * nz + x12 + 24 * ncomp
        if min(bits_tone, 24 * ncomp) > min(est_v, est_c):
            r = split_comps(S, ncomp)
            if r:
                ncomp = r
                vlc = len0v * r + est_v
                clc = len0c * r + est_c
                do_shrink = True
    if do_shrink:
        k = shrink_values(S, ncomp) * ncomp
        vlc -= len0v * k
        clc -= len0c * k
    if vlc > clc:
        S.sel_mode = 1
        best = clc
    else:
        S.sel_mode = 0
        best = vlc
    return hdr + best


def _replace_weakest(S, gi, blk, sf):
    """Make room in a full block list: remove (restore) the component with the lowest sf below
    `sf` and return its record for reuse (0x438ce0)."""
    grp = S.groups[gi]
    best = -1
    lst = grp.lst[blk]
    for j in range(7):
        if lst[j].sf < sf:
            best = j
            sf = lst[j].sf
    if best < 0:
        return None
    rec = lst.pop(best)
    restore_tone(S.spec, rec, grp.nvm1, grp.qstep)
    return rec


def all_tonal(bits_left, P, sfidx, S):
    """Allocation for frames with plenty of bits (0x4387c0): every coefficient of magnitude
    >= 2.0 becomes a tonal component (CLC coded, step 5 for sf <= 30, step 7 above) until the
    spectrum is exhausted or the bits run out; the remaining subbands get selector 1."""
    need = SB_START[S.nsub] >> 8
    cost = 22
    if S.nbands < need:
        cost = (need - S.nbands) * 3 + 22
        S.nbands = need
    S.groups = [ToneGroup(), ToneGroup()]
    S.groups[0].qstep, S.groups[1].qstep = 5, 7
    S.groups[0].nvm1 = S.groups[1].nvm1 = 3
    S.sel_mode = 1
    S.ngroups = 2
    ncomp = 0
    limit = bits_left - 200
    thr2 = (SF_BITS[18] * 2) & 0xffffffff          # 2.0
    end = SB_START[S.nsub]
    done = False
    while not done and end > 0:
        start_n = ncomp
        k = 0
        while k < end:
            if ((bits(float(S.spec[k])) << 1) & 0xffffffff) >= thr2:
                rec = S.comps[ncomp]
                ncomp += 1
                rec.pos = min(k, 0x3fc)
                sf = sf_index4(S.spec[rec.pos:rec.pos + 4])
                gi = 1 if sf > 30 else 0
                blk = k >> 6
                if len(S.groups[gi].lst[blk]) == 7:
                    if len(S.groups[1 - gi].lst[blk]) == 7:
                        ncomp -= 1
                        rec = _replace_weakest(S, gi, blk, sf)
                        if rec is None:
                            k += 1
                            continue
                    else:
                        gi = 1 - gi
                grp = S.groups[gi]
                if grp.flag[blk >> 2] == 0:
                    grp.flag[blk >> 2] = 1
                    cost += 12
                grp.lst[blk].append(rec)
                quant_tone(S.spec, rec, grp.qstep, 4)
                cost += 28 + 8 * gi
                if cost > limit or ncomp >= 64:
                    done = True
                    break
                k += 3
            k += 1
        if start_n == ncomp:
            break
    for i in range(S.nsub):
        mx = 0
        for gg in range(SB_START[i] >> 2, SB_START[i + 1] >> 2):
            v = sf_index4(S.spec[4 * gg:4 * gg + 4])
            sfidx[gg] = v
            mx = max(mx, v)
        S.sf[i] = mx
        S.sel[i] = 1 if mx else 0
        P[i] = mx * 341
    return cost


# ---------------------------------------------------------------------- allocation
def est_bits(i, sfidx, S):
    """Bit estimate of subband i in tenths of a bit (0x4385c0)."""
    sel = S.sel[i]
    thr = S.sf[i] - EST_THR[sel]
    b = EST_BITS[sel] * SB_WIDTH[i] + 60
    g0 = SB_START[i] >> 2
    for g in range(g0, g0 + (SB_WIDTH[i] >> 2)):
        if sfidx[g] < thr:
            b -= EST_SAVE[sel]
    return b


def gain_header_bits(U, S):
    """Bands needed for the coded subbands and the bits of their gain data (0x438630)."""
    nb = (SB_START[S.nsub] + 0xff) >> 8
    S.nbands = nb
    return (sum(U.gain_counts[:nb]) * 3 + nb) * 3


def drop_empty_bands(U, S):
    """Drop trailing QMF bands without gain points or tones from the band count (0x438690).
    Not done for the joint-stereo unit while intensity weighting is active."""
    if U.js and (U.weight_sent & 7) != 7:
        return 0
    n = S.nbands
    k = n - 1
    while k > 0 and U.gain_counts[k] == 0:
        if any(S.groups[j].flag[k] for j in range(S.ngroups)):
            break
        k -= 1
    S.nbands = k + 1
    return (S.ngroups + 3) * (n - (k + 1))


def allocate(U, spec):
    """allocate: choose nsub, selectors, scale factors and tonal components for one unit.
    Returns the AllocState (S) holding everything the packer needs."""
    S = AllocState(spec)
    budget = U.budget
    nsub = U.nsub_max
    sfidx = sf_index_all(S.spec)
    cnt7 = 0
    cum = 0
    top = 0
    sums = [0] * 32
    P = [0] * 32
    for i in range(32):
        s = 0
        mx = 0
        for g in range(SB_START[i] >> 2, SB_START[i + 1] >> 2):
            v = sfidx[g]
            s += v
            mx = max(mx, v)
            if v > 7:
                cnt7 += 1
        sums[i] = s
        if top < mx:
            top = mx
            if not U.js:          # a new loudest subband may extend the bandwidth
                if (SB_START[i] & 0x80) == 0 or i >= 22:
                    ok = cnt7 * 8 < budget
                else:
                    ok = cnt7 * 16 < budget
                if ok:
                    t = 28 if i <= 26 else min(i + 2, 32)
                    nsub = max(nsub, t)
        elif s < MIN_SFSUM[i] or mx < 3:
            s = 0
            mx = 0
        S.sf[i] = mx
        P[i] = cnt7
        cum += s
        S.sel[i] = cum
    S.nsub = nsub
    shift = 11 if budget * 2 < P[nsub - 1] * 12 else 10
    bits_left = budget - 3 * nsub - gain_header_bits(U, S)
    if bits_left > P[nsub - 1] * 32:
        bits_left -= all_tonal(bits_left, P, sfidx, S)
    else:
        for i in range(nsub):
            m = S.sf[i]
            v = 0
            if m >= 3:
                m = m * 256 - cum
                v = _clamp_sel(m, shift)
            S.sel[i] = v
            P[i] = m
        S.ngroups = 0
        if not U.js:
            bits_left -= extract_tones(bits_left, P, S)
        for i in range(nsub):
            m = S.sf[i]
            if m < 0:                         # subband lost tonal components: re-evaluate
                mx = 0
                s = 0
                for g in range(SB_START[i] >> 2, SB_START[i + 1] >> 2):
                    v = sf_index4(S.spec[4 * g:4 * g + 4])
                    sfidx[g] = v
                    mx = max(mx, v)
                    s += v
                S.sf[i] = mx
                if S.sel[i] == 0:
                    P[i] += 256
                else:
                    k = (S.sel[i] - 1) >> 1
                    if k < 0 or k > 4:
                        k = 4
                    d = TONE_PENALTY[k] * (mx + m)
                    if budget < 1100:
                        d *= 2
                    elif SB_WIDTH[i] <= 16 and P[i] > 0x1800:
                        d2 = 2 * d
                        d = d2 if d2 < P[i] - 0x400 else P[i] - 0x400
                    if i > 1:
                        P[i - 1] -= d >> 3
                    if i < 31:
                        P[i + 1] -= d >> 2
                    P[i] -= d
                v = 0
                if mx > 7 and s > MIN_SFSUM[i] and (P[i] >= 0x400 or budget > 0x708):
                    v = _clamp_sel(P[i], shift)
                S.sel[i] = v
            elif shift != 10:
                if SB_WIDTH[i] * m < sums[i] * 5:
                    if m >= 62:
                        continue
                    S.sf[i] = m + 1
                if S.sf[i] < 62:
                    S.sf[i] += 1
            elif SB_WIDTH[i] * m > sums[i] * 8:
                P[i] += 256
    # drop isolated subbands far below their neighbours, estimate the bits of the rest
    est = 0
    w7 = 0
    wsum = 0
    for i in range(nsub):
        if S.sel[i] == 0:
            continue
        if (i > 1 and P[i] < P[i - 1] - 0xf00) or (i < nsub - 1 and P[i] < P[i + 1] - 0x1400):
            S.sel[i] = 0
        else:
            wsum += SB_WIDTH[i]
            if S.sel[i] == 7:
                w7 += SB_WIDTH[i]
            est += est_bits(i, sfidx, S)
    k = nsub - 1
    while k > 0 and S.sel[k] == 0:
        bits_left += 3
        k -= 1
    nsub = k + 1
    S.nsub = nsub
    bits_left += drop_empty_bands(U, S)
    b10 = bits_left * 10
    off = 0
    width = None
    if b10 <= est:
        width = wsum
    elif wsum - w7 != 0:
        width = wsum - w7
    if width is not None:
        t = b10 if bits_left < 600 else bits_left * 9
        diff = t - est
        num = diff * 1024
        off = abs(num) // (width * 10)
        if (num < 0) != (width < 0):
            off = -off
        if diff < 0 and shift == 10:
            off += off >> 3
            P[0] -= off
            P[1] -= (off * 188) >> 8
            for k in range(2, 8):
                P[k] -= (off * 97) >> 8
            for k in range(8, 18):
                P[k] -= (off * 74) >> 8
            for k in range(18):
                if P[k] + off < 0 and S.sf[k] < 63:
                    S.sf[k] += 1
            for k in range(18, nsub):
                if P[k] + off < 0x400:
                    S.sf[k] = min(63, S.sf[k] - ((P[k] + off - 0x400) >> 10))
    for i in range(nsub - 1, -1, -1):
        if S.sel[i] != 0:
            S.sel[i] = _clamp_sel(P[i] + off, shift)
            b10 -= est_bits(i, sfidx, S)
            P[i] = b10
    # exact bit counts, local selector adjustment
    used = 0
    order = [-1] * 32
    for i in range(nsub):
        sel = S.sel[i]
        if sel == 0:
            continue
        if SB_WIDTH[i] > (bits_left - used) * 2 - 16:
            S.sel[i] = 0
            continue
        d = P[i] - used * 10
        if sel < 7 and d > 4000:
            S.sel[i] = sel + 1
        elif i != 0 and d < 300:
            if sel > 1:
                S.sel[i] = sel - 1
            elif S.sf[i] != 63:
                S.sf[i] += 1
        b = bitcount(S.sel[i], S.sf[i], SB_WIDTH[i], S.spec[SB_START[i]:])
        P[i] = b
        if b * 2 == ZERO_LEN2[S.sel[i]] * SB_WIDTH[i] + 12:      # everything quantised to 0
            S.sel[i] = 0
        else:
            used += b
            order[i] = S.sf[i]
    k = nsub - 1
    while k > 0 and S.sel[k] == 0:
        bits_left += 3
        k -= 1
    slack = bits_left - used
    nsub = k + 1
    S.nsub = nsub
    spent = 0
    if U.prev_count0 + U.gain_counts[0] != 0 and not U.js:
        for k in range(5):                    # transient frames: favour the lowest subbands
            if 0 <= order[k] <= 0x40:
                order[k] = 0x40
    # final pass in order of decreasing scale factor: fit the budget, spend leftovers
    while True:
        best = -1
        idx = 0
        for k in range(nsub):
            if best < order[k]:
                best = order[k]
                idx = k
        if best < 0:
            return S
        b = P[idx]
        order[idx] = -1
        w = SB_WIDTH[idx]
        x = S.spec[SB_START[idx]:]
        if spent + b > bits_left:
            zthr = ZERO_LEN2[S.sel[idx]] * (w >> 1) + 6
            k = nsub - 1
            while k > idx:
                if S.sel[k] != 0 and order[k] == -1:
                    break
                S.sel[k] = 0
                bits_left += 3
                order[k] = -1
                k -= 1
            nsub = k + 1
            S.nsub = nsub
            while spent + b > bits_left:
                S.sf[idx] += 1
                if S.sf[idx] >= 64 or b <= zthr:
                    break
                b = bitcount(S.sel[idx], S.sf[idx], w, x)
        elif S.sel[idx] < 7 and STEP_COST[S.sel[idx]] * (w >> 1) < bits_left - spent - 7:
            ns = S.sel[idx] + 1
            if slack > w * 4 and ns < 7:
                ns += 1
            b2 = bitcount(ns, S.sf[idx], w, x)
            if b2 - b <= (slack if slack > 0 else 0):
                S.sel[idx] = ns
                slack += b - b2
                b = b2
        spent += b
        if spent > bits_left:
            spent -= b
            slack += b
            S.sel[idx] = 0
