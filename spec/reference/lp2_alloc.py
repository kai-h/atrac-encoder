"""LP2 encoder core: tonal component extraction, bit allocation and quantisation of one sound
unit (0x43d1f0 and its helpers).

Inputs per channel and frame: the gain-compensated spectrum (quantised), the plain spectrum
(used for decisions) and the gain-control data.  Output: a SoundUnit that uses exactly the
number of bits it reports, never more than the 1536-bit budget.
"""
import collections
import numpy as np
from lp2_tables import (SF, SB_START, SB_SIZE, MAXQ, T1, T2, TONE_NVALS, TONE_SEL, TONE_LEN_CODE,
                        TONE_MODE, EST_BITS_CLC, MASK_SLOPE_NORMAL, MASK_SLOPE_TONAL,
                        offset_divisor, offset_step_weight, dead_zone_limit,
                        NSUB_DEFAULT, NBANDS_DEFAULT, UNIT_BITS)
from lp2_bitstream import (Component, ToneGroup, SoundUnit, mantissa_bits, gain_bits,
                           packed_bits)

f32, f64 = np.float32, np.float64
STATS = collections.Counter()      # branch statistics (informational only)


# ================================================================ scale factors
def sf_index(v):
    """sf_index (0x43eb70): 5-step binary search, result = smallest index whose
    scale factor exceeds v (0 if v < SF[0], 63 at most)."""
    v = f32(v)
    if not (SF[0] <= v):
        return 0
    e, c = 32, 16
    while c > 0:
        e = e + c if SF[e] <= v else e - c
        c >>= 1
    if e < 63 and SF[e] <= v:
        e += 1
    return e


def sf_index_tone(v):
    """sf_index_tone (0x43f0d0): first i with |v| < SF[i], else 63."""
    v = abs(f32(v))
    for i in range(63):
        if v < SF[i]:
            return i
    return 63


def group_sf(spec):
    """Scale-factor index of every 4-line group (256 groups)."""
    m = np.abs(np.asarray(spec, f32).reshape(256, 4)).max(axis=1)
    return [sf_index(x) for x in m]


def subband_sf(G, nsub):
    """Per subband, the largest group index."""
    return [max(G[SB_START[i] // 4:SB_START[i + 1] // 4]) for i in range(nsub)]


def abs_max(a):
    return f32(np.abs(np.asarray(a, f32)).max())


# ================================================================ quantisation
def dead_zone(e, sel):
    """dead_zone (0x43ebd0): magnitude below which a normalised line is quantised to 0."""
    if e == 0:
        return f32(0.0)
    return f32(f64(T1[min(max(e, 0), 12)]) / f64(T2[sel]))


def quantise_subband(mode, sel, thr, x):
    """Quantise normalised lines x with selector sel; returns (bits incl. 6-bit SF, mantissas).
    q = floor(x*(maxq+0.5) + 0.5) (computed as trunc(x*s + 31.5) - 31), clipped; |x| < thr -> 0."""
    if sel == 0:
        return 0, None
    mq = MAXQ[sel]
    x = np.asarray(x, f32)
    q = np.trunc(f64(f32(mq + 0.5)) * x.astype(f64) + 31.5).astype(np.int64) - 31
    q = np.clip(q, -mq, mq)
    q[~(thr <= np.abs(x))] = 0
    return mantissa_bits(mode, sel, q) + 6, q


def tone_mantissa(x, sfval, mq):
    """Tonal mantissa (0x43f790)."""
    v = int(np.trunc((f64(mq) + 0.5) * (f64(f32(x)) / f64(sfval)) + 31.5)) - 31
    return max(-mq, min(mq, v))


# ================================================================ tonal components
def quantise_component(spec, c):
    """Choose the component's scale factor and mantissas; returns its bit cost (incl. 12 bits
    for scale factor and position)."""
    n = TONE_NVALS[c.len_code]
    c.sfi = sf_index_tone(abs_max(spec[c.pos:c.pos + n]))
    sfv, mq = SF[c.sfi], MAXQ[c.sel]
    c.mant = [tone_mantissa(spec[c.pos + k], sfv, mq) if c.pos + k <= 1023 else 0 for k in range(n)]
    return mantissa_bits(c.mode, c.sel, c.mant) + 12


def subtract_component(spec, c):
    """Subtract the quantised component from the spectrum (in place)."""
    step = f32(f64(SF[c.sfi]) / (f64(MAXQ[c.sel]) + 0.5))
    for k in range(TONE_NVALS[c.len_code]):
        p = c.pos + k
        if p <= 0x3ff:
            spec[p] = f32(f64(spec[p]) - f64(c.mant[k]) * f64(step))


def _new_component(G, g, GA):
    c = Component()
    c.pos, c.sfi, c.len_code, c.sel, c.mode = 4 * g, GA[g], G.len_code, G.sel, G.mode
    return c


def single_tone_class(GA):
    """single_tone_class (0x43f110): 1 or 2 if the strongest 4-line group stands >= 17 scale
    steps (about 34 dB) above everything except +-2 groups around it and its QMF alias image."""
    p = 0
    for g in range(256):
        if GA[p] < GA[g]:
            p = g
    r, blk = p & 63, p >> 6
    bands = [blk, blk - 1 if r < 32 else blk + 1]     # peak band and nearest neighbour band
    centre = [r, 63 - r]                              # alias image is mirrored in the neighbour
    mx = 0
    for b in range(4):
        e = bands.index(b) if b in bands else -1
        base = 64 * b
        if e < 0:
            mx = max([mx] + GA[base:base + 64])
            continue
        u = centre[e]
        mx = max([mx] + GA[base:base + max(u - 2, 0)])
        if u + 2 < 63:
            mx = max([mx] + GA[base + u + 3:base + 64])
    v = GA[p]
    if mx < max(v - 17, 9):
        if v >= 36: return 2
        if v >= 21: return 1
    return 0


def single_tone_path(budget, ngroups_tone, width, peak, nbands, spec, GA, u):
    """single_tone_path (0x43f270): code the lines around a dominant peak (and its alias
    image in the neighbouring QMF band) as tonal components, in `ngroups_tone` groups."""
    blk, r = peak >> 6, peak & 63
    bands = [blk, blk - 1 if r < 32 else blk + 1]
    centre = [r, 63 - r]
    bits = (nbands + 6) * ngroups_tone
    go_on = True
    for _ in range(ngroups_tone):
        G = ToneGroup(TONE_SEL, TONE_LEN_CODE, TONE_MODE); u.groups.append(G)
        if not go_on:
            continue
        for e in (0, 1):
            b = bands[e]
            if b < 0 or b >= nbands:
                continue
            if not G.flags[b]:
                G.flags[b] = 1; bits += 12
            if bits >= budget:
                G.flags[b] = 0; bits -= 12; go_on = True
                continue
            c0 = centre[e] + 64 * b
            for g in range(max(c0 - width, 0), min(c0 + width, 255) + 1):
                cell = g >> 4
                if len(G.cells[cell]) >= 7:
                    continue
                c = _new_component(G, g, GA)
                cb = quantise_component(spec, c)
                if bits + cb > budget:
                    go_on = False
                    continue
                u.ncomp += 1; G.cells[cell].append(c)
                subtract_component(spec, c)
                GA[g] = sf_index(abs_max(spec[4 * g:4 * g + 4]))
                bits += cb
                STATS['single_tone_components'] += 1
    return bits


def multitone_path(budget, ngroups, nbands, mean, GA, GB, spec, u):
    """multitone_path (0x43ec90): up to two passes; every 4-line group whose plain-spectrum
    scale factor exceeds mean + 32 becomes a component (max 7 per cell, 64 in total)."""
    bits = 0; stop = False
    thr = f64(f32(f64(mean) + 32.0))
    for _ in range(2):
        if stop:
            continue
        cnt = [0] * 16; nc = 0
        for g in range(ngroups):
            if f64(GB[g]) > thr and cnt[g >> 4] < 7 and u.ncomp + nc < 64:
                cnt[g >> 4] += 1; nc += 1
        if nc <= 0:
            continue
        if nbands + bits + 6 > budget:
            stop = True; continue
        bits += nbands + 6
        G = ToneGroup(TONE_SEL, TONE_LEN_CODE, TONE_MODE); u.groups.append(G)
        for g in range(ngroups):
            cell, band = g >> 4, g >> 6
            if not (f64(GB[g]) > thr and len(G.cells[cell]) < 7 and u.ncomp < 64):
                continue
            if not G.flags[band]:
                if bits + 12 <= budget:
                    bits += 12; G.flags[band] = 1
                else:
                    stop = True
            if G.flags[band]:
                if bits + 4 * EST_BITS_CLC[G.sel] + 12 > budget:
                    stop = True
                    continue
                c = _new_component(G, g, GA)
                u.ncomp += 1; G.cells[cell].append(c)
                bits += quantise_component(spec, c)
                subtract_component(spec, c)
                s = sf_index(abs_max(spec[4 * g:4 * g + 4]))
                GB[g] += s - GA[g]; GA[g] = s
                STATS['multitone_components'] += 1
    return bits


# ================================================================ word-length decisions
def is_attack(g):
    """is_attack (0x43e760): gain levels (including the neutral level 4) span more than 2."""
    levs = [4] + [lev for lev, loc in g.pts]
    return max(levs) - min(levs) > 2


def word_lengths(slope, offs, SB, n, wmin=1, wmax=7):
    """word_lengths (0x43e7a0): word lengths from offsets, then zero the subbands that are
    absolutely quiet or masked by one of the 8 subbands below them.  Returns (W, quiet_thr)."""
    W = [min(max(int(np.trunc(f64(offs[i]) + 0.5)), wmin), wmax) for i in range(n)]
    mx = max(SB[1:n]) if n > 1 else 0
    mx = max(mx, 0)
    quiet = 6 if (n > 1 and mx >= 30) else max(mx // 6, 1)
    for i in range(n):
        if SB[i] < quiet:
            W[i] = 0
        for c in range(8):
            j = i - 1 - c
            if j >= 0 and SB[i] < SB[j] - slope[c]:
                W[i] = 0
    return W, quiet


def order_descending(keys):
    """order_descending (0x43f590): Shell sort (gaps ..., 13, 4, 1), descending; returns indices."""
    n = len(keys); K = list(keys); order = list(range(n))
    h = 1
    while n >= 1 and h <= n:
        h = 3 * h + 1
    while True:
        h //= 3
        if h < 1:
            return order
        if h >= n:
            continue
        for t in range(n - h):
            v = K[h + t]; j = t
            while j >= 0 and K[j] < v:
                K[j + h] = K[j]
                order[j + h], order[j] = order[j], order[j + h]
                j -= h
            K[j + h] = v


# ================================================================ encode_unit
def encode_unit(spec_gc, spec_plain, gain, prev_gain, bits=UNIT_BITS,
                  nsub=NSUB_DEFAULT, nbands=NBANDS_DEFAULT):
    """Build the SoundUnit of one channel for one frame."""
    A = np.array(spec_gc, f32)                # modified: tones removed, then normalised
    B = np.asarray(spec_plain, f32)
    u = SoundUnit(); u.gain = gain; u.ncomp = 0
    Q = np.zeros(1024, np.int64)
    GA, GB = group_sf(A), group_sf(B)
    peak = 0
    for g in range(256):
        if GB[peak] < GB[g]:
            peak = g

    # ---- 1. tonal components and the spectrum bit budget
    tone_class = single_tone_class(GA)
    tonal = False
    tone_bits = 0
    if tone_class > 0 and nsub > 25:
        tonal = True
        STATS['single_tone_frames'] += 1
        STATS['single_tone_class_%d' % tone_class] += 1
        if peak >= 160:
            nsub, nbands = 32, 4
            STATS['single_tone_4bands'] += 1
        elif peak >= 96:
            nsub, nbands = 30, 3
            STATS['single_tone_30sub'] += 1
        budget = bits - 3 * nsub - 8 - gain_bits(gain, nbands) - 13
        tone_bits = single_tone_path(budget, tone_class, (bits > 1200) + 1, peak, nbands,
                                         A, GA, u)
        budget -= tone_bits
    else:
        budget = bits - 3 * nsub - 8 - gain_bits(gain, nbands) - 13
    ngroups = SB_START[nsub] // 4
    mean = f32(f64(sum(GB[:ngroups])) / f64(ngroups))
    attack = -1                               # 0: attack now, 1: attack in previous frame
    if any(is_attack(gain[b]) for b in range(nbands)):
        attack = 0
    elif any(is_attack(prev_gain[b]) for b in range(nbands)):
        attack = 1
    slope = MASK_SLOPE_TONAL if tonal else MASK_SLOPE_NORMAL
    multi_bits = 0
    if tone_bits == 0:
        multi_bits = multitone_path(budget, ngroups, nbands, mean, GA, GB, A, u)
    budget -= multi_bits
    if not u.groups:
        budget += 2                           # no 2-bit tonal mode selector
    SA, SB = subband_sf(GA, nsub), subband_sf(GB, nsub)

    # ---- 2. initial word lengths from the scale-factor profile
    off = [f32((f64(SB[i]) - f64(mean)) * offset_divisor(i)) for i in range(nsub)]
    if attack != -1:
        STATS['attack_frames'] += 1
        for i in range(8):
            off[i] = f32(f64(off[i]) + 0.7)
        for i in range(8, 18):
            off[i] = f32(f64(off[i]) + 0.5)
        off[0] = max(off[0], f32(6.0))
    else:
        off[0] = max(off[0], f32(6.0))
        for i in range(1, 4):
            off[i] = max(off[i], f32(3.0))
    W, quiet = word_lengths(slope, off, SB, nsub)
    fixed = [0] * nsub                        # 1 = may not be raised later
    for i in range(nsub):
        if SA[i] < quiet:
            fixed[i] = 1; W[i] = 0
        if multi_bits == 0 and (W[i] == 0 or f64(SA[i]) < f64(mean)):
            fixed[i] = 1
    W0 = list(W)                              # upper limit for later increments

    # ---- 3. normalise each subband by its scale factor
    for i in range(nsub):
        inv = f32(1.0 / f64(SF[SA[i]]))
        a, b = SB_START[i], SB_START[i + 1]
        A[a:b] = (f64(inv) * A[a:b].astype(f64)).astype(f32)
    elim = dead_zone_limit(f64(budget) / f64(SB_START[nsub] + 256))
    STATS['elim_%d' % elim] += 1

    E = [0] * nsub                            # dead-zone index per subband
    sbits = [0] * nsub                        # bits per subband (mantissas + 6)

    def quant(i, sel, e):
        a, b = SB_START[i], SB_START[i + 1]
        return quantise_subband(0, sel, dead_zone(e, sel), A[a:b])

    def count_bits(flags):
        tot = 0
        for i in range(nsub):
            if flags[i]:
                sbits[i], q = quant(i, W[i], E[i])
                if q is not None:
                    Q[SB_START[i]:SB_START[i + 1]] = q
            tot += sbits[i]
        return tot

    total = count_bits([1] * nsub)

    # ---- 4. global offset search (bisection-like, 15 iterations max)
    offset, step = f32(2.0), f32(4.0)
    offs = list(off)
    for it in range(15):
        if offset > 0:
            offs = [f32(f64(offset) + f64(off[i])) for i in range(nsub)]
        else:
            offs = []
            for i in range(nsub):
                w = offset_step_weight(i)
                offs.append(f32(f64(off[i]) + f64(offset) * w) if w is not None
                            else f32(f64(offset) + f64(off[i])))
        W_old = list(W)
        W, _ = word_lengths(slope, offs, SB, nsub)
        W = [0 if W0[i] == 0 else W[i] for i in range(nsub)]
        total = count_bits([int(W[i] != W_old[i]) for i in range(nsub)])
        if total > budget:
            offset = f32(f64(offset) - f64(step))
        else:
            if it >= 6:
                break
            offset = f32(f64(step) + f64(offset)) if total < budget else f32(f64(offset) - f64(step))
        step = f32(f64(step) * (0.5 if it < 7 else 1.5))

    # ---- 5. if still too large: widen dead zones, then lower word lengths
    for _ in range(4):
        i = nsub - 1
        while i >= 0 and total > budget:
            if W[i] > 0 and E[i] < elim:
                E[i] += 1
                total = count_bits([1] * nsub)
                STATS['deadzone_steps'] += 1
            i -= 1
    if total > budget:
        STATS['wl_reduction_frames'] += 1
        keys = ([(SB[i] + 1) * 32 - i for i in range(nsub)] if tonal
                else [SB[i] - i // 2 for i in range(nsub)])
        order = order_descending(keys)
        k = nsub
        while total > budget and k > 0:
            k -= 1
            i = order[k]                      # lowest-priority subband first
            while W[i] > 0 and total > budget:
                W[i] -= 1
                total = count_bits([int(j == i) for j in range(nsub)])

    # ---- 6. fix the coded subband count; spend the remaining bits
    nused = nsub
    while W[nused - 1] == 0 and nused > 1:
        nused -= 1
    used = bits - (budget + 3 * (nsub - nused) - total)    # all bits of the unit so far
    sel = W[:nused]
    sfi = SA[:nused]
    margin = [f32(f64(offs[i]) - f64(W[i])) if W[i] else f32(-9.99) for i in range(nused)]
    lo = hi = margin[0]
    for m in margin[1:]:
        if not (lo <= m): lo = m
        if m > hi: hi = m
    span = f64(hi) - f64(lo)
    order, rank = [], [-1] * nused            # subbands by margin, best first (10 classes)
    for k in range(10):
        t = f32(f64(hi) - f64(k + 1) * span * 0.1)
        for i in range(nused - 1, -1, -1):
            if not (margin[i] < t) and rank[i] == -1:
                rank[i] = len(order); order.append(i)

    def commit(i, b_new, q):
        nonlocal total, used
        total += b_new - sbits[i]; used += b_new - sbits[i]; sbits[i] = b_new
        Q[SB_START[i]:SB_START[i + 1]] = q

    # 6a. one step per subband in margin order: raise W (up to W0) or narrow the dead zone
    for i in order:
        if fixed[i] or not (0 < sel[i] < 7):
            continue
        s0, e0 = sel[i], E[i]
        if s0 < W0[i]:
            sel[i] += 1; STATS['refine_wl'] += 1
        elif E[i] > 0:
            E[i] -= 1; STATS['refine_deadzone'] += 1
        else:
            if used > bits - 8:
                break
            continue
        b_new, q = quant(i, sel[i], E[i])
        if used + b_new - sbits[i] > bits:
            sel[i], E[i] = s0, e0
            continue
        commit(i, b_new, q)
        STATS['refine_steps'] += 1
        if used > bits - 8:
            break
    # 6b. remove remaining dead zones where affordable
    for i in range(nused):
        if E[i] > 0 and sel[i] > 0:
            b_new, q = quant(i, sel[i], 0)
            if used + b_new - sbits[i] <= bits:
                commit(i, b_new, q); E[i] = 0; STATS['deadzone_removed'] += 1

    def try_raise(i):
        sel[i] += 1
        b_new, q = quant(i, sel[i], E[i])
        if used + b_new - sbits[i] > bits:
            sel[i] -= 1
            return False
        commit(i, b_new, q)
        STATS['raise_steps'] += 1
        return True

    # 6c. raise word lengths towards W0 (7 passes)
    for _ in range(7):
        if used > bits - 2: break
        for i in range(nused):
            if used > bits - 2: break
            if 0 < sel[i] < W0[i] and not fixed[i]:
                if not try_raise(i):
                    fixed[i] = 1
    # 6d. raise low word lengths in scale-factor order (4 passes, limit 3..6)
    if nused > 1:
        prio = order_descending([(SB[i] + 1) * 32 - i for i in range(nused)])
        for p in range(4):
            if used > bits - 2: break
            for i in prio:
                if used > bits - 2: break
                if 0 < sel[i] < p + 3:
                    try_raise(i)
    # 6e. raise any word length below 7 (7 passes)
    for _ in range(7):
        if used > bits - 2: break
        for i in range(nused):
            if used > bits - 2: break
            if 0 < sel[i] < 7:
                try_raise(i)

    # ---- 7. scale-factor correction: if the quantised energy exceeds 1.25x the original,
    #         use the next smaller scale factor
    for i in range(nused):
        st = f32(1.0 / (f64(MAXQ[sel[i]]) + 0.5))
        e_orig = e_quant = 0.0
        for k in range(SB_START[i], SB_START[i + 1]):
            x = f64(A[k]); e_orig = e_orig + x * x
            y = f64(Q[k]) * f64(st); e_quant = e_quant + y * y
        if e_orig * 1.25 < e_quant and sfi[i] > 0:
            sfi[i] -= 1
            STATS['sf_decrements'] += 1

    STATS['units'] += 1
    STATS['units_nsub_reduced'] += nused < nsub
    STATS['tone_groups'] += len(u.groups)
    STATS['attack_state_%d' % attack] += 1
    u.nbands, u.nsub, u.mode, u.tone_mode = nbands, nused, 0, TONE_MODE
    u.sel, u.sfi, u.q = sel, sfi, Q
    u.bits = used
    if packed_bits(u) != used:
        raise AssertionError('internal bit count mismatch')
    return u
