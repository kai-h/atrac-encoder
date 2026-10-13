"""LP2 encoder: sound-unit data structures, bit counting and serialisation.

Bit counting follows Sony's functions around 0x43ce00; serialisation follows 0x43f840.  The syntax itself is the one parsed by atrac3_decoder.py.
"""
import numpy as np
from lp2_tables import QBITS, VALS_PER_CODE, ENC_CODES, CODE_LEN, SB_START, TONE_NVALS


def mantissa_bits(mode, sel, q):
    """Bits for the mantissas q with code table (mode, sel) (0x43fdd0)."""
    b = QBITS[sel]; m = (1 << b) - 1; L = CODE_LEN[(mode, sel)]
    q = np.asarray(q, np.int64) & m
    if VALS_PER_CODE[sel] == 1:
        return int(L[q].sum())
    return int(L[(q[0::2] << b) | q[1::2]].sum())


class Component:
    """One tonal component: 4 consecutive MDCT lines starting at `pos` (multiple of 4)."""
    __slots__ = ('pos', 'sfi', 'mant', 'len_code', 'sel', 'mode')


class ToneGroup:
    """One tonal-component group: common length code and quantiser, components per 64-line cell."""

    def __init__(self, sel, len_code, mode):
        self.sel, self.len_code, self.mode = sel, len_code, mode
        self.flags = [0, 0, 0, 0]               # QMF bands that carry components
        self.cells = [[] for _ in range(16)]    # components per 64-line cell (max 7 each)


class SoundUnit:
    """Everything that is serialised for one channel in one frame."""

    def __init__(self):
        self.nbands = 0          # QMF bands coded (1..4)
        self.gain = []           # GainInfo per band
        self.groups = []         # ToneGroup list (max 31)
        self.tone_mode = 1       # 2-bit tonal coding-mode selector (always 1 = CLC here)
        self.nsub = 0            # subbands coded (1..32)
        self.mode = 0            # spectrum coding mode (0 = VLC; Sony LP2 always uses VLC)
        self.sel = []            # word length per subband
        self.sfi = []            # scale-factor index per subband
        self.q = None            # 1024 mantissas (only coded subbands are meaningful)


def component_bits(u):
    n = 7 if u.groups else 5
    for G in u.groups:
        n += u.nbands + 6
        for c in range(4 * u.nbands):
            if G.flags[c >> 2]:
                n += 3
                for cp in G.cells[c]:
                    cnt = min(TONE_NVALS[cp.len_code], 1024 - cp.pos)
                    n += 12 + mantissa_bits(cp.mode, cp.sel, cp.mant[:cnt])
    return n


def spectrum_bits(u):
    n = 3 * u.nsub + 6
    for i in range(u.nsub):
        if u.sel[i] > 0:
            n += 6 + mantissa_bits(u.mode, u.sel[i], u.q[SB_START[i]:SB_START[i + 1]])
    return n


def gain_bits(gain, nbands):
    """Gain-control data: 3 bits per band + 9 bits per point."""
    return 3 * (nbands + 3 * sum(len(gain[b].pts) for b in range(nbands)))


def packed_bits(u):
    return component_bits(u) + spectrum_bits(u) + gain_bits(u.gain, u.nbands) + 8


class BitWriter:
    """MSB-first writer into a zero-filled buffer (0x440190)."""

    def __init__(self, nbytes):
        self.buf = bytearray(nbytes); self.pos = 0

    def put(self, v, n):
        for k in range(n - 1, -1, -1):
            if (v >> k) & 1:
                self.buf[self.pos >> 3] |= 0x80 >> (self.pos & 7)
            self.pos += 1


def pack_mantissas(bw, mode, sel, q):
    b = QBITS[sel]; m = (1 << b) - 1; tab = ENC_CODES[(mode, sel)]
    if VALS_PER_CODE[sel] == 1:
        for v in q:
            c, l = tab[int(v) & m]; bw.put(c, l)
    else:
        for i in range(0, len(q), 2):
            c, l = tab[((int(q[i]) & m) << b) | (int(q[i + 1]) & m)]; bw.put(c, l)


def pack_unit(u, nbytes):
    """Serialise a SoundUnit into `nbytes` bytes."""
    bw = BitWriter(nbytes)
    bw.put(0x28, 6)
    bw.put(u.nbands - 1, 2)
    for b in range(u.nbands):
        bw.put(len(u.gain[b].pts), 3)
        for lev, loc in u.gain[b].pts:
            bw.put(lev, 4); bw.put(loc, 5)
    bw.put(len(u.groups), 5)
    if u.groups:
        bw.put(u.tone_mode, 2)
    for G in u.groups:
        for b in range(u.nbands):
            bw.put(G.flags[b], 1)
        bw.put(G.len_code, 3); bw.put(G.sel, 3)
        for c in range(4 * u.nbands):
            if G.flags[c >> 2]:
                bw.put(len(G.cells[c]), 3)
                for cp in G.cells[c]:
                    bw.put(cp.sfi, 6); bw.put(cp.pos & 63, 6)
                    n = min(TONE_NVALS[G.len_code], 0x400 - cp.pos)
                    pack_mantissas(bw, G.mode, G.sel, cp.mant[:n])
    bw.put(u.nsub - 1, 5)
    bw.put(u.mode, 1)
    for i in range(u.nsub):
        bw.put(u.sel[i], 3)
    for i in range(u.nsub):
        if u.sel[i] > 0:
            bw.put(u.sfi[i], 6)
    for i in range(u.nsub):
        if u.sel[i] > 0:
            pack_mantissas(bw, u.mode, u.sel[i], u.q[SB_START[i]:SB_START[i + 1]])
    if bw.pos > 8 * nbytes:
        raise ValueError('sound unit exceeds %d bytes' % nbytes)
    return bytes(bw.buf)
