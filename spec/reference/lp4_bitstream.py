"""Sound unit packing of the LP4 encoder (0x437490 header/gain/tones, 0x4376b0 tonal block,
0x434610 spectral codes).  Bits are written MSB first; the frame buffer starts zeroed."""
import numpy as np
from lp4_tables import SB_START, QMASK, VLC
from lp4_alloc import quantise


class BitWriter:
    def __init__(self, buf, bitpos):
        self.buf = buf
        self.pos = bitpos

    def put(self, value, n):
        for k in range(n - 1, -1, -1):
            if (value >> k) & 1:
                self.buf[self.pos >> 3] |= 0x80 >> (self.pos & 7)
            self.pos += 1


def pack_unit(buf, byte_off, S, gain, js_header=None):
    """Write one sound unit starting at byte byte_off of buf (a bytearray, zero filled).

    S: lp4_alloc.AllocState; gain: list of 4 lp4_gain.BandGain (points of this frame);
    js_header: (weight_index, [matrix selectors x4]) for the second unit of a joint-stereo
    frame, else None.  Returns the bit position after the unit."""
    if js_header is not None:
        w, mat = js_header
        buf[byte_off] = ((w * 4 + mat[0]) * 4 + mat[1]) & 0xff
        byte_off += 1
        head = (mat[3] + mat[2] * 4) * 16 + 0x0c          # 2 matrix selectors + 2-bit id 3
    else:
        head = 0xa0                                       # 6-bit id 0x28
    buf[byte_off] = (S.nbands + head - 1) & 0xff
    bw = BitWriter(buf, byte_off * 8 + 8)
    for b in range(S.nbands):
        G = gain[b]
        bw.put(G.count, 3)
        for k in range(G.count):
            bw.put(G.lev[k] * 32 + G.loc[k], 9)
    bw.put(S.ngroups, 5)
    if S.ngroups:
        _pack_tones(bw, S)
    bw.put(S.nsub * 2 - 2, 6)                             # (nsub - 1) and coding mode 0 = VLC
    for i in range(S.nsub):
        bw.put(S.sel[i], 3)
    for i in range(S.nsub):
        if S.sel[i]:
            bw.put(S.sf[i], 6)
    for i in range(S.nsub):
        sel = S.sel[i]
        if not sel:
            continue
        q, _ = quantise(S.spec[SB_START[i]:SB_START[i + 1]], sel, S.sf[i])
        tab = VLC[(0, sel)]
        if sel == 1:
            q = q & 3
            for j in range(0, len(q), 2):
                code, ln = tab[int(q[j]) * 4 + int(q[j + 1])]
                bw.put(code, ln)
        else:
            for idx in (q & QMASK[sel]).tolist():
                code, ln = tab[idx]
                bw.put(code, ln)
    return bw.pos


def _pack_tones(bw, S):
    bw.put(S.sel_mode, 2)
    for gi in range(S.ngroups):
        g = S.groups[gi]
        v = g.flag[0]
        for k in range(1, S.nbands):
            v = v * 2 + g.flag[k]
        bw.put(v, S.nbands)
        bw.put(g.qstep + g.nvm1 * 8, 6)
        tab = VLC[(S.sel_mode, g.qstep)]
        for blk in range(S.nbands * 4):
            if not g.flag[blk >> 2]:
                continue
            bw.put(len(g.lst[blk]), 3)
            for rec in g.lst[blk]:
                bw.put(rec.sf, 6)
                bw.put(rec.pos & 0x3f, 6)
                for k in range(g.nvm1 + 1):
                    code, ln = tab[rec.q[k]]
                    bw.put(code, ln)
