"""ATRAC3 reference decoder (executable companion to the spec).

Straightforward float64 implementation of the bitstream syntax and synthesis described in
spec/ATRAC3.md. Not optimised. Usage:  python atrac3_decoder.py in.at3 out.wav
"""
import sys, struct, wave, collections
import numpy as np
from at3_tables import CODES, MAX_QUANT, SUBBAND_START

STATS = collections.Counter()   # feature usage counters (informational only)
SF = [2.0 ** ((i - 15) / 3.0) for i in range(64)]
# joint stereo: matrix selector -> (coef for channel-0 signal, coef for channel-1 signal); 2 is invalid
MATRIX = {0: (0.0, 1.0), 1: (1.0, 1.0), 3: (0.5, 0.5)}
# weight index (flag*8 + level) -> (w0, w1); level 7 = no weighting (2, 2)
WEIGHTS = {}
for _lv in range(8):
    _a = 2.0 * _lv / 7.0
    _b = 2.0 * (2.0 - (_lv / 7.0) ** 2) ** 0.5
    WEIGHTS[_lv], WEIGHTS[8 + _lv] = (_a, _b), (_b, _a)


class BitReader:
    def __init__(self, data):
        self.data, self.pos = data, 0

    def read(self, n):
        v = 0
        for _ in range(n):
            byte = self.data[self.pos >> 3] if (self.pos >> 3) < len(self.data) else 0
            v = (v << 1) | ((byte >> (7 - (self.pos & 7))) & 1)
            self.pos += 1
        return v

    def read_code(self, table):
        code = ''
        while code not in table:
            code += '1' if self.read(1) else '0'
            if len(code) > 8: raise ValueError('bad codeword')
        return table[code]


def read_mantissas(br, selector, mode, count):
    """Return `count` integer mantissas for a quantiser selector (1..7)."""
    table = CODES[(mode, selector)]
    out = []
    if selector == 1:
        while len(out) < count:
            out.extend(br.read_code(table))
        return out[:count]
    for _ in range(count):
        out.append(br.read_code(table))
    return out


class SoundUnit:
    """Per-channel decoder state."""
    def __init__(self):
        self.prev_gain = [[] for _ in range(4)]      # gain points of previous frame per QMF band
        self.overlap = np.zeros((4, 256))            # second half of previous IMDCT per band
        self.qmf_delay = [np.zeros(46), np.zeros(46), np.zeros(46)]


def decode_unit(br, js_second):
    """Parse one sound unit. Returns (bands_coded, gain, spectrum[1024])."""
    if js_second:
        if br.read(2) != 3: raise ValueError('JS sound unit id')
    elif br.read(6) != 0x28:
        raise ValueError('sound unit id')
    nbands = br.read(2) + 1
    gain = []
    for b in range(nbands):
        n = br.read(3)
        pts = [(br.read(4), br.read(5)) for _ in range(n)]
        for k in range(1, n):
            if pts[k][1] <= pts[k - 1][1]: raise ValueError('gain locations not increasing')
        gain.append(pts)
    gain += [[] for _ in range(4 - nbands)]
    spec = np.zeros(1024)
    # tonal components
    ngroups = br.read(5)
    if ngroups:
        sel_mode = br.read(2)
        if sel_mode == 2: raise ValueError('tonal coding mode 2')
        ncomp = 0
        for _ in range(ngroups):
            flags = [br.read(1) for _ in range(nbands)]
            nvalues = br.read(3) + 1
            qstep = br.read(3)
            if qstep <= 1: raise ValueError('tonal quant step')
            mode = br.read(1) if sel_mode == 3 else sel_mode
            for b in range(nbands):
                if not flags[b]: continue
                for g in range(4):
                    for _ in range(br.read(3)):
                        sfi = br.read(6)
                        pos = b * 256 + g * 64 + br.read(6)
                        cnt = min(nvalues, 1024 - pos)
                        m = read_mantissas(br, qstep, mode, cnt)
                        scale = SF[sfi] / MAX_QUANT[qstep]
                        for k in range(cnt): spec[pos + k] += m[k] * scale
                        ncomp += 1
        STATS['tonal_components'] += ncomp
        STATS['tonal_mode_sel_%d' % sel_mode] += 1
    # spectrum
    nsub = br.read(5) + 1
    mode = br.read(1)
    STATS['units'] += 1
    STATS['clc_units'] += mode
    STATS['gain_points'] += sum(len(g) for g in gain)
    sels = [br.read(3) for _ in range(nsub)]
    sfs = [br.read(6) if s else 0 for s in sels]
    for i in range(nsub):
        if not sels[i]: continue
        a, b = SUBBAND_START[i], SUBBAND_START[i + 1]
        m = read_mantissas(br, sels[i], mode, b - a)
        spec[a:b] += np.array(m) * (SF[sfs[i]] / MAX_QUANT[sels[i]])
    return nbands, gain, spec


# ---------------------------------------------------------------- synthesis
def _imdct_matrix():
    n = np.arange(512)[:, None]; k = np.arange(256)[None, :]
    return np.cos(np.pi / 256.0 * (n + 0.5 + 128) * (k + 0.5))

IMDCT = _imdct_matrix()

def _window():
    w = np.zeros(512)
    for i in range(256):
        j = 255 - i
        wi = np.sin(((i + 0.5) / 256.0 - 0.5) * np.pi) + 1.0
        wj = np.sin(((j + 0.5) / 256.0 - 0.5) * np.pi) + 1.0
        w[i] = w[511 - i] = wi / (0.5 * (wi * wi + wj * wj))
    return w

WINDOW = _window()


def gain_curve(points):
    """Per-sample gain over 256 samples for one set of gain points (applied to overlap region)."""
    g = np.ones(256)
    if not points: return g
    levels = [p[0] for p in points] + [4]
    pos = 0
    for i, (lev, loc) in enumerate(points):
        start = loc * 8
        cur = 2.0 ** (4 - lev)
        inc = 2.0 ** ((lev - levels[i + 1]) / 8.0)
        g[pos:start] = cur
        for k in range(8):
            if start + k < 256: g[start + k] = cur * inc ** k
        pos = start + 8
    return g


def synth_band(unit, band, coefs, gain_now, gain_next):
    if band & 1: coefs = coefs[::-1]
    y = IMDCT @ coefs * WINDOW
    scale_next = 2.0 ** (4 - gain_next[0][0]) if gain_next else 1.0
    out = (y[:256] * scale_next + unit.overlap[band]) * gain_curve(gain_now)
    unit.overlap[band] = y[256:]
    return out


QMF_HALF = [-0.00001461907, -0.00009205479, -0.000056157569, 0.00030117269, 0.0002422519,
            -0.00085293897, -0.0005205574, 0.0020340169, 0.00078333891, -0.0042153862,
            -0.00075614988, 0.0078402944, -0.000061169922, -0.01344162, 0.0024626821,
            0.021736089, -0.007801671, -0.034090221, 0.01880949, 0.054326009,
            -0.043596379, -0.099384367, 0.13207909, 0.46424159]
QMF_WIN = np.array(QMF_HALF + QMF_HALF[::-1]) * 2.0


def iqmf(lo, hi, delay):
    n = len(lo)
    tmp = np.concatenate([delay, np.empty(2 * n)])
    tmp[46::2] = lo + hi
    tmp[47::2] = lo - hi
    out = np.empty(2 * n)
    for j in range(n):
        p = tmp[2 * j: 2 * j + 48]
        out[2 * j] = np.dot(p[1::2], QMF_WIN[1::2])
        out[2 * j + 1] = np.dot(p[0::2], QMF_WIN[0::2])
    delay[:] = tmp[2 * n: 2 * n + 46]
    return out


class Decoder:
    def __init__(self, channels, frame_bytes, joint_stereo):
        self.ch, self.fb, self.js = channels, frame_bytes, joint_stereo
        self.units = [SoundUnit() for _ in range(channels)]
        # joint-stereo parameters in effect: (weight index 0..15, matrix selectors[4], last coded band)
        self.js_prev = self.js_cur = (None, [3, 3, 3, 3], 4)

    def decode_frame(self, data):
        per = self.fb // self.ch
        parsed = []
        if not self.js:
            for c in range(self.ch):
                parsed.append(decode_unit(BitReader(data[c * per:(c + 1) * per]), False))
        else:
            parsed.append(decode_unit(BitReader(data), False))
            rev = data[::-1]
            k = 0
            while rev[k] == 0xF8: k += 1
            br = BitReader(rev[k:])
            wflag, widx = br.read(1), br.read(3)
            mat = [br.read(2) for _ in range(4)]
            parsed.append(decode_unit(br, True))
            if 2 in mat: raise ValueError('matrix selector 2 is invalid')
            self.js_next = (wflag * 8 + widx, mat, parsed[1][0] - 1)
        bands_out = []
        for c, (nb, gain, spec) in enumerate(parsed):
            u = self.units[c]
            bands = np.zeros((4, 256))
            for b in range(4):
                bands[b] = synth_band(u, b, spec[b * 256:(b + 1) * 256], u.prev_gain[b], gain[b])
            u.prev_gain = gain
            bands_out.append(bands)
        if self.js:
            self.apply_js(bands_out)
        out = []
        for c in range(self.ch):
            b = bands_out[c]
            d = self.units[c].qmf_delay
            low = iqmf(b[0], b[1], d[0])
            high = iqmf(b[3], b[2], d[1])
            out.append(iqmf(low, high, d[2]))
        return np.stack(out, 1)

    def apply_js(self, bands):
        """Inverse matrixing + channel weighting on the QMF-band time signals (in place)."""
        (pw, pm, pnb), (cw, cm, cnb) = self.js_prev, self.js_cur
        L, R = bands
        for b in range(4):
            p0, p1 = MATRIX[pm[b]]
            c0, c1 = MATRIX[cm[b]]
            pw0, pw1 = WEIGHTS[pw] if b > pnb else (2.0, 2.0)
            cw0, cw1 = WEIGHTS[cw] if b > cnb else (2.0, 2.0)
            start = 0
            if pw0 != cw0 or pm[b] != cm[b]:
                for n in range(8):
                    t = (n - 8) / 8.0
                    s = (c0 + t * (c0 - p0)) * L[b, n] + (c1 + t * (c1 - p1)) * R[b, n]
                    R[b, n] = (L[b, n] - s) * (cw1 + t * (cw1 - pw1))
                    L[b, n] = (cw0 + t * (cw0 - pw0)) * s
                start = 8
            s = c0 * L[b, start:] + c1 * R[b, start:]
            R[b, start:] = (L[b, start:] - s) * cw1
            L[b, start:] = s * cw0
        self.js_prev, self.js_cur = self.js_cur, self.js_next


def read_at3(path):
    data = open(path, 'rb').read()
    assert data[:4] == b'RIFF' and data[8:12] == b'WAVE'
    pos, fmt, body = 12, None, None
    while pos + 8 <= len(data):
        cid, size = data[pos:pos + 4], struct.unpack('<I', data[pos + 4:pos + 8])[0]
        if cid == b'fmt ': fmt = data[pos + 8:pos + 8 + size]
        if cid == b'data': body = data[pos + 8:pos + 8 + size]
        pos += 8 + size + (size & 1)
    tag, ch, sr, _, align = struct.unpack('<HHIIH', fmt[:14])
    assert tag == 0x270
    js = struct.unpack('<H', fmt[24:26])[0] if len(fmt) >= 26 else 0
    return ch, align, js, body


if __name__ == '__main__':
    ch, align, js, body = read_at3(sys.argv[1])
    dec = Decoder(ch, align, js)
    frames = [dec.decode_frame(body[i:i + align]) for i in range(0, len(body) - align + 1, align)]
    pcm = np.round(np.clip(np.concatenate(frames), -32767, 32767)).astype('<i2')
    with wave.open(sys.argv[2], 'wb') as w:
        w.setnchannels(ch); w.setsampwidth(2); w.setframerate(44100); w.writeframes(pcm.tobytes())
