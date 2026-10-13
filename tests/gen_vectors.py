"""Generate the regression-test WAV files (16-bit, 44.1 kHz).

Standard library only, and bit-identical on every platform: the signals use nothing but
IEEE basic arithmetic (no libm, no platform RNG), so CI on Linux, macOS, Windows and
WebAssembly encodes exactly the same input.

Usage:  python tests/gen_vectors.py OUTDIR
"""
import os, struct, sys

SR = 44100
PI = 3.141592653589793


class Rng:
    """xorshift32: uniform floats in [-1, 1)."""

    def __init__(self, seed):
        self.s = seed

    def u(self):
        s = self.s
        s ^= (s << 13) & 0xffffffff
        s ^= s >> 17
        s ^= (s << 5) & 0xffffffff
        self.s = s
        return s / 2147483648.0 - 1.0

    def n(self):
        """Roughly Gaussian (Irwin-Hall, 4 terms), unit variance."""
        return (self.u() + self.u() + self.u() + self.u()) * 0.8660254037844386


def sin_turns(t):
    """sin(2*pi*t) from a Taylor polynomial (accurate to ~1e-9)."""
    t -= int(t) if t >= 0 else int(t) - 1           # [0, 1)
    if t >= 0.5:
        return -sin_turns(t - 0.5)
    if t > 0.25:
        t = 0.5 - t
    x = 2.0 * PI * t                                # [0, pi/2]
    x2 = x * x
    return x * (1 - x2 / 6 * (1 - x2 / 20 * (1 - x2 / 42 * (1 - x2 / 72 * (1 - x2 / 110)))))


class Osc:
    """Phase-accumulating sine oscillator with a per-sample frequency."""

    def __init__(self, phase=0.0):
        self.ph = phase

    def __call__(self, hz):
        v = sin_turns(self.ph)
        self.ph += hz / SR
        if self.ph >= 1.0:
            self.ph -= 1.0
        return v


def add(values):
    """Left-to-right float sum (built-in sum() is compensated since Python 3.12)."""
    s = 0.0
    for v in values:
        s = s + v
    return s


def save(path, chans):
    n = len(chans[0])
    out = bytearray()
    for i in range(n):
        for c in chans:
            v = int(round(c[i] * 32767.0))
            out += struct.pack('<h', max(-32768, min(32767, v)))
    hdr = b'RIFF' + struct.pack('<I', 36 + len(out)) + b'WAVEfmt ' + struct.pack(
        '<IHHIIHH', 16, 1, len(chans), SR, SR * 2 * len(chans), 2 * len(chans), 16) + \
        b'data' + struct.pack('<I', len(out))
    with open(path, 'wb') as f:
        f.write(hdr + out)


def chirp_noise(n):
    r, a, b, am = Rng(1), Osc(), Osc(0.3), Osc()
    L, R = [], []
    for i in range(n):
        t = i / SR
        l = 0.3 * a(200 + 3000 * t) * (1 + 0.5 * am(3))
        if (i % 22050) < 2205:
            l += 0.2 * r.n()
        L.append(l)
        R.append(0.25 * b(5000 - 1500 * t) + 0.05 * r.n())
    return [L, R]


def transients(n):
    r, o = Rng(2), Osc()
    x = []
    env = 0.0
    for i in range(n):
        if i % 11025 == 4410:
            env = 0.8
        x.append(0.05 * o(440) + env * r.n())
        env *= 0.9933
    return [x, [0.7 * v for v in x]]


def harmonics(n):
    hs = [Osc(0.1 * h) for h in range(12)]
    fs = [Osc(0.05 * k) for k in range(4)]
    L, R = [], []
    for i in range(n):
        L.append(add(0.12 / (h + 1) * hs[h](523.25 * (h + 1)) for h in range(12)))
        R.append(add(0.1 * fs[k](f) for k, f in enumerate((1000, 2500, 7000, 12000))))
    return [L, R]


def single_tone(n):
    a, b = Osc(), Osc()
    return [[0.5 * a(3000) for _ in range(n)], [0.4 * b(9000) for _ in range(n)]]


def correlated(n):
    r, o = Rng(3), Osc()
    L, R = [], []
    for i in range(n):
        m = 0.3 * o(330) + 0.1 * r.n()
        L.append(m)
        R.append(0.9 * m + 0.02 * r.n())
    return [L, R]


def noise(n):
    r = Rng(4)
    return [[0.5 * r.n() for _ in range(n)], [0.5 * r.n() for _ in range(n)]]


def silence_burst(n):
    r = Rng(5)
    return [[0.0 if i < n // 2 else 0.3 * r.n() for i in range(n)] for _ in range(2)]


def onesided(n):
    os_ = [Osc() for _ in range(29)]
    L, R = [], []
    for i in range(n):
        s = 0.15 * add(os_[k](220 * (k + 1)) / (k + 1) for k in range(29))
        L.append(s if i < n // 2 else 0.0)
        R.append(0.0 if i < n // 2 else s)
    return [L, R]


def comb(n):
    os_ = [Osc(0.01 * k * k % 1.0) for k in range(1, 60)]
    L = []
    for i in range(n):
        f0 = 95 + 10 * i / SR
        L.append(0.02 * add(o(f0 * (k + 1)) for k, o in enumerate(os_)))
    return [L, [0.8 * L[(i - 37) % n] for i in range(n)]]


def quiet(n):
    a, b, c, d = Osc(), Osc(), Osc(), Osc()
    L, R = [], []
    for i in range(n):
        q = 0.001 * a(1234.5) + 0.0003 * b(7777) + 0.0001 * c(15000)
        L.append(q)
        R.append(0.5 * q + 0.0002 * d(300))
    return [L, R]


def square(n):
    o = Osc()
    L = [0.99 if o(60 + 2000 * i / SR) >= 0 else -0.99 for i in range(n)]
    return [L, [-L[(i - 11) % n] for i in range(n)]]


def clicks(n):
    o = Osc()
    c = [0.01 * o(80) for _ in range(n)]
    for k in range(n // 5160):
        p = 2205 + 5160 * k
        amp = 0.95 if k % 2 else -0.7
        for j in range(30):
            if p + j < n:
                c[p + j] += amp
            amp *= 0.8465
    return [c, [c[(i - 2000) % n] for i in range(n)]]


def panned(n):
    r, lfo = Rng(6), Osc()
    L, R = [], []
    for i in range(n):
        g = 0.5 + 0.5 * lfo(0.7)
        v = 0.3 * r.n()
        L.append(g * v)
        R.append((1 - g) * v + 0.05 * r.n())
    return [L, R]


def music(n, seed=7):
    """A short synthetic arrangement: chords with envelopes, bass, kick, hi-hat."""
    r = Rng(seed)
    notes = [261.63, 329.63, 392.0, 440.0, 349.23, 293.66]
    voices = [Osc() for _ in range(3)]
    bass, kick = Osc(), Osc()
    L, R = [], []
    beat = SR // 4
    env = kenv = henv = 0.0
    kf = 0.0
    for i in range(n):
        if i % beat == 0:
            env = 1.0
            chord = (i // beat) % len(notes)
            if (i // beat) % 2 == 0:
                kenv, kf = 1.0, 160.0
        if i % (beat // 2) == beat // 4:
            henv = 1.0
        tone = add(voices[k](notes[(chord + 2 * k) % len(notes)] * (1 + k)) for k in range(3))
        b = bass(notes[chord] / 4)
        k = kick(kf)
        h = r.n()
        mono = 0.12 * env * tone + 0.25 * b + 0.6 * kenv * k
        L.append(mono + 0.1 * henv * h)
        R.append(0.8 * mono + 0.07 * henv * h + 0.05 * env * tone)
        env *= 0.99993
        kenv *= 0.9995
        kf = 40 + (kf - 40) * 0.9993
        henv *= 0.997
    return [L, R]


def mono(sig):
    return lambda n: sig(n)[:1]


VECTORS = [                        # name, generator, samples
    ('chirp_noise', chirp_noise, 88200),
    ('transients', transients, 88200),
    ('harmonics', harmonics, 66150),
    ('single_tone', single_tone, 66150),
    ('correlated', correlated, 88200),
    ('noise', noise, 66150),
    ('silence_burst', silence_burst, 66150),
    ('onesided', onesided, 88200),
    ('comb', comb, 66150),
    ('quiet', quiet, 66150),
    ('square', square, 66150),
    ('clicks', clicks, 88200),
    ('panned', panned, 88200),
    ('music', music, 132300),
    ('music_mono', mono(music), 132300),
    ('transients_mono', mono(transients), 66150),
    ('harmonics_mono', mono(harmonics), 44100),
    ('len_6145', chirp_noise, 6145),          # just above at3tool's minimum input
    ('len_9285', music, 9285),                # 955 + N is a multiple of 1024
    ('len_10171', transients, 10171),         # (955 + N) mod 1024 == 886: one extra flush frame
    ('len_10171_mono', mono(transients), 10171),
    ('len_100', music, 100),                  # shorter than one frame
]


def main_into(out):
    os.makedirs(out, exist_ok=True)
    for name, gen, n in VECTORS:
        save(os.path.join(out, name + '.wav'), gen(n))


if __name__ == '__main__':
    main_into(sys.argv[1] if len(sys.argv) > 1 else 'vectors')
