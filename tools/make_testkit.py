"""Generate the hardware-encoder test kit (testkit/): sync-marked 16-bit 44.1 kHz stereo WAVs plus a
manifest with exact sample positions, for recording on MiniDisc decks and analysing the ripped ATRAC data.

Every track:  0.5 s silence | sync A (up-chirp) | 0.5 s silence | probe | 0.5 s silence |
              sync B (down-chirp) | 1.0 s silence, padded to a whole number of CD sectors (588 samples).
Deterministic: re-running produces identical files.   Usage: python tools/make_testkit.py [outdir]
"""
import os, sys, json, wave
import numpy as np

SR = 44100
SECTOR = 588
rng = np.random.default_rng(20261013)


def db(x): return 10.0 ** (x / 20.0)
def silence(sec): return np.zeros((int(round(sec * SR)), 2))
def t_axis(sec): return np.arange(int(round(sec * SR))) / SR


def fade(x, ms=5.0):
    n = min(int(SR * ms / 1000), len(x) // 2)
    if n:
        w = 0.5 - 0.5 * np.cos(np.pi * np.arange(n) / n)
        x[:n] *= w[:, None] if x.ndim == 2 else w
        x[-n:] *= w[::-1, None] if x.ndim == 2 else w[::-1]
    return x


def st(l, r=None): return np.stack([l, l if r is None else r], 1)


def tone(f, sec, peak_db, phase=0.0): return db(peak_db) * np.sin(2 * np.pi * f * t_axis(sec) + phase)


def shaped_noise(sec, rms_db, lo=0.0, hi=SR / 2, slope=0.0, notch=None):
    """Gaussian noise with spectrum limited to [lo, hi] Hz, power slope in dB/octave, optional notch (f, octaves)."""
    n = int(round(sec * SR))
    spec = np.fft.rfft(rng.standard_normal(n))
    f = np.fft.rfftfreq(n, 1 / SR)
    g = ((f >= lo) & (f <= hi)).astype(float)
    if slope:
        g *= np.where(f > 0, (np.maximum(f, 20.0) / 1000.0) ** (slope / 6.0206), 0.0)
    if notch:
        fc, octs = notch
        g *= ~((f > fc * 2 ** (-octs / 2)) & (f < fc * 2 ** (octs / 2)))
    x = np.fft.irfft(spec * g, n)
    return x / np.sqrt(np.mean(x ** 2)) * db(rms_db)


def chirp(up=True, sec=0.4, f0=300.0, f1=6000.0, peak_db=-12.0):
    t = t_axis(sec)
    ph = 2 * np.pi * (f0 * t + 0.5 * (f1 - f0) / sec * t ** 2)
    x = np.sin(ph) * db(peak_db)
    n = int(0.01 * SR)
    w = np.ones_like(x); w[:n] = 0.5 - 0.5 * np.cos(np.pi * np.arange(n) / n); w[-n:] = w[:n][::-1]
    x *= w
    return x if up else x[::-1]


SYNC_A, SYNC_B = chirp(True), chirp(False)


class Probe:
    """Collects named segments of a probe signal and records their positions."""
    def __init__(self): self.parts, self.segs, self.pos = [], [], 0
    def add(self, name, x, desc=''):
        x = np.asarray(x, float)
        if x.ndim == 1: x = st(x)
        self.parts.append(x); self.segs.append(dict(name=name, start=self.pos, length=len(x), desc=desc))
        self.pos += len(x)
    def gap(self, sec=0.25): self.add('gap', silence(sec))
    def signal(self): return np.concatenate(self.parts) if self.parts else np.zeros((0, 2))


# ----------------------------------------------------------------------------------- probes
def p_calibration(p):
    p.add('tone1k_-20', fade(st(tone(1000, 6, -20))), '1 kHz, -20 dBFS peak, both channels: path gain')
    p.gap(); p.add('left_only_1k', fade(st(tone(1000, 3, -20), np.zeros(3 * SR))), 'L only: channel map')
    p.gap(); p.add('right_only_2k', fade(st(np.zeros(3 * SR), tone(2000, 3, -20))), 'R only: channel map')
    p.gap(); p.add('white_-20rms', fade(st(shaped_noise(6, -20))), 'white noise, -20 dBFS RMS, identical channels')


def p_level_steps(p):
    for k in range(31):
        lv = -1 - 3 * k
        p.add('tone1k_%ddB' % lv, fade(st(tone(1000, 1.0, lv))), '1 kHz tone at %d dBFS peak' % lv)


def p_noise_steps(p):
    for lv in range(-16, -82, -6):
        p.add('pink_%ddB' % lv, fade(st(shaped_noise(1.5, lv, slope=-3))), 'pink noise %d dBFS RMS' % lv)


def p_transients(p):
    bed = lambda s: tone(440, s, -40)
    x = bed(8.0)
    for k in range(16):                                  # decaying noise bursts (castanet-like)
        s = int((0.25 + 0.5 * k) * SR) + int(rng.integers(0, 1024))
        n = np.arange(3000)
        x[s:s + 3000] += db(-3 - 2 * (k % 8)) * np.exp(-n / 200.0) * rng.standard_normal(3000) / 4.5
    p.add('bursts_on_tone', fade(st(x)), 'decaying noise bursts at -3..-17 dB over a -40 dBFS 440 Hz bed')
    p.gap()
    x = np.zeros(int(6 * SR))
    for k in range(24):                                  # single impulses of different heights
        s = int((0.2 + 0.24 * k) * SR) + int(rng.integers(0, 1024))
        x[s] = db(-1 - 3 * (k % 12)) * (1 if k % 2 else -1)
    p.add('impulses', st(x), 'single-sample impulses, -1..-34 dBFS, random sub-frame positions')
    p.gap()
    x = np.zeros(int(8 * SR))
    for k, f in enumerate([100, 300, 1000, 3000, 6000, 9000, 12000, 15000] * 2):
        s = int((0.2 + 0.5 * k) * SR) + int(rng.integers(0, 1024))
        n = int(0.25 * SR)
        x[s:s + n] = db(-6) * np.sin(2 * np.pi * f * np.arange(n) / SR) * np.exp(-np.arange(n) / (0.08 * SR))
    p.add('tone_onsets', st(x), 'hard-onset decaying tone bursts 100 Hz-15 kHz (pre-echo / gain control)')


def p_attacks(p):
    x = np.zeros(int(16 * SR))
    for k in range(40):                                  # silence -> loud noise, varied band and offset
        s = int((0.15 + 0.38 * k) * SR) + int(rng.integers(0, 1024))
        lo, hi = [(0, 22050), (0, 2000), (2000, 6000), (6000, 12000), (12000, 20000)][k % 5]
        n = int(0.2 * SR)
        x[s:s + n] = shaped_noise(0.2, -16, lo, hi)[:n] * np.exp(-np.arange(n) / (0.06 * SR))
    p.add('attacks', st(x), 'silence to loud band noise, 40 onsets at random positions within a frame')


def p_tones(p):
    fs = [100, 440, 1000, 3000, 5000, 5450, 5512.5, 5575, 7000, 10000, 10950, 11025, 11100, 13000,
          15000, 16000, 16450, 16537.5, 16625, 18000, 20000]
    for f in fs:
        p.add('tone_%g' % f, fade(st(tone(f, 1.5, -12))), 'pure tone %g Hz at -12 dBFS (QMF edges 5512.5/11025/16537.5)' % f)


def p_sweep(p):
    sec = 20.0
    t = t_axis(sec); f0, f1 = 20.0, 21000.0
    k = np.log(f1 / f0) / sec
    x = db(-20) * np.sin(2 * np.pi * f0 * (np.exp(k * t) - 1) / k)
    p.add('log_sweep', fade(st(x)), 'log sweep 20 Hz-21 kHz, -20 dBFS (bandwidth per mode)')


def p_harmonic(p):
    t = t_axis(6.0)
    x = sum(np.sin(2 * np.pi * 110 * h * t) / h for h in range(1, 40) if 110 * h < 20000) * db(-10) / 1.8
    p.add('saw110', fade(st(x)), 'harmonic complex 110 Hz, 1/n amplitudes')
    p.gap()
    ph = 2 * np.pi * np.cumsum(220 * (1 + 0.01 * np.sin(2 * np.pi * 5 * t))) / SR
    x = sum(np.sin(h * ph) / h ** 1.5 for h in range(1, 30)) * db(-10) / 1.3
    p.add('vibrato220', fade(st(x)), 'harmonic complex 220 Hz with 5 Hz vibrato')


def p_masking(p):
    for mf in (1000.0, 4000.0):
        for df in (-0.5, -0.2, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0):
            pf = mf * 2 ** df
            for rel in (-30, -45, -60):
                x = tone(mf, 1.5, -10) + tone(pf, 1.5, -10 + rel, 0.3)
                p.add('mask_%d_%d_%d' % (mf, pf, rel), fade(st(x)),
                      'masker %g Hz -10 dBFS + probe %.1f Hz at %d dB below' % (mf, pf, rel))
    for rel in (-20, -30, -40):
        x = shaped_noise(1.5, -16, 900, 1110) + tone(1000 * 2 ** 0.3, 1.5, -16 + rel)
        p.add('noisemask_%d' % rel, fade(st(x)), 'narrow-band noise masker at 1 kHz + probe 1231 Hz, %d dB' % rel)


def p_noise_spectra(p):
    p.add('white', fade(st(shaped_noise(3, -20))), 'white noise -20 dBFS RMS')
    p.add('pink', fade(st(shaped_noise(3, -20, slope=-3))), 'pink noise')
    p.add('brown', fade(st(shaped_noise(3, -20, slope=-6))), 'brown noise')
    for lo, hi in [(20, 2000), (2000, 5512), (5512, 11025), (11025, 16537), (16537, 22000)]:
        p.add('band_%d_%d' % (lo, hi), fade(st(shaped_noise(2, -20, lo, hi))), 'band noise %d-%d Hz' % (lo, hi))
    for fc in (1000, 4000, 10000):
        p.add('notch_%d' % fc, fade(st(shaped_noise(3, -20, notch=(fc, 1 / 3)))), 'white noise with 1/3-octave notch at %d Hz' % fc)


def p_tone_in_noise(p):
    for snr in (20, 10, 0, -10, -20):
        x = shaped_noise(2, -30) + tone(1000, 2, -30 + snr + 3)
        p.add('tone_snr_%+d' % snr, fade(st(x)), '1 kHz tone at %+d dB SNR in -30 dBFS white noise' % snr)


def p_stereo(p):
    n = lambda s: shaped_noise(s, -20)
    z = lambda s: np.zeros(int(round(s * SR)))
    p.add('noise_L', fade(st(n(3), z(3))), 'noise left only')
    p.add('noise_R', fade(st(z(3), n(3))), 'noise right only')
    a = n(3); p.add('noise_centre', fade(st(a, a)), 'identical noise both channels')
    a = n(3); p.add('noise_antiphase', fade(st(a, -a)), 'noise, right inverted')
    a, b = n(3), n(3); p.add('noise_wide', fade(st(a, b)), 'independent noise per channel')
    sec = 8.0; t = t_axis(sec); a, b = n(sec), n(sec)
    rho = np.cos(np.pi * t / sec)                        # correlation +1 -> -1
    p.add('corr_sweep', fade(st(a, rho * a + np.sqrt(1 - rho ** 2) * b)), 'inter-channel correlation sweep +1 -> -1')
    th = 0.5 * np.pi * t / sec; a = n(sec)
    p.add('pan_noise', fade(st(np.cos(th) * a, np.sin(th) * a)), 'noise panned L -> R')
    s = tone(1000, sec, -14)
    p.add('pan_tone', fade(st(np.cos(th) * s, np.sin(th) * s)), '1 kHz tone panned L -> R')
    p.add('two_tones', fade(st(tone(1000, 3, -14), tone(1500, 3, -14))), 'L 1 kHz, R 1.5 kHz')


def p_quiet(p):
    p.add('digital_zero', silence(2), 'digital silence')
    lsb = 1 / 32768.0
    x = (rng.integers(0, 2, (2 * SR, 2)) + rng.integers(0, 2, (2 * SR, 2)) - 1) * lsb
    p.add('tpdf_1lsb', x, 'TPDF dither noise +-1 LSB')
    p.add('noise_-80', fade(st(shaped_noise(2, -80))), 'white noise -80 dBFS RMS')
    p.add('tone_-90', fade(st(tone(1000, 2, -90))), '1 kHz at -90 dBFS')
    x = np.concatenate([np.zeros(SR), shaped_noise(1, -20), np.zeros(SR)])
    p.add('silence_burst_silence', st(x), 'silence, 1 s noise, silence (hard edges)')


def p_pressure(p):
    p.add('white_-16', fade(st(shaped_noise(4, -16))), 'full-band white noise, -16 dBFS RMS')
    t = t_axis(4)
    fs = np.exp(rng.uniform(np.log(60), np.log(18000), 60))
    x = sum(np.sin(2 * np.pi * f * t + rng.uniform(0, 6.3)) for f in fs)
    x *= db(-3) / np.max(np.abs(x))
    p.add('multitone60', fade(st(x)), '60 random sine tones 60 Hz-18 kHz, dense, -3 dBFS peak')
    p.add('hf_noise', fade(st(shaped_noise(4, -16, slope=+3))), 'treble-heavy noise (+3 dB/oct), -16 dBFS RMS')


def p_music(p):
    """A short synthetic 'song': bass, chords with decaying piano-like partials, hi-hat and snare hits."""
    sec, bpm = 24.0, 112
    n = int(sec * SR); x = np.zeros((n, 2)); beat = 60.0 / bpm
    def note(f, start, dur, amp, pan):
        s, m = int(start * SR), int(dur * SR)
        if s >= n: return
        m = min(m, n - s); tt = np.arange(m) / SR
        y = sum(np.sin(2 * np.pi * f * h * tt * (1 + 0.0004 * h * h)) * np.exp(-tt * (1.5 + 0.8 * h)) / h for h in range(1, 12))
        y *= amp * (1 - np.exp(-tt / 0.003))
        x[s:s + m, 0] += y * np.cos(pan); x[s:s + m, 1] += y * np.sin(pan)
    chords = [(220.0, 277.18, 329.63), (196.0, 246.94, 293.66), (174.61, 220.0, 261.63), (196.0, 246.94, 329.63)]
    for bar in range(int(sec / (4 * beat))):
        ch = chords[bar % 4]
        for b in range(4):
            tb = (bar * 4 + b) * beat
            note(ch[0] / 2, tb, beat, 0.22, 0.785)
            for k, f in enumerate(ch): note(f, tb + 0.01 * k, 2 * beat, 0.09, 0.5 + 0.3 * k)
            for h in range(2):
                s = int((tb + h * beat / 2) * SR); m = min(int(0.05 * SR), n - s)
                if m > 0: x[s:s + m] += (shaped_noise(0.05, -18, 7000, 20000)[:m] * np.exp(-np.arange(m) / 300.0))[:, None] * [0.6, 1.0]
            if b % 2 == 1:
                s = int(tb * SR); m = min(int(0.2 * SR), n - s)
                if m > 0: x[s:s + m] += (shaped_noise(0.2, -10, 150, 9000)[:m] * np.exp(-np.arange(m) / 1800.0))[:, None]
    x *= db(-1.5) / np.max(np.abs(x))
    p.add('synth_song', fade(x), 'synthetic pop-like excerpt (bass, piano-like chords, hats, snare)')


TRACKS = [
    ('01_calibration', p_calibration), ('02_level_steps', p_level_steps), ('03_noise_steps', p_noise_steps),
    ('04_transients', p_transients), ('05_attacks', p_attacks), ('06_tones', p_tones), ('07_sweep', p_sweep),
    ('08_harmonic', p_harmonic), ('09_masking', p_masking), ('10_noise_spectra', p_noise_spectra),
    ('11_tone_in_noise', p_tone_in_noise), ('12_stereo', p_stereo), ('13_quiet', p_quiet),
    ('14_pressure', p_pressure), ('15_music', p_music),
]


def build(outdir):
    os.makedirs(outdir, exist_ok=True)
    manifest = dict(sample_rate=SR, sync=dict(f0=300.0, f1=6000.0, seconds=0.4, peak_dbfs=-12.0,
                                              a='up-chirp', b='down-chirp'), tracks=[])
    total = 0
    for name, fn in TRACKS:
        p = Probe(); fn(p)
        probe = p.signal()
        peak = np.max(np.abs(probe))
        assert peak <= db(-0.5), (name, 20 * np.log10(peak))
        parts = [silence(0.5), st(SYNC_A), silence(0.5)]
        a_start = int(round(0.5 * SR)); probe_start = a_start + len(SYNC_A) + int(round(0.5 * SR))
        parts += [probe, silence(0.5)]
        b_start = probe_start + len(probe) + int(round(0.5 * SR))
        parts += [st(SYNC_B), silence(1.0)]
        x = np.concatenate(parts)
        pad = (-len(x)) % SECTOR
        x = np.concatenate([x, np.zeros((pad, 2))])
        pcm = np.round(np.clip(x, -1, 32767 / 32768) * 32768).astype('<i2')
        path = os.path.join(outdir, name + '.wav')
        with wave.open(path, 'wb') as w:
            w.setnchannels(2); w.setsampwidth(2); w.setframerate(SR); w.writeframes(pcm.tobytes())
        manifest['tracks'].append(dict(file=name + '.wav', samples=len(x), sync_a=a_start, sync_b=b_start,
                                       probe_start=probe_start, probe_samples=len(probe),
                                       peak_dbfs=round(20 * np.log10(peak), 2), segments=p.segs))
        total += len(x)
        print('%-18s %6.1f s  peak %6.2f dBFS  %3d segments' % (name, len(x) / SR, 20 * np.log10(peak), len(p.segs)))
    with open(os.path.join(outdir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1)
    print('total %.1f s (%.1f min)' % (total / SR, total / SR / 60))


if __name__ == '__main__':
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'testkit'))
