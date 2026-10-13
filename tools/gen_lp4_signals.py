"""Extra deterministic test signals aimed at rarely taken paths of the LP4 encoder model
(one-sided stereo, dense harmonics, very quiet tones, clipping squares, clicks, odd lengths)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_signals import save, SR


def signals(dur=3.0):
    rng = np.random.default_rng(4321)
    n = int(SR * dur)
    t = np.arange(n) / SR
    out = {}
    # left only (right digital silence), then right only: one-sided stereo decisions
    saw = sum(np.sin(2 * np.pi * 220 * k * t) / k for k in range(1, 30)) * 0.15
    x = np.zeros((n, 2))
    x[:n // 2, 0] = saw[:n // 2]
    x[n // 2:, 1] = saw[n // 2:]
    out['onesided'] = x
    # dense harmonic comb (many tonal components per 64-coefficient block), slowly gliding
    f0 = 95 + 10 * t
    ph = 2 * np.pi * np.cumsum(f0) / SR
    comb = sum(np.sin(k * ph + k * k) for k in range(1, 60)) * 0.02
    out['comb'] = np.stack([comb, np.roll(comb, 37) * 0.8], 1)
    # very quiet sines (-60..-80 dBFS): lots of spare bits
    q = 0.001 * np.sin(2 * np.pi * 1234.5 * t) + 0.0003 * np.sin(2 * np.pi * 7777 * t) \
        + 0.0001 * np.sin(2 * np.pi * 15000 * t)
    out['quiet'] = np.stack([q, 0.5 * q + 0.0002 * np.sin(2 * np.pi * 300 * t)], 1)
    # clipped square waves with a sweep: more bits wanted than available
    sq = np.sign(np.sin(2 * np.pi * (60 + 2000 * t) * t)) * 0.99
    out['square'] = np.stack([sq, -np.roll(sq, 11)], 1)
    # isolated clicks in silence and over a low tone (gain control in all bands)
    c = 0.01 * np.sin(2 * np.pi * 80 * t)
    for k in range(25):
        p = int((0.05 + 0.117 * k) * SR)
        c[p:p + 30] += (0.95 if k % 2 else -0.7) * np.exp(-np.arange(30) / 6.0)
    out['clicks'] = np.stack([c, np.roll(c, 2000)], 1)
    # panned noise bursts with changing inter-channel balance (stereo weights)
    g = 0.5 + 0.5 * np.sin(2 * np.pi * 0.7 * t)
    nz = 0.3 * rng.standard_normal(n)
    hp = np.convolve(nz, [1, -0.9], 'same')
    out['panned'] = np.stack([hp * g, hp * (1 - g) + 0.01 * rng.standard_normal(n)], 1)
    # DC step and low-frequency rumble, length leaving a long final block (>= 886 samples)
    m = n + 900
    tt = np.arange(m) / SR
    d = 0.2 * (tt > 1.0) + 0.05 * np.sin(2 * np.pi * 25 * tt)
    out['dc_odd'] = np.stack([d, d * -0.5], 1)
    # exact multiple of the block size after the tool's 69-sample lead-in
    m = 1024 * 100 + 69
    tt = np.arange(m) / SR
    s = 0.2 * np.sin(2 * np.pi * 440 * tt) * np.sin(2 * np.pi * 3 * tt)
    out['exact_len'] = np.stack([s, s * 0.3], 1)
    # harmonic comb with one partial per 8-bin subband (f0 = 8 MDCT bins = 172.27 Hz)
    f0 = 8 * SR / 2048.0
    comb8 = sum(np.sin(2 * np.pi * f0 * k * t + 0.3 * k * k) for k in range(1, 40)) * 0.02
    out['comb8'] = np.stack([comb8, comb8 * 0.9], 1)
    # loud tones at the top of the coded spectrum (around 16.5 kHz, coefficients ~1020)
    hf = 0.7 * np.sin(2 * np.pi * 16560 * t) + 0.2 * np.sin(2 * np.pi * 16700 * t + 1)
    out['hf_tone'] = np.stack([hf, hf * 0.8], 1)
    # broadband bursts of very different spectral shapes (allocation under pressure)
    bb = rng.standard_normal(n)
    env = (np.sin(2 * np.pi * 5 * t) > 0) * 0.8 + 0.05
    tone_mix = sum(0.05 * np.sin(2 * np.pi * f * t) for f in (311, 1777, 4123, 9001, 13013))
    out['pressure'] = np.stack([np.clip(bb * env * 0.6 + tone_mix, -1, 1),
                                np.clip(np.roll(bb, 333) * env * 0.5 - tone_mix, -1, 1)], 1)
    return out


if __name__ == '__main__':
    d = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(d, exist_ok=True)
    for k, v in signals().items():
        save(os.path.join(d, k + '.wav'), v)
        print(k)
