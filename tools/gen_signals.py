"""Generate deterministic 16-bit 44.1 kHz stereo test signals that exercise ATRAC3 features."""
import sys, os, wave, numpy as np

SR = 44100

def save(path, x):
    x = np.clip(np.round(x * 32767.0), -32768, 32767).astype('<i2')
    with wave.open(path, 'wb') as w:
        w.setnchannels(x.shape[1]); w.setsampwidth(2); w.setframerate(SR)
        w.writeframes(x.tobytes())

def signals(dur=3.0):
    rng = np.random.default_rng(1234)
    t = np.arange(int(SR * dur)) / SR
    out = {}
    # chirps + amplitude modulation + noise bursts (L), descending chirp + noise (R)
    l = 0.3 * np.sin(2 * np.pi * (200 + 3000 * t) * t) * (1 + 0.5 * np.sin(2 * np.pi * 3 * t))
    l += 0.2 * rng.standard_normal(len(t)) * ((t % 0.5) < 0.05)
    r = 0.25 * np.sin(2 * np.pi * (5000 - 1500 * t) * t) + 0.05 * rng.standard_normal(len(t))
    out['mix'] = np.stack([l, r], 1)
    # sharp transients over quiet tone (gain control)
    x = 0.05 * np.sin(2 * np.pi * 440 * t)
    for k in range(int(dur * 4)):
        s = int((0.1 + 0.25 * k) * SR)
        n = np.arange(2000)
        x[s:s + 2000] += 0.8 * np.exp(-n / 150.0) * rng.standard_normal(2000)
    out['transients'] = np.stack([x, x * 0.7], 1)
    # stable tones (tonal components), harmonics
    x = sum(0.12 / (h + 1) * np.sin(2 * np.pi * 523.25 * (h + 1) * t + h) for h in range(12))
    y = sum(0.1 * np.sin(2 * np.pi * f * t) for f in (1000, 2500, 7000, 12000))
    out['tones'] = np.stack([x, y], 1)
    # highly correlated stereo (joint stereo)
    m = 0.3 * np.sin(2 * np.pi * 330 * t) + 0.1 * rng.standard_normal(len(t))
    out['corr'] = np.stack([m, 0.9 * m + 0.02 * rng.standard_normal(len(t))], 1)
    # white noise near full scale, and silence-then-burst
    out['noise'] = 0.5 * rng.standard_normal((len(t), 2)).clip(-2, 2)
    x = np.zeros((len(t), 2)); x[len(t) // 2:] = 0.3 * rng.standard_normal((len(t) - len(t) // 2, 2))
    out['silence_burst'] = x
    return out

if __name__ == '__main__':
    d = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(d, exist_ok=True)
    for k, v in signals().items():
        save(os.path.join(d, k + '.wav'), v)
        print(k)
