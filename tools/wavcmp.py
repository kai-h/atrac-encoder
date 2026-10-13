"""Compare two WAV files: best alignment offset, SNR, max abs diff."""
import sys, wave, numpy as np

def load(p):
    w = wave.open(p, 'rb')
    ch, n = w.getnchannels(), w.getnframes()
    x = np.frombuffer(w.readframes(n), dtype='<i2').astype(np.float64).reshape(-1, ch)
    return x

def compare(a, b, maxlag=4096):
    A, B = load(a), load(b)
    best = None
    for lag in range(-maxlag, maxlag + 1, 1):
        if lag >= 0: x, y = A[lag:], B
        else: x, y = A, B[-lag:]
        n = min(len(x), len(y), 44100)
        if n < 1000: continue
        x0, y0 = x[2048:n], y[2048:n]
        err = np.sum((x0 - y0) ** 2)
        if best is None or err < best[1]: best = (lag, err)
    lag = best[0]
    x, y = (A[lag:], B) if lag >= 0 else (A, B[-lag:])
    n = min(len(x), len(y))
    d = x[:n] - y[:n]
    snr = 10 * np.log10(np.sum(x[:n] ** 2) / max(np.sum(d ** 2), 1e-9))
    return lag, snr, np.max(np.abs(d)), n

if __name__ == '__main__':
    lag, snr, mx, n = compare(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4096)
    print('lag=%d  SNR=%.2f dB  maxdiff=%d  n=%d' % (lag, snr, mx, n))
