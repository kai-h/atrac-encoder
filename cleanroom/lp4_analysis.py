"""Chapter 8 §§5 and 8: independent LP4 QMF and MDCT kernels."""
import numpy as np

from .tables import load

T = load('lp4')


class QMF:
    def __init__(self):
        self.input_history = np.zeros(46, np.float32)
        self.vector_history = np.zeros((23, 4), np.float32)

    def process(self, pcm):
        if np.shape(pcm) != (1024,):
            raise ValueError('QMF requires exactly 1024 samples')
        h = np.concatenate([self.input_history, np.asarray(pcm, np.float32)])
        vectors = np.empty((279, 4), np.float32)
        vectors[:23] = self.vector_history
        result = np.empty((256, 4), np.float32)
        q1, q2 = T['QMF1'], T['QMF2']
        for i in range(256):
            base = 40 + 4 * i
            a = q1[11] * h[base + 4:base + 8]
            b = q1[11] * h[base + 6:base + 10]
            for k in range(11):
                offset = -40 + 4 * k
                a = a + q1[k] * h[base + offset:base + offset + 4]
                b = b + q1[k] * h[base + offset + 2:base + offset + 6]
            u = np.array([a[2] + a[0], b[2] + b[0], a[3] + a[1], b[3] + b[1]], np.float32)
            vectors[23 + i] = [u[0] + u[2], u[1] + u[3], u[2] - u[0], u[3] - u[1]]
            w = q2[23] * vectors[23 + i]
            for k in range(23):
                w = w + q2[k] * vectors[i + k]
            result[i] = [w[0] + w[1], w[1] - w[0], w[3] - w[2], w[2] + w[3]]
        self.input_history = h[-46:].copy()
        self.vector_history = vectors[-23:].copy()
        return result


class MDCT:
    def __init__(self):
        self.state = np.zeros((128, 4), np.float32)

    def process(self, previous_samples, first_levels):
        if np.shape(previous_samples) != (256, 4) or len(first_levels) != 4:
            raise ValueError('MDCT requires four 256-sample bands and four initial gain levels')
        y = np.empty((256, 4), np.float32)
        next_state = np.empty_like(self.state)
        for j, p in enumerate(T['MDCT_PERM']):
            r0, r1, r2, r3 = T['MDCT_TW'][p]
            a, b = previous_samples[128 + p], previous_samples[127 - p]
            s = self.state[j] * T['GAIN_POW2'][first_levels]
            next_state[j] = r1 * (a * r0 - b)
            t = r2 * (b * r0 + a)
            d = r3 * (t - s)
            y[2 * j], y[2 * j + 1] = (t + s) + d, d
        factors = T['MDCT_FFT_TW']
        for j in range(64):
            a, b, c, d = y[4 * j:4 * j + 4].copy()
            f = factors[2 * j]
            e3, e0 = f * (b - d), f * (a - c)
            y[4 * j:4 * j + 4] = [(a + c) + e3, (b + d) + e0, e0, e3]
        for n in (4, 8, 16, 32, 64):
            for k in range(256 // (2 * n)):
                start = 2 * n * k
                low, high = y[start:start + n].copy(), y[start + n:start + 2 * n].copy()
                high_new = factors[n // 2 - 1 + n * k] * (low - high)
                y[start:start + n] = (low + high) + high_new[::-1]
                y[start + n:start + 2 * n] = high_new
        low, high = y[:128].copy(), y[128:].copy()
        y[:128], y[128:] = low + high, T['MDCT_SQRT_HALF'] * (low - high)
        result = np.empty((4, 256), np.float32)
        for band in range(4):
            z = y[:128, band] + y[255:127:-1, band]
            tail = y[255:127:-1, band]
            result[band, :128] = tail if band % 2 else z
            result[band, 255:127:-1] = z if band % 2 else tail
        self.state = next_state
        return result
