"""Chapter 8 §§4, 10.5–10.6, and 11: bit-based scale indices and quantisers."""
import numpy as np

from .bitstream import Tone, codewords
from .lp4_stereo import T, f32
from .lp4_gain import float_bits

MAGIC = 12582912.0


def scale_index4(values):
    bits = np.asarray(values, np.float32).view(np.uint32)
    magnitude = max((int(v) << 1) & 0xffffffff for v in bits)
    exponent, mantissa = magnitude >> 24, magnitude & 0xffffff
    index = 3 * exponent - 364
    if mantissa > 0x965fe9:
        index += 1
    if mantissa < 0x428a30:
        index -= 1
    return index if 0 <= index < 64 else 0


def quantiser_scale(scale, selector):
    if not 0 <= scale < 64 or not 1 <= selector <= 7:
        raise ValueError('Invalid quantiser scale/selector')
    bits = int(T['QSCALE_BITS'][3 * selector - scale % 3]) - ((scale // 3) << 23)
    return float(np.array(bits, np.uint32).view(np.float32))


def quantise(values, scale, selector):
    product = np.asarray(values, np.float32) * np.float32(quantiser_scale(scale, selector))
    rounded = product + np.float32(MAGIC)
    return rounded.view(np.uint32).astype(np.int64) - 0x4b400000


def bit_cost(values, scale, selector):
    q = quantise(values, scale, selector)
    return 6 + sum(length for _, length in codewords(q, selector, 0, T))


def signed_low16(value):
    bits = float_bits(value) & 0xffff
    return bits - 65536 if bits >= 32768 else bits


def extract_tone(spectrum, position, selector):
    values = spectrum[position:position + 4].astype(np.float64)
    lo = scale_index4(values)
    hi = min(lo + 3, 64) - 1
    best, scale = 1e10, lo
    for candidate in range(hi, lo - 1, -1):
        sc = quantiser_scale(candidate, selector)
        error, accumulator = 0.0, 0
        for x in values:
            y = sc * float(x)
            q = signed_low16(y + MAGIC)
            accumulator |= q
            if abs(y - q) > error:
                error = f32(abs(y - q))
        if sc * best > error:
            best, scale = error / sc, candidate
            while accumulator % 2 == 0 and scale < 61:
                accumulator >>= 1
                scale += 3
    sc = quantiser_scale(scale, selector)
    inverse = f32(1 / sc)
    q = [signed_low16(sc * float(x) + MAGIC) for x in values]
    spectrum[position:position + 4] = [f32(float(x) - value * inverse) for x, value in zip(values, q)]
    masked = [value & int(T['QMASK'][selector]) for value in q]
    tone = Tone(position, scale, masked)
    cost = 12 + sum(length for _, length in codewords(masked, selector, 0, T))
    return tone, cost


def restore_tone(spectrum, tone, selector):
    inverse = 1 / quantiser_scale(tone.scale, selector)
    width = [0, 2, 3, 3, 4, 4, 5, 6][selector]
    for j, value in enumerate(tone.values):
        if value & (1 << (width - 1)):
            value -= 1 << width
        spectrum[tone.position + j] = f32(value * inverse + float(spectrum[tone.position + j]))
