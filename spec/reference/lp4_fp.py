"""Floating-point helpers for the LP4 encoder model.

Sony's encoder mixes two kinds of arithmetic and the model has to follow both exactly:
  * SSE single precision (QMF, envelope, MDCT, bit counting, final quantisation): every
    operation rounds to float32.  Modelled with numpy float32 arrays/scalars.
  * x87 code with the precision control set to 53 bits (joint stereo, gain control, tone
    coding, rate split): intermediates are IEEE doubles and are rounded to float32 only when
    the program stores them into a float variable.  Modelled with Python floats plus f32().
"""
import struct
import numpy as np

F32 = np.float32
MAGIC = 12582912.0          # 1.5 * 2**23: adding it rounds a float32 to an integer (bits hold the result)


def f32(x):
    """Round a Python float to the nearest float32 value (ties to even), returned as float."""
    return float(np.float32(x))


def bits(x):
    """IEEE single bit pattern of x (rounded to float32 first)."""
    return struct.unpack('<I', struct.pack('<f', x))[0]


def from_bits(u):
    return struct.unpack('<f', struct.pack('<I', u & 0xffffffff))[0]


def low16(x):
    """Signed low 16 bits of the float32 pattern of x: the integer produced by the MAGIC trick."""
    v = bits(x) & 0xffff
    return v - 0x10000 if v & 0x8000 else v
