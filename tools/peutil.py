"""Tiny helpers for reading data from psp_at3tool.exe by virtual address."""
import struct, pefile, os
_PE_PATH = os.path.join(os.path.dirname(__file__), '..', 'private', 'psp_at3tool.exe')
pe = pefile.PE(_PE_PATH, fast_load=True)
BASE = pe.OPTIONAL_HEADER.ImageBase

def read(va, n):
    return pe.get_data(va - BASE, n)

def u32s(va, n):
    return list(struct.unpack('<%dI' % n, read(va, 4 * n)))

def i32s(va, n):
    return list(struct.unpack('<%di' % n, read(va, 4 * n)))

def i16s(va, n):
    return list(struct.unpack('<%dh' % n, read(va, 2 * n)))

def u8s(va, n):
    return list(read(va, n))

def f32s(va, n):
    return list(struct.unpack('<%df' % n, read(va, 4 * n)))

def f64s(va, n):
    return list(struct.unpack('<%dd' % n, read(va, 8 * n)))
