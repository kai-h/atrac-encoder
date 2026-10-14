"""ATRAC3.md §4: shared sound-unit and frame serialisation."""
from dataclasses import dataclass, field

@dataclass
class Tone:
    position: int
    scale: int
    values: list


@dataclass
class ToneGroup:
    selector: int = 7
    value_count: int = 4
    mode: int = 1
    bands: set = field(default_factory=set)
    cells: dict = field(default_factory=dict)


class Bits:
    def __init__(self):
        self.bytes = bytearray()
        self.count = 0

    def put(self, value, width):
        if width < 0 or value < 0 or value >= 1 << width:
            raise ValueError(f'Value {value} does not fit an unsigned {width}-bit field')
        for bit in range(width - 1, -1, -1):
            if self.count % 8 == 0:
                self.bytes.append(0)
            self.bytes[-1] |= ((value >> bit) & 1) << (7 - self.count % 8)
            self.count += 1


def codewords(values, selector, mode, tables):
    codes = tables['ENC_CODES'] if 'ENC_CODES' in tables else tables['VLC']
    width = [0, 2, 3, 3, 4, 4, 5, 6][selector]
    mask = (1 << width) - 1
    if selector == 1:
        if len(values) % 2:
            raise ValueError('Selector 1 requires complete mantissa pairs')
        indices = [((int(values[i]) & 3) << 2) | (int(values[i + 1]) & 3)
                   for i in range(0, len(values), 2)]
    else:
        indices = [int(v) & mask for v in values]
    for index in indices:
        code, length = (int(v) for v in codes[mode, selector][index])
        if length <= 0:
            raise ValueError('Unspecified mantissa code')
        yield code, length


def pack_unit(nbands, gains, groups, selectors, scales, mantissas, tables, joint_header=None):
    bits = Bits()
    if joint_header is None:
        bits.put(0x28, 6)
    else:
        weight, matrices = joint_header
        bits.put(weight, 4)
        for matrix in matrices:
            if matrix == 2:
                raise ValueError('Invalid joint-stereo selector 2')
            bits.put(matrix, 2)
        bits.put(3, 2)
    if not 1 <= nbands <= 4 or not 1 <= len(selectors) <= 32:
        raise ValueError('Invalid band or subband count')
    bits.put(nbands - 1, 2)
    for band in range(nbands):
        points = gains[band].points
        if len(points) > 7 or any(a[1] >= b[1] for a, b in zip(points, points[1:])):
            raise ValueError('Invalid gain locations/count')
        bits.put(len(points), 3)
        for level, location in points:
            bits.put(level, 4)
            bits.put(location, 5)
    bits.put(len(groups), 5)
    if groups:
        modes = {g.mode for g in groups}
        tone_mode = next(iter(modes)) if len(modes) == 1 else 3
        bits.put(tone_mode, 2)
        count = 0
        for group in groups:
            if not group.bands.issubset(set(range(nbands))):
                raise ValueError('Tonal band lies outside the unit header')
            for band in range(nbands):
                bits.put(int(band in group.bands), 1)
            bits.put(group.value_count - 1, 3)
            bits.put(group.selector, 3)
            if tone_mode == 3:
                bits.put(group.mode, 1)
            for band in range(nbands):
                if band not in group.bands:
                    continue
                for cell in range(4 * band, 4 * band + 4):
                    tones = group.cells.get(cell, [])
                    bits.put(len(tones), 3)
                    for tone in tones:
                        count += 1
                        bits.put(tone.scale, 6)
                        bits.put(tone.position % 64, 6)
                        for code, width in codewords(tone.values, group.selector, group.mode, tables):
                            bits.put(code, width)
        if count > 64:
            raise ValueError('More than 64 tonal components')
    bits.put(len(selectors) - 1, 5)
    bits.put(0, 1)  # Both encoder chapters specify VLC spectrum coding.
    for selector in selectors:
        bits.put(int(selector), 3)
    for selector, scale in zip(selectors, scales):
        if selector:
            bits.put(int(scale), 6)
    for selector, values in zip(selectors, mantissas):
        if selector:
            for code, width in codewords(values, selector, 0, tables):
                bits.put(code, width)
    return bytes(bits.bytes), bits.count


def independent_frame(units, unit_bytes):
    out = bytearray()
    for unit in units:
        if len(unit) > unit_bytes:
            raise ValueError('Sound unit exceeds its fixed bit budget')
        out.extend(unit)
        out.extend(bytes(unit_bytes - len(unit)))
    return bytes(out)


def joint_frame(unit0, unit1, size=192):
    if len(unit0) + len(unit1) > size:
        raise ValueError('Joint-stereo sound units overlap')
    return unit0 + bytes(size - len(unit0) - len(unit1)) + unit1[::-1]
