"""Small syntax inspector for debugging encoder conformance, ATRAC3.md §4.

This parses only codec fields; it performs no audio synthesis.
"""
from .tables import load

COMMON = load('common')


class Reader:
    def __init__(self, data):
        self.data, self.position = data, 0

    def read(self, count):
        if self.position + count > 8 * len(self.data):
            raise ValueError('Bitstream underflow')
        value = 0
        for _ in range(count):
            value = 2 * value + ((self.data[self.position // 8] >> (7 - self.position % 8)) & 1)
            self.position += 1
        return value

    def mantissas(self, selector, mode, count):
        codes = COMMON['CODES'][mode, selector]
        result = []
        while len(result) < count:
            code = ''
            while code not in codes:
                code += str(self.read(1))
                if len(code) > 8:
                    raise ValueError('Invalid mantissa code')
            value = codes[code]
            result.extend(value.tolist() if hasattr(value, 'tolist') else [value])
        return result


def unit(data, joint=False):
    reader = Reader(data)
    header = None
    if joint:
        weight = reader.read(4)
        header = weight, [reader.read(2) for _ in range(4)]
    if reader.read(2 if joint else 6) != (3 if joint else 0x28):
        raise ValueError('Invalid sound-unit ID')
    bands = reader.read(2) + 1
    gains = []
    for _ in range(bands):
        gains.append([(reader.read(4), reader.read(5)) for _ in range(reader.read(3))])
    groups = []
    count = reader.read(5)
    tone_mode = reader.read(2) if count else None
    for _ in range(count):
        flags = [reader.read(1) for _ in range(bands)]
        values, selector = reader.read(3) + 1, reader.read(3)
        mode = reader.read(1) if tone_mode == 3 else tone_mode
        records = []
        for band in range(bands):
            if flags[band]:
                for cell in range(4):
                    for _ in range(reader.read(3)):
                        scale = reader.read(6)
                        position = 256 * band + 64 * cell + reader.read(6)
                        q = reader.mantissas(selector, mode, min(values, 1024 - position))
                        records.append((position, scale, q))
        groups.append({'selector': selector, 'values': values, 'mode': mode, 'flags': flags, 'records': records})
    nsub = reader.read(5) + 1
    mode = reader.read(1)
    selectors = [reader.read(3) for _ in range(nsub)]
    scales = [reader.read(6) if selector else None for selector in selectors]
    mantissas = []
    edges = COMMON['SUBBAND_START']
    for i, selector in enumerate(selectors):
        mantissas.append(reader.mantissas(selector, mode, int(edges[i + 1] - edges[i])) if selector else [])
    return {'header': header, 'bands': bands, 'gains': gains, 'groups': groups, 'selectors': selectors,
            'scales': scales, 'mantissas': mantissas, 'bits': reader.position}
