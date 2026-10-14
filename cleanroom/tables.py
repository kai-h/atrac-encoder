"""Read the portable data format defined in spec/tables/README.md."""
import json
from pathlib import Path

import numpy as np

SPEC = Path(__file__).resolve().parents[1] / 'spec'


def decode(value):
    if isinstance(value, list):
        return [decode(v) for v in value]
    if not isinstance(value, dict):
        return value
    kind = value.get('type')
    if kind in ('binary32', 'binary64'):
        width = 32 if kind == 'binary32' else 64
        data = np.array([int(v, 16) for v in value['bits']], dtype=f'uint{width}')
        return data.view(f'float{width}').reshape(value['shape'])
    if kind == 'integer':
        return np.array(value['values'], dtype=np.int64).reshape(value['shape'])
    if kind == 'map':
        return {(tuple(e['key']) if isinstance(e['key'], list) else e['key']): decode(e['value'])
                for e in value['entries']}
    return {k: decode(v) for k, v in value.items()}


def load(mode):
    document = json.loads((SPEC / 'tables' / f'{mode}.json').read_text())
    if document['format_version'] != 1:
        raise ValueError('Unsupported specification table version')
    return decode(document['tables'])
