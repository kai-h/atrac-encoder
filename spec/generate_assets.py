"""Generate/check the portable spec assets. Needs Python and NumPy, no Sony binary.

Run from any directory: python spec/generate_assets.py --write | --check
Optional maintainer check: --verify-sony (needs tools/at3emu.py and private binary).
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import wave

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'reference'))
import at3_tables
import lp2_tables
import lp4_tables
import lp2_analysis
from lp2_alloc import order_descending
import atrac3_encoder_lp2 as lp2
import atrac3_encoder_lp4 as lp4


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def floating(value, width=32):
    """Bit patterns, not decimal approximations; row-major flattening."""
    a = np.asarray(value, dtype=np.float32 if width == 32 else np.float64)
    u = a.view(np.uint32 if width == 32 else np.uint64).reshape(-1)
    return {'type': f'binary{width}', 'shape': list(a.shape),
            'bits': [f'{int(v):0{width // 4}x}' for v in u]}


def portable(value):
    if isinstance(value, dict):
        # Explicit keys retain tuple-valued codebook identifiers and string bit codes.
        return {'type': 'map', 'entries': [
            {'key': list(k) if isinstance(k, tuple) else k, 'value': portable(v)}
            for k, v in sorted(value.items())]}
    if isinstance(value, (float, np.floating)):
        return floating(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (list, tuple, np.ndarray)):
        a = np.asarray(value)
        if a.dtype.kind == 'f':
            return floating(a)
        if a.dtype.kind in 'iu':
            return {'type': 'integer', 'shape': list(a.shape),
                    'values': [int(v) for v in a.reshape(-1)]}
        return [portable(v) for v in value]
    raise TypeError(f'Unsupported export value: {type(value)}')


def tables(module):
    return {name: portable(value) for name, value in vars(module).items()
            if name.isupper() and not name.startswith('_')}


def pcm_signals():
    """Integer-only recipes: no RNG version or transcendental-function dependency."""
    n = 8192
    i = np.arange(n, dtype=np.int64)
    tri = ((i % 64) - 32) * 768
    # A deterministic 32-bit LCG, emitting the signed upper 16 bits.
    state = 0x12345678
    noise = []
    for _ in range(n * 2):
        state = (1664525 * state + 1013904223) & 0xffffffff
        noise.append((state >> 16) - 32768)
    noise = np.array(noise, dtype=np.int16).reshape(n, 2)
    transient = noise.copy()
    transient[:2048] //= 128
    transient[4096:6144] //= 64
    transient[::257, 0] = 32767
    transient[::263, 1] = -32768
    pan = np.stack([tri, tri], axis=1)
    pan[:2048, 1] = 0
    pan[2048:4096, 0] = 0
    pan[4096:6144, 1] //= 16
    pan[6144:, 1] = -pan[6144:, 0]
    # Quantised smooth sinusoid, specified by a fixed integer lookup table.
    lut = np.array([0, 6243, 12246, 17778, 22627, 26607, 29564, 31385,
                    32000, 31385, 29564, 26607, 22627, 17778, 12246, 6243,
                    0, -6243, -12246, -17778, -22627, -26607, -29564, -31385,
                    -32000, -31385, -29564, -26607, -22627, -17778, -12246, -6243])
    tone = lut[i % 32]
    yield 'silence', np.zeros((n, 2), np.int16), 'Silence and minimal unit headers.'
    yield 'tone', np.stack([tone, tone], axis=1), 'Correlated periodic tone and tonal extraction.'
    yield 'quiet', np.stack([(i % 3) - 1, tone // 1024], axis=1), 'Near-silence and weak tones.'
    yield 'transients', transient, 'Noise, attacks/releases, clicks, and full-range int16 PCM.'
    yield 'stereo', pan, 'Left/right dominance, imbalance, correlation, and opposite polarity.'
    yield 'clipping', np.stack([np.where(i % 32 < 16, 32767, -32768),
                               noise[:, 1]], axis=1), 'Full-scale square wave and broadband pressure.'
    for remainder in (0, 885, 886, 1023):
        length = 7 * 1024 - 955 + remainder
        yield f'tail_{remainder}', noise[:length], f'(N + 955) mod 1024 = {remainder}.'


def wav_bytes(pcm):
    buf = io.BytesIO()
    with wave.open(buf, 'wb') as w:
        w.setparams((2, 2, 44100, len(pcm), 'NONE', 'not compressed'))
        w.writeframes(np.asarray(pcm, dtype='<i2').tobytes())
    return buf.getvalue()


def data_chunk(at3):
    pos = 12
    while pos + 8 <= len(at3):
        size = int.from_bytes(at3[pos + 4:pos + 8], 'little')
        if at3[pos:pos + 4] == b'data':
            return at3[pos + 8:pos + 8 + size]
        pos += 8 + size + (size & 1)
    raise ValueError('Missing data chunk')


def gain_info(points):
    return [{'points_level_location': [list(p) for p in g.pts], 'peak': floating(g.peak)}
            for g in points]


def stages(pcm):
    """Three unprimed raw-library calls. Float checkpoints precede allocation."""
    a = lp2_analysis.ChannelAnalysis()
    lp2_calls = []
    for k in range(3):
        A, B, gain, prev = a.process(pcm[k * 1024:(k + 1) * 1024, 0].astype(np.float32))
        lp2_calls.append({'call': k, 'qmf_bands': floating(np.stack([h[512:] for h in a.hist])),
                          'spectrum_A': floating(A), 'spectrum_B': floating(B),
                          'gain': gain_info(gain), 'previous_gain': gain_info(prev)})
    enc = lp4.Encoder(66)
    lp4_calls = []
    original_qmf, original_mdct = lp4.qmf_analysis, lp4.mdct
    qmf_out, mdct_out = [], []

    def qmf(*args):
        out = original_qmf(*args)
        qmf_out.append(floating(out[0]))
        return out

    def mdct(*args):
        out = original_mdct(*args)
        mdct_out.append(floating(out[0]))
        return out

    try:
        lp4.qmf_analysis, lp4.mdct = qmf, mdct
        for k in range(3):
            qmf_out.clear(); mdct_out.clear()
            frame = enc.encode(pcm[k * 1024:(k + 1) * 1024])
            lp4_calls.append({'call': k, 'qmf_before_stereo': list(qmf_out),
                              'spectrum': list(mdct_out),
                              'gain': [[[int(g.lev[j]), int(g.loc[j])] for j in range(g.count)]
                                       for u in enc.units for g in u.gain],
                              'sent_matrix': list(enc.stereo.dec_prev),
                              'sent_weight': enc.stereo.w_prev,
                              'budgets': [u.budget for u in enc.units],
                              'frame_hex': frame.hex() if frame is not None else None})
    finally:
        lp4.qmf_analysis, lp4.mdct = original_qmf, original_mdct
    return {'format_version': 1, 'input': 'transients.wav', 'lead_in_samples': 0,
            'lp2_channel': 0, 'lp2': lp2_calls, 'lp4': lp4_calls}


def build_assets():
    assets = {}
    common = tables(at3_tables)
    t2, t4 = tables(lp2_tables), tables(lp4_tables)
    # Constants encoded in helper functions or arithmetic expressions, rather than tables.
    t2['OFFSET_FACTORS'] = floating([lp2_tables.offset_divisor(i) for i in range(32)], 64)
    t2['NEGATIVE_OFFSET_WEIGHTS'] = floating([lp2_tables.offset_step_weight(i) or 1.0
                                             for i in range(32)], 64)
    t2['DENSITY_THRESHOLDS'] = floating([1.15, 1.1, 1.05, 1.0, 0.95, 0.9], 64)
    t2['DOUBLE_CONSTANTS'] = {str(v): floating(v, 64) for v in
                             (1.5, 0.5, 0.7, 0.1, 1.25, 1e-20, 10.0, 31.5)}
    t2['FLOAT_CONSTANTS'] = {'release_ratio': floating(1.85), 'unused_margin': floating(-9.99)}
    t4['DOUBLE_CONSTANTS'] = {name: floating(v, 64) for name, v in
                             [('side_ratio', 3.94), ('side_slew', 0.05),
                              ('side_limit', 5 / 3), ('magic', 12582912.0)]}
    t4['FLOAT_CONSTANTS'] = {'attack_ratio': floating(1.6), 'energy_seed': floating(0.1)}
    for name, value in [('common', common), ('lp2', t2), ('lp4', t4)]:
        assets[f'tables/{name}.json'] = json_bytes({'format_version': 1, 'tables': value})
    ordering = []
    for keys in ([0] * 29, [3, 1, 3, 2, 1, 3, 2, 0],
                 [i % 4 for i in range(32)], [(i * 7) % 5 for i in range(29)]):
        ordering.append({'keys': keys, 'indices': order_descending(keys)})
    assets['vectors/ordering.json'] = json_bytes({'format_version': 1, 'cases': ordering})
    # Hand-specified cases distinguish interval filling from summing suffix sums twice.
    assets['vectors/gain-profiles.json'] = json_bytes({'format_version': 1, 'cases': [
        {'attacks': [[3, 5], [7, 2]], 'releases': [],
         'profile': [5] * 4 + [2] * 4 + [0] * 25, 'points_level_location': [[9, 3], [6, 7]]},
        {'attacks': [], 'releases': [[12, 4], [4, 2]],
         'profile': [0] * 4 + [2] * 8 + [4] * 21, 'points_level_location': [[0, 3], [2, 11]]},
        {'attacks': [[3, 5], [7, 2]], 'releases': [[12, 4], [4, 2]],
         'profile': [5] * 4 + [4] * 4 + [2] * 4 + [4] * 21,
         'points_level_location': [[5, 3], [4, 7], [2, 11]]}]})
    cases = []
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        for name, pcm, description in pcm_signals():
            pcm = np.asarray(pcm, np.int16)
            wav = wav_bytes(pcm)
            source = root / 'input.wav'
            source.write_bytes(wav)
            assets[f'vectors/{name}.wav'] = wav
            modes = []
            for mode, size in [('lp2', 384), ('lp4', 192)]:
                dst = root / 'output.at3'
                if mode == 'lp2':
                    lp2.encode_file(str(source), str(dst))
                else:
                    lp4.encode_wav(str(source), str(dst), 66)
                encoded = dst.read_bytes()
                raw = data_chunk(encoded)
                D = (len(pcm) + 955 + 1023) // 1024
                r = (len(pcm) + 955) % 1024
                count = D if 0 < r < 886 else D + 1
                assert len(raw) == count * size, (name, mode, len(raw), count)
                assets[f'vectors/{name}.{mode}.at3'] = encoded
                assets[f'vectors/{name}.{mode}.frames'] = raw
                modes.append({'mode': mode, 'frame_bytes': size, 'frame_count': count,
                              'at3': f'{name}.{mode}.at3', 'raw_frames': f'{name}.{mode}.frames',
                              'frame_sha256': [sha(raw[i:i + size]) for i in range(0, len(raw), size)]})
            cases.append({'name': name, 'description': description, 'pcm': f'{name}.wav',
                          'samples_per_channel': len(pcm), 'final_data_block_samples': r or 1024,
                          'empty_analysis_calls': 2 if 0 < r < 886 else 3, 'modes': modes})
            if name == 'transients':
                assets['vectors/stages.json'] = json_bytes(stages(pcm))
    assets['vectors/cases.json'] = json_bytes({'format_version': 1, 'sample_rate': 44100,
                                             'channels': 2, 'pcm_bits': 16, 'cases': cases})
    assets['assets-manifest.json'] = json_bytes({
        'format_version': 1, 'generator': 'generate_assets.py',
        'provenance': 'Synthetic PCM; expected outputs/checkpoints from spec/reference. '
                      'Sony verification, when available, is recorded separately.',
        'files': {name: {'bytes': len(data), 'sha256': sha(data)}
                  for name, data in sorted(assets.items())}})
    return assets


def verify_sony():
    """Compare stored complete files; don't silently overwrite expected outputs."""
    manifest = json.loads((HERE / 'assets-manifest.json').read_text())
    for name, info in manifest['files'].items():
        if sha((HERE / name).read_bytes()) != info['sha256']:
            raise ValueError(f'Asset hash mismatch: {name}')
    sys.path.insert(0, str(HERE.parent / 'tools'))
    from at3emu import Emu
    results = []
    with tempfile.TemporaryDirectory() as d:
        for case in json.loads((HERE / 'vectors/cases.json').read_text())['cases']:
            for mode in case['modes']:
                output = Path(d) / 'sony.at3'
                output.unlink(missing_ok=True)
                code = Emu(quiet=True).run_main(['-e', str(HERE / 'vectors' / case['pcm']),
                                               str(output), '-br',
                                               '132' if mode['mode'] == 'lp2' else '66'])
                expected = HERE / 'vectors' / mode['at3']
                if not output.exists() or output.read_bytes() != expected.read_bytes():
                    raise ValueError(f'Sony mismatch: {expected.name} (exit {code})')
                results.append({'file': f'vectors/{mode["at3"]}', 'sha256': sha(expected.read_bytes())})
                print(f'Sony identical: {expected.name}', flush=True)
    report = {'format_version': 1, 'oracle': 'psp_at3tool.exe 3.0.0.0 under tools/at3emu.py',
              'scope': 'Complete non-looping .at3 files; intermediate checkpoints are model-derived.',
              'files': results}
    (HERE / 'vectors/sony-verification.json').write_bytes(json_bytes(report))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    group.add_argument('--verify-sony', action='store_true')
    args = parser.parse_args()
    if args.verify_sony:
        verify_sony()
        return
    assets = build_assets()
    mismatches = []
    for name, data in assets.items():
        dst = HERE / name
        if args.write:
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(data)
        elif not dst.exists() or dst.read_bytes() != data:
            mismatches.append(name)
    if mismatches:
        sys.exit('Assets differ: ' + ', '.join(mismatches))
    print(f'{"Wrote" if args.write else "Verified"} {len(assets)} assets; '
          '20 complete LP2/LP4 files, raw frames, exact tables, ordering and stage checkpoints.')


if __name__ == '__main__':
    main()
