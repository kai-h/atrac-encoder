"""Check frozen assets and optionally independent encoder output. Standard library only."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outputs', type=Path, help='Directory containing expected output basenames')
    parser.add_argument('--raw', action='store_true', help='Compare .frames instead of .at3')
    args = parser.parse_args()
    if args.raw and args.outputs is None:
        parser.error('--raw requires --outputs')
    manifest = json.loads((HERE / 'assets-manifest.json').read_text())
    failures = []
    for name, info in manifest['files'].items():
        path = HERE / name
        if not path.exists():
            failures.append(f'Missing asset: {name}')
            continue
        data = path.read_bytes()
        if len(data) != info['bytes'] or hashlib.sha256(data).hexdigest() != info['sha256']:
            failures.append(f'Asset hash/size mismatch: {name}')
    if failures:
        sys.exit('\n'.join(failures))
    print(f'Asset integrity: {len(manifest["files"])} files passed.')
    compared = 0
    if args.outputs is not None:
        cases = json.loads((HERE / 'vectors/cases.json').read_text())['cases']
        for case in cases:
            for mode in case['modes']:
                name = mode['raw_frames' if args.raw else 'at3']
                path = args.outputs / name
                if not path.exists():
                    failures.append(f'Missing output: {path}')
                    continue
                expected = (HERE / 'vectors' / name).read_bytes()
                actual = path.read_bytes()
                if actual != expected:
                    offset = next((i for i, (a, b) in enumerate(zip(actual, expected)) if a != b),
                                  min(len(actual), len(expected)))
                    position = f'byte {offset}'
                    if args.raw:
                        position += f' (frame {offset // mode["frame_bytes"]}, '
                        position += f'byte {offset % mode["frame_bytes"]} within frame)'
                    failures.append(f'{name}: differs at {position}; '
                                    f'actual/expected lengths {len(actual)}/{len(expected)}')
                else:
                    compared += 1
        if failures:
            sys.exit('\n'.join(failures))
        print(f'Encoder outputs: {compared} files passed.')


if __name__ == '__main__':
    main()
