"""Regression test for the compiled encoder/decoder (standard library only).

Generates the test vectors (tests/gen_vectors.py), encodes each in every mode it supports,
decodes the result, and compares SHA-256 hashes with tests/expected.txt.  The expected
encoder output was verified byte-identical to the Python reference models (and so to
Sony's at3tool) with tools/verify_c.py before it was recorded.  Also checks the RIFF
parser against oversized chunk sizes (riff_checks), and with --wasm that the streaming
encoder fed in pieces (wasm_cli.mjs -split) matches the whole-file encoder byte for byte.

Usage:
  python tests/run_tests.py [--bin build/atrac3]        native command-line tool
  python tests/run_tests.py --wasm build-wasm           WebAssembly module, run with node
  python tests/run_tests.py --update                    rewrite expected.txt (after verify_c.py!)
"""
import argparse, hashlib, os, shutil, struct, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen_vectors  # noqa: E402

EXPECTED = os.path.join(HERE, 'expected.txt')
MODES = {2: [132, 105, 66], 1: [66, 52]}


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def riff_checks(tool, tmp, got):
    """Oversized chunk sizes once wrapped the 32-bit chunk walk in src/files.c, hanging the
    WebAssembly build; a streamed WAV's 0xFFFFFFFF sizes must still read as "to end of file".
    Returns a list of failures."""
    fails = []
    def run(args):
        try:
            return subprocess.run(tool + args, capture_output=True, text=True, timeout=30)
        except subprocess.TimeoutExpired:
            return None
    fmt = struct.pack('<4sIHHIIHH', b'fmt ', 16, 1, 2, 44100, 176400, 4, 16)
    for size in (0xfffffff7, 0xfffffff8):
        body = b'WAVE' + fmt + b'junk' + struct.pack('<I', size) + bytes(64)
        path = os.path.join(tmp, 'junk_%08x.wav' % size)
        with open(path, 'wb') as f:
            f.write(b'RIFF' + struct.pack('<I', len(body)) + body)
        for mode in ('-e', '-d'):
            r = run([mode, path, path + '.out'])
            if r is None or r.returncode == 0 or 'missing data chunk' not in r.stderr:
                fails.append('%s %s: %s' % (mode, os.path.basename(path), 'hung' if r is None else 'wrong result'))
    name = next(v[0] for v in gen_vectors.VECTORS if not v[0].endswith('_mono'))
    with open(os.path.join(tmp, name + '.wav'), 'rb') as f:
        wav = bytearray(f.read())
    wav[4:8] = wav[40:44] = b'\xff' * 4      # RIFF and data sizes as ffmpeg writes to a pipe
    path = os.path.join(tmp, 'streamed.wav')
    with open(path, 'wb') as f:
        f.write(wav)
    key = '%s.%d.at3' % (name, MODES[2][0])
    r = run(['-e', '-br', str(MODES[2][0]), path, path + '.at3'])
    if r is None or r.returncode != 0 or sha(path + '.at3') != got[key]:
        fails.append('-e streamed.wav: %s' % ('hung' if r is None else 'differs from ' + key))
    return fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', default=None, help='path of the atrac3 executable')
    ap.add_argument('--wasm', default=None, help='directory holding atrac3.mjs (Emscripten build)')
    ap.add_argument('--update', action='store_true')
    a = ap.parse_args()
    if a.wasm:
        tool = ['node', os.path.join(HERE, 'wasm_cli.mjs'), os.path.abspath(a.wasm)]
    else:
        exe = a.bin or os.path.join(os.path.dirname(HERE), 'build', 'atrac3')
        if not os.path.exists(exe) and os.path.exists(exe + '.exe'):
            exe += '.exe'
        tool = [os.path.abspath(exe)]

    tmp = tempfile.mkdtemp()
    try:
        gen_vectors.main_into(tmp)
        got, split_bad = {}, []
        for name, gen, n in gen_vectors.VECTORS:
            wav = os.path.join(tmp, name + '.wav')
            ch = 1 if name.endswith('_mono') else 2
            for kbps in MODES[ch]:
                at3 = os.path.join(tmp, '%s.%d.at3' % (name, kbps))
                dec = os.path.join(tmp, '%s.%d.wav' % (name, kbps))
                subprocess.run(tool + ['-e', '-br', str(kbps), wav, at3], check=True, stderr=subprocess.DEVNULL)
                subprocess.run(tool + ['-d', at3, dec], check=True)
                got[os.path.basename(at3)] = sha(at3)
                got[os.path.basename(dec)] = sha(dec)
                if a.wasm:
                    subprocess.run(tool + ['-e', '-split', '-br', str(kbps), wav, at3 + '.split'], check=True)
                    if sha(at3 + '.split') != got[os.path.basename(at3)]:
                        split_bad.append(os.path.basename(at3))
        riff_fails = riff_checks(tool, tmp, got)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if a.update:
        with open(EXPECTED, 'w', newline='\n') as f:
            for k in sorted(got):
                f.write('%s  %s\n' % (got[k], k))
        print('wrote %d hashes to %s' % (len(got), EXPECTED))
        return
    want = {}
    with open(EXPECTED) as f:
        for line in f:
            h, k = line.split()
            want[k] = h
    bad = [k for k in sorted(want) if got.get(k) != want[k]] + sorted(set(got) - set(want))
    for k in bad:
        print('FAIL', k, got.get(k, 'missing'))
    print('%d/%d outputs match' % (len(want) - len(bad), len(want)))
    for k in riff_fails:
        print('FAIL', k)
    print('RIFF chunk-size checks: %s' % ('%d failed' % len(riff_fails) if riff_fails else 'ok'))
    if a.wasm:
        for k in split_bad:
            print('FAIL streaming', k)
        n_split = sum(len(MODES[1 if v[0].endswith('_mono') else 2]) for v in gen_vectors.VECTORS)
        print('streaming encoder in pieces: %d/%d match' % (n_split - len(split_bad), n_split))
    sys.exit(1 if bad or riff_fails or split_bad else 0)


if __name__ == '__main__':
    main()
