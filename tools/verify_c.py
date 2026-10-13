"""Check the C implementation (build/atrac3) against the Python reference models.

For every input WAV and every mode it supports, the C encoder's .at3 must be byte-identical
to the reference model's.  The C decoder's output is compared with the reference decoder's
(after the same delay trimming); the two differ only in floating-point summation order.

Usage:  .venv/bin/python tools/verify_c.py [--bin build/atrac3] [--no-decode] wav-or-dir ...
"""
import argparse, glob, os, struct, subprocess, sys, tempfile, wave
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'spec', 'reference'))
import atrac3_encoder_lp2 as lp2       # noqa: E402
import atrac3_encoder_lp4 as lp4       # noqa: E402
import atrac3_decoder as dec          # noqa: E402

MODES = {2: [132, 105, 66], 1: [66, 52]}


def read_pcm(path):
    with wave.open(path, 'rb') as w:
        ch = w.getnchannels()
        return np.frombuffer(w.readframes(w.getnframes()), '<i2').reshape(-1, ch)


def ref_lp4_frames(pcm, kbps):
    """at3tool framing around the algorithm-0 model (atrac3_encoder_lp4.encode_wav)."""
    pcm = np.concatenate([np.zeros((955, 2), np.int16), pcm])
    enc = lp4.Encoder(kbps)
    frames = []
    for i in range(0, len(pcm), 1024):
        f = enc.encode(pcm[i:i + 1024])
        if f is not None:
            frames.append(f)
    frames += enc.flush()
    return frames[1:]


def reference(path, kbps):
    """The reference model's .at3 file for one input and mode."""
    pcm = read_pcm(path)
    n, ch = pcm.shape
    if kbps == 132 or (kbps == 66 and ch == 1):
        frames = lp2.encode_pcm(pcm)
        return lp2.at3_header(ch, n, sum(map(len, frames))) + b''.join(frames)
    if ch == 1:                       # 52 kbps mono: unit 0 of the 105 kbps stereo encoder
        frames = [f[:152] for f in ref_lp4_frames(np.repeat(pcm, 2, 1), 105)]
        align, js = 152, 0
    else:
        frames = ref_lp4_frames(pcm, kbps)
        align, js = len(frames[0]), int(kbps == 66)
    body = b''.join(frames)
    fmt = struct.pack('<HHIIHHHHIHHHH', 0x270, ch, 44100, int(44100 * align / 1024 + 0.5), align,
                      0, 14, 1, ch * 2048, js, js, 1, 0)
    fact = struct.pack('<III', n, 1024, 1024)
    return (b'RIFF' + struct.pack('<I', 4 + 8 + len(fmt) + 8 + len(fact) + 8 + len(body)) +
            b'WAVE' + b'fmt ' + struct.pack('<I', len(fmt)) + fmt +
            b'fact' + struct.pack('<I', len(fact)) + fact + b'data' + struct.pack('<I', len(body)) + body)


def ref_decode(at3_path, nsamples):
    ch, align, js, body = dec.read_at3(at3_path)
    d = dec.Decoder(ch, align, js)
    out = np.concatenate([d.decode_frame(body[i:i + align]) for i in range(0, len(body) - align + 1, align)])
    out = np.round(np.clip(out, -32767, 32767)).astype(np.int16)
    return out[1093:1093 + nsamples]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', default=os.path.join(ROOT, 'build', 'atrac3'))
    ap.add_argument('--no-decode', action='store_true')
    ap.add_argument('inputs', nargs='+')
    a = ap.parse_args()
    wavs = []
    for p in a.inputs:
        wavs += sorted(glob.glob(os.path.join(p, '*.wav'))) if os.path.isdir(p) else [p]
    tmp = tempfile.mkdtemp()
    bad = 0
    frames_total = 0
    for w in wavs:
        pcm = read_pcm(w)
        for kbps in MODES[pcm.shape[1]]:
            out = os.path.join(tmp, 'c.at3')
            subprocess.run([a.bin, '-e', '-br', str(kbps), w, out], check=True, stderr=subprocess.DEVNULL)
            c = open(out, 'rb').read()
            r = reference(w, kbps)
            align = struct.unpack('<H', r[32:34])[0]
            nfr = (len(r) - 80) // align
            frames_total += nfr
            if c == r:
                status = 'identical (%d frames)' % nfr
            else:
                bad += 1
                k = next((i for i in range(min(len(c), len(r))) if c[i] != r[i]), min(len(c), len(r)))
                status = 'MISMATCH at byte %d (frame %d), sizes %d/%d' % (k, (k - 80) // align, len(c), len(r))
            line = '%-28s %3d kbps  %s' % (os.path.basename(w)[:28], kbps, status)
            if not a.no_decode and c == r:
                wav_out = os.path.join(tmp, 'c.wav')
                subprocess.run([a.bin, '-d', out, wav_out], check=True)
                cd = read_pcm(wav_out)
                rd = ref_decode(out, len(pcm))
                if cd.shape != rd.shape:
                    line += '  decode SHAPE %s vs %s' % (cd.shape, rd.shape)
                    bad += 1
                else:
                    diff = np.abs(cd.astype(int) - rd.astype(int))
                    line += '  decode max|d|=%d (%d samples differ)' % (diff.max(), (diff > 0).sum())
                    bad += diff.max() > 1
            print(line, flush=True)
    print('%d failures, %d frames compared' % (bad, frames_total))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
