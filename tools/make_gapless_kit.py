"""Generate the gapless upload test set (testkit/gapless/): two continuous signals, a steady tone and
a log sweep, each cut into five tracks off the 1024-sample ATRAC3 frame grid.  Uploaded with Web
MiniDisc Pro's Gapless option, the track changes must be inaudible; tools/gapless_check.py compares
an MZ-RH1 rip with the frames the atrac-encoder plugin sent.

Deterministic: re-running produces identical files.   Usage: python tools/make_gapless_kit.py [outdir]
"""
import os, sys, json, wave
import numpy as np

SR = 44100
PEAK_DB = -12.0
FADE = 0.02            # seconds, at the very start and end only, so a click can only come from a join

# Track lengths in samples.  Each set mixes odd and even per-track frame counts, because an LP4
# sound group holds two frames and a track with an odd count may be padded on the disc.
SETS = {
    'tone': dict(desc='440 Hz sine in both channels, continuous across all five tracks',
                 lengths=[369117, 313591, 420863, 301207, 353825]),
    'sweep': dict(desc='log sweep from 20 Hz to 20 kHz in both channels, continuous across all five tracks',
                  lengths=[341333, 390001, 289911, 407777, 332307]),
}


def frames_sent(lengths):
    """Frames per track that the plugin returns for a gapless upload.  After S samples the
    streaming encoder has returned floor((955 + S) / 1024) - 2 frames; finish() then completes
    at3tool's frame count for the whole stream (spec/ATRAC3.md section 6.4)."""
    t = sum(lengths) + 955
    whole = -(-t // 1024) + (0 if 0 < t % 1024 < 886 else 1)
    ends = np.cumsum(lengths)
    emitted = [max(0, (955 + int(s)) // 1024 - 2) for s in ends[:-1]] + [whole]
    return [int(n) for n in np.diff([0] + emitted)]


def signal(name, n):
    t = np.arange(n) / SR
    if name == 'tone':
        x = np.sin(2 * np.pi * 440.0 * t)
    else:
        f0, f1, dur = 20.0, 20000.0, n / SR
        k = np.log(f1 / f0)
        x = np.sin(2 * np.pi * f0 * dur / k * (np.exp(t / dur * k) - 1))
    m = int(FADE * SR)
    w = 0.5 - 0.5 * np.cos(np.pi * np.arange(m) / m)
    x[:m] *= w
    x[-m:] *= w[::-1]
    return x * 10 ** (PEAK_DB / 20)


def build(outdir):
    os.makedirs(outdir, exist_ok=True)
    manifest = dict(sample_rate=SR, peak_dbfs=PEAK_DB, sets={})
    for name, s in SETS.items():
        lengths = s['lengths']
        pcm = np.round(signal(name, sum(lengths)) * 32768).astype('<i2')
        counts = frames_sent(lengths)
        tracks, start, sent = [], 0, 0
        for k, (n, f) in enumerate(zip(lengths, counts)):
            path = os.path.join(outdir, '%s-%d.wav' % (name, k + 1))
            with wave.open(path, 'wb') as w:
                w.setnchannels(2); w.setsampwidth(2); w.setframerate(SR)
                w.writeframes(np.repeat(pcm[start:start + n], 2).tobytes())
            # how far before its audio this track's first frame starts (decoded i = source i - 1093)
            early = (start + 1093 - 1024 * sent) / SR * 1000 if k else 0.0
            tracks.append(dict(file=os.path.basename(path), start=start, samples=n, frames=f,
                               mark_early_ms=round(early, 1)))
            print('%-11s %5.2f s  %3d frames (%s)%s' % (os.path.basename(path), n / SR, f,
                  'odd' if f % 2 else 'even', '  mark %.1f ms early' % early if k else ''))
            start += n
            sent += f
        manifest['sets'][name] = dict(desc=s['desc'], tracks=tracks)
    with open(os.path.join(outdir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1)


if __name__ == '__main__':
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'testkit', 'gapless'))
