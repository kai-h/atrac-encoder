"""Check MZ-RH1 rips of the gapless test set (tools/make_gapless_kit.py) against the frames the
atrac-encoder Web MiniDisc Pro plugin sent.

Give the ripped tracks of one upload in track order (ATRAC3 in WAV, or unencrypted OMA).  The mode
comes from the frame size; the set, and whether the upload was gapless, from the frames themselves.
The source is re-encoded with the command-line tool to reproduce exactly what was sent, and each
track is compared frame by frame: on a disc that plays gaplessly every track is identical.
Usage:  python tools/gapless_check.py RIP1 RIP2 ... [--kit testkit/gapless] [--bin build/atrac3]
"""
import os, sys, json, struct, subprocess, tempfile, wave, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kit_align import _riff_chunks, _oma
from make_gapless_kit import frames_sent

MODES = {384: (132, 'LP2'), 192: (66, 'LP4')}


def read_rip(path):
    """Return (frame bytes, ATRAC3 data) of a ripped track."""
    d = open(path, 'rb').read()
    if d[:4] == b'RIFF' and d[8:12] == b'WAVE':
        ch = dict(_riff_chunks(d))
        tag, _, _, _, align = struct.unpack('<HHIIH', ch[b'fmt '][:14])
        if tag != 0x270:
            sys.exit('%s: not ATRAC3 (WAV format 0x%x)' % (path, tag))
        return align, ch[b'data']
    if d[:3] in (b'ea3', b'EA3'):
        _, align, _, body = _oma(d)
        return align, body
    sys.exit('%s: not an ATRAC3 WAV or OMA file' % path)


def encode(tool, pcm, kbps, tmp):
    """Frames at3tool (and so the command-line tool) writes for 16-bit stereo PCM bytes."""
    src, out = os.path.join(tmp, 'in.wav'), os.path.join(tmp, 'out.at3')
    with wave.open(src, 'wb') as w:
        w.setnchannels(2); w.setsampwidth(2); w.setframerate(44100); w.writeframes(pcm)
    subprocess.run([tool, '-e', '-br', str(kbps), src, out], check=True, stderr=subprocess.DEVNULL)
    return open(out, 'rb').read()[80:]


def expected(kit, tracks, kbps, align, tool, tmp):
    """{'gapless': [...], 'separate': [...]}: the frames sent for each track of either upload."""
    pcm = []
    for t in tracks:
        with wave.open(os.path.join(kit, t['file'])) as w:
            pcm.append(w.readframes(w.getnframes()))
    joined, counts = encode(tool, b''.join(pcm), kbps, tmp), frames_sent([t['samples'] for t in tracks])
    assert len(joined) == sum(counts) * align, 'frame count model disagrees with the encoder'
    gapless, pos = [], 0
    for n in counts:
        gapless.append(joined[pos:pos + n * align])
        pos += n * align
    return dict(gapless=gapless, separate=[encode(tool, p, kbps, tmp) for p in pcm])


def frames(data, align):
    return [data[i:i + align] for i in range(0, len(data) - align + 1, align)]


def leading_match(got, want, align):
    g, w = frames(got, align), frames(want, align)
    return next((i for i in range(min(len(g), len(w))) if g[i] != w[i]), min(len(g), len(w)))


def describe(extra, sent):
    if all(not any(f) for f in extra):
        return 'all zero bytes'
    if all(f in sent for f in extra):
        return 'copies of frames that were sent'
    return 'new content, starting %s' % extra[0][:8].hex(' ')


def compare(got, want, align, sent):
    if got == want:
        return 'identical'
    g, w = frames(got, align), frames(want, align)
    same = leading_match(got, want, align)
    msg = []
    if len(got) % align:
        msg.append('%d stray bytes after the last whole frame' % (len(got) % align))
    if same == len(w) and len(g) > len(w):
        msg.append('all %d frames identical, then %d extra (%s)' % (len(w), len(g) - len(w), describe(g[len(w):], sent)))
    elif same == len(g) and len(g) < len(w):
        msg.append('first %d frames identical, the last %d are missing' % (len(g), len(w) - len(g)))
    elif g[1:1 + len(w)] == w:
        msg.append('1 extra frame at the start (%s), then all %d identical' % (describe(g[:1], sent), len(w)))
    else:
        msg.append('first %d frames identical, then they differ (%d ripped, %d sent)' % (same, len(g), len(w)))
    return '; '.join(msg)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('rips', nargs='+')
    ap.add_argument('--kit', default=os.path.join(HERE, '..', 'testkit', 'gapless'))
    ap.add_argument('--bin', default=os.path.join(HERE, '..', 'build', 'atrac3'))
    a = ap.parse_args(argv)
    manifest = json.load(open(os.path.join(a.kit, 'manifest.json')))
    rips = [read_rip(p) for p in a.rips]
    align = rips[0][0]
    if align not in MODES or any(r[0] != align for r in rips):
        sys.exit('frame sizes %s: expected all 384 (LP2) or all 192 (LP4)' % sorted({r[0] for r in rips}))
    kbps, mode = MODES[align]

    with tempfile.TemporaryDirectory() as tmp:
        cands = []
        for name, s in manifest['sets'].items():
            for kind, want in expected(a.kit, s['tracks'], kbps, align, a.bin, tmp).items():
                score = sum(leading_match(r[1], w, align) for r, w in zip(rips, want))
                cands.append((score, name, kind, s['tracks'], want))
    score, name, kind, tracks, want = max(cands, key=lambda c: c[0])
    if score == 0:
        sys.exit('the rips match no upload of the test set in %s' % mode)

    print('%s (%d kbps, %d-byte frames): %s set, %s upload, %d of %d tracks ripped' %
          (mode, kbps, align, name, kind, len(rips), len(tracks)))
    sent = {f for w in want for f in frames(w, align)}
    ok = True
    for k, (path, (_, data)) in enumerate(zip(a.rips, rips)):
        w = want[k] if k < len(want) else b''
        result = compare(data, w, align, sent)
        ok &= result == 'identical'
        mark = '  (mark %.1f ms early)' % tracks[k]['mark_early_ms'] if kind == 'gapless' and k else ''
        print('  %d  %-28s sent %3d, ripped %3d frames: %s%s' %
              (k + 1, os.path.basename(path)[:28], len(w) // align, len(data) // align, result, mark))
    if len(rips) < len(tracks):
        print('  (tracks %d-%d not given)' % (len(rips) + 1, len(tracks)))
    if kind == 'separate':
        print('Each track was encoded as a stream of its own, so expect a gap of 25 to 48 ms at every join.')
    elif ok:
        print('The disc holds exactly the frames that were sent: the joins are continuous in the data.')
    else:
        print('The disc does not hold exactly what was sent; see the differences above.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
