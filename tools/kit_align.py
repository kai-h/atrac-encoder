"""Find test-kit tracks in a rip of a MiniDisc recording and measure the capture path.

Accepts 16-bit PCM WAV, ATRAC3 in WAV (format 0x270), ATRAC3 in OMA (unencrypted EA3), or anything
FFmpeg decodes (e.g. ATRAC1 .aea). Locates each track by its sync chirps (spacing identifies the
track), then reports offset, clock ratio, per-channel gain, polarity, channel mapping and SNR against
the source WAV.  Usage:  python tools/kit_align.py rip1 [rip2 ...] [--kit testkit] [--save outdir]
"""
import sys, os, json, wave, struct, subprocess, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spec', 'reference'))
sys.path.insert(0, HERE)
from make_testkit import SYNC_A, SYNC_B, SR


# ------------------------------------------------------------------------------------- decoding
def _decode_atrac3(ch, align, js, body):
    import atrac3_decoder as R
    dec = R.Decoder(ch, align, js)
    frames = [dec.decode_frame(body[k:k + align]) for k in range(0, len(body) - align + 1, align)]
    return np.concatenate(frames) / 32768.0


def _riff_chunks(d):
    pos = 12
    while pos + 8 <= len(d):
        cid, size = d[pos:pos + 4], struct.unpack('<I', d[pos + 4:pos + 8])[0]
        yield cid, d[pos + 8:pos + 8 + size]
        pos += 8 + size + (size & 1)


def _oma(d):
    pos = 0
    if d[:3] == b'ea3':                                   # ID3v2-style tag with an 'ea3' identifier
        size = (d[6] << 21) | (d[7] << 14) | (d[8] << 7) | d[9]
        pos = 10 + size
    h = d[pos:pos + 96]
    if h[:3] != b'EA3' or h[5] != 96:
        raise ValueError('no EA3 header')
    if struct.unpack('>H', h[6:8])[0] != 0xFFFF:
        raise ValueError('encrypted OMA file; export it decrypted')
    codec, params = h[32], int.from_bytes(h[33:36], 'big')
    if codec != 0:
        raise ValueError('OMA codec %d is not ATRAC3' % codec)
    return 2, (params & 0x3FF) * 8, (params >> 17) & 1, d[pos + 96:]


def load(path):
    """Return float PCM in [-1, 1), shape (n, channels)."""
    d = open(path, 'rb').read()
    if d[:4] == b'RIFF' and d[8:12] == b'WAVE':
        ch = dict(_riff_chunks(d))
        tag, nch, sr, _, align = struct.unpack('<HHIIH', ch[b'fmt '][:14])
        if tag == 1:
            return np.frombuffer(ch[b'data'], '<i2').reshape(-1, nch) / 32768.0
        if tag == 0x270:
            fmt = ch[b'fmt ']
            js = struct.unpack('<H', fmt[24:26])[0] if len(fmt) >= 26 else 0
            return _decode_atrac3(nch, align, js, ch[b'data'])
    if d[:3] in (b'ea3', b'EA3'):
        return _decode_atrac3(*_oma(d))
    raw = subprocess.run(['ffmpeg', '-v', 'error', '-i', path, '-f', 's16le', '-acodec', 'pcm_s16le', '-'],
                         capture_output=True, check=True)
    probe = subprocess.run(['ffprobe', '-v', 'error', '-select_streams', 'a:0', '-show_entries',
                            'stream=channels', '-of', 'csv=p=0', path], capture_output=True, text=True)
    nch = int(probe.stdout.strip() or 2)
    return np.frombuffer(raw.stdout, '<i2').reshape(-1, nch) / 32768.0


# ------------------------------------------------------------------------------------- sync search
def _xcorr_norm(x, tmpl):
    n = len(x) + len(tmpl)
    nfft = 1 << (n - 1).bit_length()
    c = np.fft.irfft(np.fft.rfft(x, nfft) * np.conj(np.fft.rfft(tmpl, nfft)), nfft)[:len(x)]
    e = np.convolve(x * x, np.ones(len(tmpl)), 'full')[len(tmpl) - 1:len(tmpl) - 1 + len(x)]
    return c / np.sqrt(np.maximum(e, 1e-12) * np.dot(tmpl, tmpl))


def find_sync(mono, tmpl, thresh=0.6):
    c = _xcorr_norm(mono, tmpl)
    hits, k, w = [], 0, len(tmpl)
    idx = np.nonzero(np.abs(c) > thresh)[0]
    while k < len(idx):
        lo = idx[k]; grp = idx[(idx >= lo) & (idx < lo + w)]
        best = grp[np.argmax(np.abs(c[grp]))]
        if 0 < best < len(c) - 1:                          # parabolic sub-sample refinement
            y0, y1, y2 = np.abs(c[best - 1:best + 2]); den = y0 - 2 * y1 + y2
            frac = 0.5 * (y0 - y2) / den if den else 0.0
        else:
            frac = 0.0
        hits.append((best + frac, c[best]))
        k = np.searchsorted(idx, lo + w)
    return hits


# ------------------------------------------------------------------------------------- analysis
def _resample(x, a, sync_a, ratio, start, n):
    """Rip samples at source positions start..start+n-1, compensating offset and clock ratio."""
    pos = a + (np.arange(start, start + n) - sync_a) * ratio
    if pos[0] < 0 or pos[-1] > len(x) - 2:
        return None
    return np.stack([np.interp(pos, np.arange(len(x)), x[:, c]) for c in range(x.shape[1])], 1)


def _fit(y, refs):
    """Best of the candidate references for one rip channel: (index, gain, snr_db)."""
    best = None
    for k, s in enumerate(refs):
        g = np.dot(y, s) / max(np.dot(s, s), 1e-20)
        err = np.sum((y - g * s) ** 2)
        if best is None or err < best[3]:
            best = (k, g, 10 * np.log10(np.sum((g * s) ** 2) / max(err, 1e-30)), err)
    return best[:3]


def analyse(path, kit, manifest, save=None):
    x = load(path)
    mono = x.mean(axis=1)
    A, B = find_sync(mono, SYNC_A), find_sync(mono, SYNC_B)
    dist = {t['sync_b'] - t['sync_a']: t for t in manifest['tracks']}
    names = ['L', 'R', '(L+R)/2']
    results = []
    for a, ca in A:
        for b, cb in B:
            gap = b - a
            t = min(dist, key=lambda d: abs(d - gap))
            if abs(t - gap) > 0.002 * t:                    # > 2000 ppm: not this pair
                continue
            tr = dist[t]
            ratio = gap / t
            src = load(os.path.join(kit, tr['file']))
            p0, n = tr['probe_start'], tr['probe_samples']
            seg = _resample(x, a, tr['sync_a'], ratio, p0, n)
            if seg is None:
                continue
            ref = src[p0:p0 + n]
            r = dict(file=os.path.basename(path), track=tr['file'][:-4], offset=a - tr['sync_a'],
                     clock_ppm=(ratio - 1) * 1e6, sync_corr=(round(float(ca), 3), round(float(cb), 3)), ch=[])
            segs = {sg['name']: sg for sg in tr['segments']}
            if 'tone1k_-20' in segs:                        # calibration track: path checks
                sg = segs['tone1k_-20']; lo, hi = sg['start'] + 4410, sg['start'] + sg['length'] - 4410
                rs = ref[lo:hi]
                for c in range(seg.shape[1]):
                    k, g, snr = _fit(seg[lo:hi, c], [rs[:, 0], rs[:, 1], rs.mean(axis=1)])
                    r['ch'].append(dict(source=names[k], gain_db=20 * np.log10(abs(g)), polarity='+' if g > 0 else '-', snr_db=snr))
                def energy(name):
                    sg = segs[name]; return np.sum(seg[sg['start']:sg['start'] + sg['length']] ** 2, axis=0)
                eL, eR = energy('left_only_1k'), energy('right_only_2k')
                r['channel_map'] = ('L->ch0 R->ch1' if seg.shape[1] == 2 and eL[0] > 10 * eL[1] and eR[1] > 10 * eR[0] else
                                    'SWAPPED' if seg.shape[1] == 2 and eL[1] > 10 * eL[0] and eR[0] > 10 * eR[1] else
                                    'mono/mixed')
            else:                                           # other tracks: whole-probe waveform fit
                for c in range(seg.shape[1]):
                    k, g, snr = _fit(seg[:, c], [ref[:, 0], ref[:, 1], ref.mean(axis=1)])
                    r['ch'].append(dict(source=names[k], gain_db=20 * np.log10(abs(g)), polarity='+' if g > 0 else '-', snr_db=snr))
            results.append(r)
            if save:
                os.makedirs(save, exist_ok=True)
                np.save(os.path.join(save, '%s__%s.npy' % (os.path.basename(path), tr['file'][:-4])),
                        _resample(x, a, tr['sync_a'], ratio, 0, tr['samples'] - 2))
            break
    return results


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('rips', nargs='+')
    ap.add_argument('--kit', default=os.path.join(HERE, '..', 'testkit'))
    ap.add_argument('--save', help='directory for aligned per-track excerpts (.npy)')
    a = ap.parse_args(argv)
    manifest = json.load(open(os.path.join(a.kit, 'manifest.json')))
    print('Track 01 (calibration): gain/SNR on the 1 kHz tone, plus channel map. Other tracks: whole-probe waveform\n'
          'fit; its SNR reflects the codec (noise is not waveform-coded), not the capture path.\n')
    print('%-28s %-16s %12s %8s  %s' % ('rip', 'track', 'offset', 'ppm', 'per channel: source gain polarity SNR'))
    for path in a.rips:
        res = analyse(path, a.kit, manifest, a.save)
        if not res:
            print('%-28s no kit track found' % os.path.basename(path))
        for r in res:
            chs = '  '.join('%s %+.3f dB %s %5.1f dB' % (v['source'], v['gain_db'], v['polarity'], v['snr_db'])
                            for v in r['ch'])
            extra = '  [%s]' % r['channel_map'] if 'channel_map' in r else ''
            print('%-28s %-16s %12.2f %8.1f  %s%s' % (r['file'][:28], r['track'], r['offset'], r['clock_ppm'], chs, extra))


if __name__ == '__main__':
    main(sys.argv[1:])
