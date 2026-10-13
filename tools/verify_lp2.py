"""Verify spec/reference/atrac3_encoder_lp2.py against the emulated Sony encoder (algorithm 1).

For every test signal the WAV is encoded twice -- by psp_at3tool.exe under emulation
(`-e -br 132`, or `-br 66` for the mono variant) and by the Python reference encoder -- and
the two .at3 files are compared frame by frame.

  .venv/bin/python tools/verify_lp2.py              # gen_signals.py set, stereo LP2
  .venv/bin/python tools/verify_lp2.py --mono       # ... plus 66 kbps mono (left channel)
  .venv/bin/python tools/verify_lp2.py --extra      # ... plus edge-case signals (clipping, DC, ...)
  .venv/bin/python tools/verify_lp2.py --music      # ... plus the PSP SDK music WAVs
  .venv/bin/python tools/verify_lp2.py --tables     # check lp2_tables against the binary only
  .venv/bin/python tools/verify_lp2.py file.wav ... # specific files
"""
import sys, os, tempfile, wave, argparse, time
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'spec', 'reference'))
import numpy as np

SDK_WAVE = os.path.join(HERE, '..', 'private', 'psp_sdk_660', 'usr', 'local', 'psp', 'devkit',
                        'data', 'sound', 'wave')


def save_wav(path, x):
    x = np.clip(np.round(x), -32768, 32767).astype('<i2')
    if x.ndim == 1:
        x = x[:, None]
    with wave.open(path, 'wb') as w:
        w.setnchannels(x.shape[1]); w.setsampwidth(2); w.setframerate(44100)
        w.writeframes(x.tobytes())


def extra_signals(d):
    """Edge cases: full-scale clipping, near-silence, clicks, DC, sweep, odd length."""
    rng = np.random.default_rng(99)
    n = 88200; t = np.arange(n) / 44100.0
    out = {}
    out['x_clip'] = np.stack([32767 * np.sign(np.sin(2 * np.pi * 997 * t)),
                              np.clip(rng.standard_normal(n) * 30000, -32768, 32767)], 1)
    out['x_quiet'] = np.stack([rng.integers(-1, 2, n), 3 * np.sin(2 * np.pi * 440 * t)], 1)
    c = np.zeros(n); c[::3001] = 32767
    out['x_clicks'] = np.stack([c, -c], 1)
    out['x_dc'] = np.stack([np.where((t * 4).astype(int) % 2 == 0, 20000.0, -20000.0),
                            np.full(n, 12345.0)], 1)
    ph = 2 * np.pi * 20.0 * 2.0 / np.log(1000.0) * (np.exp(t / 2.0 * np.log(1000.0)) - 1)
    sw = 29200 * np.sin(ph)
    out['x_sweep'] = np.stack([sw, -np.roll(sw, 100)], 1)
    seg = np.minimum((t * 3).astype(int), 2)  # weak isolated tones (tonal class 1 and 2)
    x = np.choose(seg, [40.0, 150.0, 400.0]) * np.sin(2 * np.pi * np.choose(seg, [1000.0, 3300.0, 9000.0]) * t)
    out['x_weaktone'] = np.stack([x, 0.5 * x], 1)
    m = 1024 * 20 - 955 + 950          # final block >= 886 samples: three flush calls
    out['x_len'] = np.stack([np.sin(np.arange(m) * 0.05) * 9000, np.cos(np.arange(m) * 0.031) * 9000], 1)
    paths = []
    for k, v in out.items():
        p = os.path.join(d, k + '.wav'); save_wav(p, v); paths.append(p)
    return paths


def frames_of(path):
    data = open(path, 'rb').read()
    i = data.find(b'data')
    hdr, body = data[:i + 8], data[i + 8:]
    align = int.from_bytes(data[32:34], 'little')
    return hdr, [body[k:k + align] for k in range(0, len(body), align)]


def run_one(job):
    wav, outdir, br = job
    from at3emu import Emu
    import atrac3_encoder_lp2 as E
    base = os.path.splitext(os.path.basename(wav))[0] + '_%d' % br
    sony, mine = os.path.join(outdir, base + '_sony.at3'), os.path.join(outdir, base + '_py.at3')
    t0 = time.time()
    Emu(quiet=True).run_main(['-e', wav, sony, '-br', str(br)])
    t1 = time.time()
    E.encode_file(wav, mine)
    t2 = time.time()
    if not os.path.exists(sony) or os.path.getsize(sony) == 0:
        return base, None
    hs, fs = frames_of(sony)
    hm, fm = frames_of(mine)
    same = sum(a == b for a, b in zip(fs, fm))
    first_bad = next((k for k, (a, b) in enumerate(zip(fs, fm)) if a != b), None)
    return base, dict(frames=len(fs), frames_py=len(fm), same=same, header=hs == hm,
                      first_bad=first_bad, t_emu=t1 - t0, t_py=t2 - t1)


def check_tables():
    from peutil import f32s, i32s
    import lp2_tables as T
    checks = [('SF', T.SF, f32s(0x48eae8, 64)), ('T1', T.T1, f32s(0x48ec68, 13)),
              ('QMF_WIN', T.QMF_WIN, f32s(0x490fc0, 48)), ('MDCT_WIN', T.MDCT_WIN, f32s(0x4907c0, 512)),
              ('PRE_COS', T.PRE_COS, f32s(0x48f7c0, 128)), ('PRE_SIN', T.PRE_SIN, f32s(0x48f9c0, 128)),
              ('FFT_COS', T.FFT_COS, f32s(0x4903c0, 64)), ('FFT_SIN', T.FFT_SIN, f32s(0x4904c0, 64)),
              ('POST1', T.POST1, f32s(0x48fbc0, 128)), ('POST2', T.POST2, f32s(0x48fdc0, 128)),
              ('POST3', T.POST3, f32s(0x48ffc0, 128)), ('POST4', T.POST4, f32s(0x4901c0, 128)),
              ('T2', T.T2, f32s(0x48ec9c, 8))]
    ok = True
    for name, mine, ref in checks:
        a = np.asarray(mine, np.float32).view(np.uint32)
        b = np.asarray(ref, np.float32).view(np.uint32)
        bad = int(np.sum(a != b)); ok &= bad == 0
        print('  %-9s %4d entries, %d differ (bitwise)' % (name, len(b), bad))
    ints = [('SB_START', T.SB_START, i32s(0x48e980, 33)), ('MAXQ', T.MAXQ, i32s(0x48ea24, 8)),
            ('QBITS', T.QBITS, i32s(0x48ea04, 8)), ('TONE_NVALS', T.TONE_NVALS, i32s(0x48ea88, 8)),
            ('BITREV7', T.BITREV7, i32s(0x4905c0, 128)), ('GAIN_EXP', T.GAIN_LEVEL_EXP, i32s(0x48eaa8, 16)),
            ('MASK_N', T.MASK_SLOPE_NORMAL, i32s(0x48ec28, 8)), ('MASK_T', T.MASK_SLOPE_TONAL, i32s(0x48ec48, 8)),
            ('EST_CLC', T.EST_BITS_CLC, i32s(0x48ebe8 + 32, 8)), ('VPC', T.VALS_PER_CODE, i32s(0x48ed40, 8)),
            ('CELL_BAND', T.CELL_BAND, i32s(0x48ea48, 16))]
    for name, mine, ref in ints:
        same = list(mine) == list(ref); ok &= same
        print('  %-9s %s' % (name, 'identical' if same else 'DIFFERENT'))
    # code tables, as laid out in the binary (0x48edc0 VLC, 0x48f2c0 CLC)
    from peutil import u32s
    for mode, base in ((0, 0x48edc0), (1, 0x48f2c0)):
        k = 0
        for sel in range(1, 8):
            n = (1 << T.QBITS[sel]) ** T.VALS_PER_CODE[sel]
            raw = u32s(base + 8 * k, 2 * n); k += n
            ref = [(raw[2 * i], raw[2 * i + 1]) for i in range(n)]
            mine = [c if c[1] else (0, 0) for c in T.ENC_CODES[(mode, sel)]]
            ref = [c if c[1] else (0, 0) for c in ref]
            if mine != ref:
                ok = False; print('  code table (%d,%d) DIFFERENT' % (mode, sel))
    print('  code tables %s' % ('identical' if ok else 'see above'))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    ap.add_argument('--mono', action='store_true', help='also test 66 kbps mono (left channel)')
    ap.add_argument('--extra', action='store_true', help='add edge-case signals')
    ap.add_argument('--music', action='store_true', help='add PSP SDK music WAVs')
    ap.add_argument('--tables', action='store_true', help='only check tables against the binary')
    ap.add_argument('-j', type=int, default=max(1, min(8, mp.cpu_count() - 1)))
    ap.add_argument('--keep', help='directory for the WAV/AT3 files (default: temporary)')
    a = ap.parse_args()
    print('tables vs psp_at3tool.exe:')
    tables_ok = check_tables()
    if a.tables:
        return 0 if tables_ok else 1
    d = a.keep or tempfile.mkdtemp(prefix='verify_lp2_')
    os.makedirs(d, exist_ok=True)
    wavs = list(a.files)
    if not wavs:
        import gen_signals
        for k, v in gen_signals.signals().items():
            p = os.path.join(d, k + '.wav'); gen_signals.save(p, v); wavs.append(p)
    if a.extra:
        wavs += extra_signals(d)
    if a.music:
        # the SDK WAVs carry `smpl` loop points, which make at3tool use its PSP loop framing
        # (spec ATRAC3.md 9.2, not implemented in the reference encoder): re-save the PCM only
        for name in ('bgm01_1.wav', 'bgm01_2.wav', 'overload_short.wav'):
            p = os.path.join(SDK_WAVE, name)
            if os.path.exists(p):
                with wave.open(p) as f:
                    x = np.frombuffer(f.readframes(f.getnframes()), '<i2').reshape(-1, f.getnchannels())
                q = os.path.join(d, 'sdk_' + name); save_wav(q, x); wavs.append(q)
    jobs = [(w, d, 132) for w in wavs]
    if a.mono:
        for w in list(wavs):
            with wave.open(w) as f:
                x = np.frombuffer(f.readframes(f.getnframes()), '<i2').reshape(-1, f.getnchannels())
            p = os.path.join(d, os.path.splitext(os.path.basename(w))[0] + '_mono.wav')
            save_wav(p, x[:, 0]); jobs.append((p, d, 66))
    print('\n%-28s %7s %9s %8s %6s  %s' % ('signal', 'frames', 'identical', '%', 'header', 'time emu/py'))
    tot = same = 0
    with mp.Pool(a.j) as pool:
        for base, r in pool.imap(run_one, jobs):
            if r is None:
                print('%-28s  (rejected by at3tool, e.g. too short)' % base); continue
            tot += r['frames']; same += r['same']
            extra = '' if r['frames'] == r['frames_py'] else '  frame count %d vs %d' % (r['frames'], r['frames_py'])
            if r['first_bad'] is not None:
                extra += '  first difference in frame %d' % r['first_bad']
            print('%-28s %7d %9d %7.2f%% %6s  %5.1fs/%4.1fs%s' % (base, r['frames'], r['same'],
                  100.0 * r['same'] / max(r['frames'], 1), 'same' if r['header'] else 'DIFF',
                  r['t_emu'], r['t_py'], extra))
    print('\ntotal: %d / %d frames identical (%.3f%%)' % (same, tot, 100.0 * same / max(tot, 1)))
    print('files in', d)
    return 0 if same == tot and tables_ok else 1


if __name__ == '__main__':
    sys.exit(main())
