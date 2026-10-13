"""Verify spec/reference/atrac3_encoder_lp4.py against Sony's emulated ATRAC3 'algorithm 0' encoder.

For every WAV and bit rate the Sony tool (psp_at3tool.exe under tools/at3emu.py) encodes the file
while the PCM blocks it hands to the library (a3e0_encode / a3e0_flush) are captured; the model
is driven with exactly the same calls and the output frames are compared byte for byte.  The
model's own WAV front end (encode_wav, which re-implements the tool's framing and container) is
also compared against the complete .at3 file.

Usage: python tools/verify_lp4.py [wav ...]            (default: signals from tools/gen_signals.py)
       options: --rates 66,105   --extra (also the signals of tools/gen_lp4_signals.py)
"""
import sys, os, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'spec', 'reference'))
import numpy as np
from unicorn.x86_const import UC_X86_REG_ESP
from at3emu import Emu
import atrac3_encoder_lp4 as M

A3E0_ENCODE, A3E0_FLUSH = 0x436ce0, 0x439f30
RET_ENCODE, RET_FLUSH = 0x404e2d, 0x4051bc      # instructions after the two call sites


def sony_encode(wav, kbps, out):
    """Run the tool; return (library calls, library output frames, .at3 file bytes).
    A call is the int16 PCM block given to a3e0_encode, or None for an a3e0_flush call."""
    calls, outs, pending = [], [], []
    e = Emu(quiet=True)

    def collect(emu):                     # runs right after a3e0_encode / a3e0_flush returned
        buf, pn = pending.pop()
        n = emu.r32(pn)
        if n:
            outs.append(emu.read(buf, n))

    def on_encode(emu):
        sp = emu.mu.reg_read(UC_X86_REG_ESP)
        pcm, n = emu.r32(sp + 8), emu.r32(emu.r32(sp + 16))
        calls.append(np.frombuffer(emu.read(pcm, n), '<i2').reshape(-1, 2).copy() if n else
                     np.zeros((0, 2), np.int16))
        pending.append((emu.r32(sp + 12), emu.r32(sp + 20)))

    def on_flush(emu):
        sp = emu.mu.reg_read(UC_X86_REG_ESP)
        calls.append(None)
        pending.append((emu.r32(sp + 8), emu.r32(sp + 12)))
    e.hook(A3E0_ENCODE, on_encode)
    e.hook(A3E0_FLUSH, on_flush)
    e.hook(RET_ENCODE, collect)
    e.hook(RET_FLUSH, collect)
    e.run_main(['-e', wav, out, '-br', str(kbps)])
    return calls, outs, open(out, 'rb').read()


def frames_of(at3):
    pos = 12
    while pos + 8 <= len(at3):
        cid, size = at3[pos:pos + 4], int.from_bytes(at3[pos + 4:pos + 8], 'little')
        if cid == b'data':
            return at3[pos + 8:pos + 8 + size]
        pos += 8 + size + (size & 1)
    raise ValueError('no data chunk')


def has_loop(wav):
    d = open(wav, 'rb').read()
    pos = 12
    while pos + 8 <= len(d):
        cid, size = d[pos:pos + 4], int.from_bytes(d[pos + 4:pos + 8], 'little')
        if cid == b'smpl':
            return True
        pos += 8 + size + (size & 1)
    return False


def model_encode(calls, kbps):
    enc = M.Encoder(kbps)
    frames = []
    for pcm in calls:
        if pcm is None:
            f = enc._call(0, None)          # one a3e0_flush call
        else:
            f = enc.encode(pcm)
        if f is not None:
            frames.append(f)
    return frames, enc.frame_bytes


def main(argv):
    rates = [66, 105]
    wavs = []
    extra = False
    i = 0
    while i < len(argv):
        if argv[i] == '--extra':
            extra = True
            i += 1
        elif argv[i] == '--rates':
            rates = [int(r) for r in argv[i + 1].split(',')]
            i += 2
        else:
            wavs.append(argv[i])
            i += 1
    tmp = tempfile.mkdtemp(prefix='verify_lp4_')
    if not wavs or extra:
        import gen_signals, gen_lp4_signals
        sets = ([gen_signals.signals()] if not wavs else []) + ([gen_lp4_signals.signals()] if extra else [])
        for sigs in sets:
            for name, sig in sigs.items():
                p = os.path.join(tmp, name + '.wav')
                gen_signals.save(p, sig)
                wavs.append(p)
    print('%-28s %5s %7s %10s %9s %s' % ('signal', 'kbps', 'frames', 'identical', 'percent', 'whole .at3 file'))
    all_ok = True
    for wav in wavs:
        for kbps in rates:
            t0 = time.time()
            sony_at3 = os.path.join(tmp, 'sony.at3')
            calls, sony, data = sony_encode(wav, kbps, sony_at3)
            mine, fb = model_encode(calls, kbps)
            ns = len(sony)
            same = sum(1 for k in range(min(ns, len(mine))) if mine[k] == sony[k])
            mine_at3 = os.path.join(tmp, 'mine.at3')
            M.encode_wav(wav, mine_at3, kbps)
            if has_loop(wav):
                file_res = 'n/a (smpl loop: tool loop mode not modelled)'
                file_ok = True
            else:
                file_ok = open(mine_at3, 'rb').read() == data
                file_res = 'identical' if file_ok else 'DIFFERS'
            ok = same == ns == len(mine) and file_ok
            all_ok &= ok
            print('%-28s %5d %7d %10d %8.2f%% %s   (%.1fs)' % (os.path.basename(wav), kbps, ns, same,
                  100.0 * same / max(ns, 1), file_res, time.time() - t0))
    print('ALL IDENTICAL' if all_ok else 'MISMATCHES FOUND')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
