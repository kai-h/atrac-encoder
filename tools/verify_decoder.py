"""Verify spec/reference/atrac3_decoder.py against the emulated Sony decoder, frame by frame.
Also report which bitstream features the stream exercises."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'spec', 'reference'))
import numpy as np
from unicorn.x86_const import UC_X86_REG_ESP
from at3emu import Emu
import atrac3_decoder as R

def sony_frames(path):
    post, cur = [], {}
    e = Emu(quiet=True)
    def on_qmf(emu): cur['p'] = emu.r32(emu.mu.reg_read(UC_X86_REG_ESP) + 4)
    def after(emu): post.append(np.frombuffer(emu.read(cur['p'], 4096), dtype='<f4').astype(np.float64))
    e.hook(0x43abd0, on_qmf); e.hook(0x43baa9, after)
    e.run_main(['-d', '-repeat', '1', path, os.devnull])
    return post

def stats_hook():
    st = collections.Counter()
    orig = R.decode_unit
    def wrapped(br, js_second):
        nb, gain, spec = orig(br, js_second)
        st['units'] += 1
        st['gain_points'] += sum(len(g) for g in gain)
        st['bands_%d' % nb] += 1
        return nb, gain, spec
    R.decode_unit = wrapped
    return st

if __name__ == '__main__':
    for path in sys.argv[1:]:
        post = sony_frames(path)
        R.STATS.clear()
        ch, align, js, body = R.read_at3(path)
        dec = R.Decoder(ch, align, js)
        nframes = len(body) // align
        worst, f = -999.0, 0
        mats, wts = collections.Counter(), collections.Counter()
        for f, i in enumerate(range(0, len(body) - align + 1, align)):
            out = dec.decode_frame(body[i:i + align])
            if js:
                mats.update(dec.js_cur[1]); wts[dec.js_cur[0]] += 1
            per = len(post) // nframes          # Sony decodes mono as duplicated stereo
            for c in range(ch):
                S = post[per * f + c]; e2 = np.sum((S - out[:, c]) ** 2)
                if np.sum(S ** 2) > 1.0:
                    worst = max(worst, 10 * np.log10(e2 / np.sum(S ** 2) + 1e-300))
        print('%-28s frames=%4d worst=%7.1f dB  %s  matrix=%s weights=%s' % (
            os.path.basename(path), f + 1, worst, dict(R.STATS), dict(mats), dict(wts)))
