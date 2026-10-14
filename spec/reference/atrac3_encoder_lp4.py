"""Reference model of Sony's ATRAC3 'algorithm 0' encoder (psp_at3tool.exe, a3e0_*).

This is the encoder used for MiniDisc LP4 (66 kbps joint stereo, 192-byte frames) and for
105 kbps stereo (304-byte frames).  The model is written to reproduce Sony's output byte for
byte; see spec/ATRAC3-encoder-LP4.md for the description.

Usage:
    enc = Encoder(66)                       # or 105
    for block in pcm_blocks:                # int16 array [n, 2], n <= 1024
        frame = enc.encode(block)           # bytes, or None (first call primes the encoder)
    frames = enc.flush()                    # trailing frames

    python atrac3_encoder_lp4.py in.wav out.at3 [66|105]    (framing identical to at3tool -e)
"""
import sys, struct, wave
import numpy as np
from lp4_tables import (IS_JOINT, BITRATE, BLOCK_ALIGN, CHANNEL_BYTES0, CHANNEL_BYTES1,
                        CHANNEL_HERTZ0, CHANNEL_HERTZ1, SB_START)
from lp4_analysis import qmf_analysis, envelope, mdct, QMF_STATE_LEN
from lp4_gain import BandGain, gain_band
from lp4_stereo import JointStereo, split_bits
from lp4_alloc import allocate, AllocState, UnitInfo
from lp4_bitstream import pack_unit

F32 = np.float32
FRAME = 1024


def _bandwidth_subbands(hz, rate=44100):
    """Number of coded subbands for a bandwidth in Hz (0x434580)."""
    idx = (hz << 11) // rate
    i = 0
    while i < 0x1f and not idx <= SB_START[i + 1]:
        i += 1
    return i + 1, idx


class Unit:
    """Per-channel encoder state (one 0x2d70-byte block in Sony's state)."""
    def __init__(self, budget, js, nsub_max):
        self.budget = budget
        self.js = js
        self.nsub_max = nsub_max
        self.qmf_state = np.zeros(QMF_STATE_LEN, F32)
        self.prev_time = np.zeros((256, 4), F32)      # previous frame's band samples
        self.mdct_state = np.zeros((128, 4), F32)
        self.gain = [BandGain() for _ in range(4)]
        self.prev_count0 = 0


class Encoder:
    def __init__(self, kbps=66):
        row = BITRATE.index(kbps)
        self.frame_bytes = BLOCK_ALIGN[row]
        self.joint = bool(IS_JOINT[row])
        cb = [CHANNEL_BYTES0[row], CHANNEL_BYTES1[row]]
        hz = [CHANNEL_HERTZ0[row], CHANNEL_HERTZ1[row]]
        self.bytes_per_ch = sum(cb) // 2
        self.units = []
        limit = -1
        for c in range(2):
            js = self.joint and c == 1
            budget = cb[c] * 8 - (8 if js else 0) - 19
            nsub, idx = _bandwidth_subbands(hz[c])
            if js:
                js = self.units[0].budget              # Sony stores unit 0's budget as the flag
                limit = max(1, (idx + 0xff) >> 8)
            self.units.append(Unit(budget, js, nsub))
        self.stereo = JointStereo(limit) if self.joint else None
        self.counter = 0

    # ------------------------------------------------------------------ public API
    def encode(self, pcm):
        """Encode one block of up to 1024 stereo int16 samples (a3e0_encode).  A block shorter
        than 1024 marks the end of the input.  Returns a frame (bytes) or None."""
        pcm = np.asarray(pcm, dtype=np.int16).reshape(-1, 2)
        return self._call(len(pcm) * 4, pcm)

    def flush(self):
        """Emit the frames still pending after the last block (a3e0_flush)."""
        out = []
        while True:
            f = self._call(0, None)
            if f is None:
                return out
            out.append(f)

    def _call(self, nbytes, pcm):
        c = self.counter
        self.counter = c + 1
        if nbytes != 4096 and self.counter > 0:
            self.counter = -3 if nbytes < 0xdd8 else -4
        if self.counter == 0:
            return None
        x = np.zeros((FRAME, 2), F32)
        if nbytes:
            x[:len(pcm)] = pcm
        frame = self.encode_frame(x)
        return None if self.counter == 1 else frame

    # ------------------------------------------------------------------ one frame
    def encode_frame(self, x):
        """0x436ea0: analyse 1024 new samples per channel, emit one frame."""
        X = []
        for c, u in enumerate(self.units):
            Xc, u.qmf_state = qmf_analysis(x[:, c].copy(), u.qmf_state)
            X.append(Xc)
        js = self.stereo
        if js:
            js.analyse(X[0], X[1])
            for b in range(4):
                js.convert(b, X[0], X[1])
        specs = []
        for c, u in enumerate(self.units):
            env = envelope(X[c])
            prev = u.prev_time.copy()
            for b in (3, 2, 1, 0):
                n = gain_band(u.gain, b, env[b], prev, bool(u.js), js.dec_prev if js else None,
                              js.dec if js else None, u)
                u.gain[b].count = n
                u.gain[b].lev[n] = 4
            spec, u.mdct_state = mdct(np.concatenate([X[c], prev]), u.mdct_state,
                                      [u.gain[b].lev[0] for b in range(4)])
            u.prev_time = X[c]
            specs.append(spec.reshape(-1))
        if js:
            u0, u1 = self.units
            u0.budget = split_bits(specs[0], specs[1], u0.nsub_max, u1.nsub_max, self.bytes_per_ch,
                                   u1.js, js.dec_prev[0], u1.gain[0].count)
        buf = bytearray(self.frame_bytes + 8)
        offset = 0
        for c, u in enumerate(self.units):
            info = UnitInfo(u.budget, u.js, u.nsub_max, [g.count for g in u.gain], u.prev_count0,
                            js.w_prev if js else 7)
            if u.budget < 40:
                S = AllocState(specs[c])
                S.nbands = S.nsub = 1
                S.ngroups = 0
                u.gain[0].count = 0
            else:
                S = allocate(info, specs[c])
            hdr = (js.w_prev, js.dec_prev) if (js and u.js) else None
            end = pack_unit(buf, offset, S, u.gain, hdr)
            if js:
                if c == 0:
                    offset = (end + 7) >> 3
                    self.units[1].budget = (self.bytes_per_ch * 2 - offset) * 8 - 0x1b
                else:
                    buf[offset:self.frame_bytes] = buf[offset:self.frame_bytes][::-1]
            else:
                offset += self.bytes_per_ch
        return bytes(buf[:self.frame_bytes])


# ---------------------------------------------------------------------- CLI
def encode_wav(src, dst, kbps=66):
    """Encode a 44.1 kHz stereo 16-bit WAV exactly like `at3tool -e <src> <dst> -br <kbps>`.

    Framing of the tool (cli_encode_loop 0x403290): the first block handed to the library holds
    955 zero samples followed by the first 69 input samples (1024 minus the 69-sample codec
    delay), then blocks of 1024 samples follow, the last one possibly short; afterwards the
    encoder is flushed.  The library itself emits nothing for its first call and the tool
    additionally discards the first frame it receives, so the file starts with the frame
    produced by the third call.  The container is RIFF/WAVE with a 32-byte 'fmt ' chunk (tag
    0x270) and a 'fact' chunk (total input samples, 1024, 1024)."""
    with wave.open(src, 'rb') as w:
        assert w.getnchannels() == 2 and w.getsampwidth() == 2
        nsamples = w.getnframes()
        pcm = np.frombuffer(w.readframes(nsamples), '<i2').reshape(-1, 2)
    pcm = np.concatenate([np.zeros((955, 2), np.int16), pcm])
    enc = Encoder(kbps)
    frames = []
    for i in range(0, len(pcm), FRAME):
        f = enc.encode(pcm[i:i + FRAME])
        if f is not None:
            frames.append(f)
    frames += enc.flush()
    frames = frames[1:]
    body = b''.join(frames)
    ba = enc.frame_bytes
    js = 1 if enc.joint else 0
    fmt = struct.pack('<HHIIHHHHIHHHH', 0x270, 2, 44100, int(44100 * ba / 1024 + 0.5), ba, 0, 14, 1,
                      0x1000, js, js, 1, 0)
    fact = struct.pack('<III', nsamples, 1024, 1024)
    with open(dst, 'wb') as f:
        f.write(b'RIFF' + struct.pack('<I', 4 + 8 + len(fmt) + 8 + len(fact) + 8 + len(body)) + b'WAVE')
        f.write(b'fmt ' + struct.pack('<I', len(fmt)) + fmt)
        f.write(b'fact' + struct.pack('<I', len(fact)) + fact)
        f.write(b'data' + struct.pack('<I', len(body)) + body)
    return frames


if __name__ == '__main__':
    encode_wav(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 66)
