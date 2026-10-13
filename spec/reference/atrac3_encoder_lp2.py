"""ATRAC3 encoder, Sony "algorithm 1" (MiniDisc LP2, 132 kbps stereo; also 66 kbps mono).

Executable companion to spec/drafts/encoder_lp2.md.  Produces output byte-identical to
psp_at3tool.exe 3.0.0.0 (`at3tool -e -br 132 in.wav out.at3`, or `-br 66` with mono input).

Usage:  python atrac3_encoder_lp2.py in.wav out.at3
        (stereo input -> 132 kbps LP2, mono input -> 66 kbps mono)

Structure:
  LP2Encoder.encode_call   one library call
    lp2_analysis           QMF split, gain control, MDCT
    lp2_alloc              tonal components, bit allocation
    lp2_bitstream          sound-unit serialisation
  encode_file              the tool's framing, priming and flushing, RIFF container
"""
import sys, struct, wave
import numpy as np
from lp2_tables import FRAME_SAMPLES, UNIT_BYTES, UNIT_BITS, PRIME_ZEROS, SAMPLE_RATE
from lp2_analysis import ChannelAnalysis
from lp2_alloc import encode_unit
from lp2_bitstream import pack_unit

END_THRESHOLD = 0x376      # 886: a call whose block had fewer samples ends the stream two calls later


class LP2Encoder:
    """Library-level encoder: two independent channels, 1024 samples per call, 192 bytes each."""

    def __init__(self, channels=2):
        self.channels = channels
        self.ana = [ChannelAnalysis() for _ in range(channels)]
        self.counts = [FRAME_SAMPLES] * 3      # samples per call, newest first (init: 1024)
        self.calls = 0

    def encode_call(self, block, nsamples):
        """block: int16 array (nsamples, channels), nsamples <= 1024 (zero-padded internally).
        Returns (frame bytes, end flag).  The frame of call k describes the input of call k-1
        (one block of internal delay)."""
        self.counts = [nsamples] + self.counts[:2]      # [this call, k-1, k-2]
        end = self.counts[2] < END_THRESHOLD
        pcm = np.zeros((FRAME_SAMPLES, self.channels), np.float32)
        if nsamples:
            pcm[:nsamples] = np.asarray(block, np.int16)[:nsamples]
        out = b''
        for c in range(self.channels):
            spec_gc, spec_plain, gain, prev_gain = self.ana[c].process(pcm[:, c])
            unit = encode_unit(spec_gc, spec_plain, gain, prev_gain, UNIT_BITS)
            out += pack_unit(unit, UNIT_BYTES)
        self.calls += 1
        return out, end


def at3_header(channels, total_samples, data_bytes):
    """RIFF/WAVE header as written by at3tool for ATRAC3 without loop points (spec §3)."""
    align = UNIT_BYTES * channels
    fmt = struct.pack('<HHIIHHH', 0x270, channels, SAMPLE_RATE,
                      int(SAMPLE_RATE * align / 1024.0 + 0.5), align, 0, 14)
    fmt += struct.pack('<HIHHHH', 1, channels * 2048, 0, 0, 1, 0)
    fact = struct.pack('<III', total_samples, 1024, 1024)
    body = b'WAVE' + b'fmt ' + struct.pack('<I', len(fmt)) + fmt + \
           b'fact' + struct.pack('<I', len(fact)) + fact + b'data' + struct.pack('<I', data_bytes)
    return b'RIFF' + struct.pack('<I', len(body) + data_bytes) + body


def encode_pcm(pcm):
    """pcm: int16 array (N, channels), channels 1 (-> 66 kbps mono) or 2 (-> 132 kbps LP2).
    Returns the list of output frames exactly as at3tool writes them."""
    pcm = np.asarray(pcm, np.int16)
    n, ch = pcm.shape
    # Mono: the tool duplicates the signal into both channels of the stereo encoder and keeps
    # unit 0.  Channels are encoded independently, so encoding channel 0 alone is equivalent.
    enc = LP2Encoder(1 if ch == 1 else 2)
    stream = np.concatenate([np.zeros((PRIME_ZEROS, ch), np.int16), pcm])
    frames = []
    pos = 0
    while pos < len(stream):                       # data calls
        blk = stream[pos:pos + FRAME_SAMPLES]
        frames.append(enc.encode_call(blk, len(blk))[0])
        pos += FRAME_SAMPLES
    while True:                                    # flush calls until the end flag
        f, end = enc.encode_call(None, 0)
        frames.append(f)
        if end:
            break
    return frames[2:]      # call 0 yields nothing (internal delay); the tool drops call 1's frame


def read_wav(path):
    if b'smpl' in open(path, 'rb').read(4096):
        print('note: input has a smpl chunk; at3tool would use its PSP loop framing, '
              'which this reference encoder does not implement (loop points ignored)', file=sys.stderr)
    with wave.open(path, 'rb') as w:
        assert w.getsampwidth() == 2 and w.getframerate() == SAMPLE_RATE, '16-bit 44.1 kHz only'
        ch = w.getnchannels()
        data = np.frombuffer(w.readframes(w.getnframes()), '<i2').reshape(-1, ch)
    return data


def encode_file(wav_in, at3_out):
    pcm = read_wav(wav_in)
    if pcm.shape[1] not in (1, 2):
        raise ValueError('mono or stereo input only')
    if len(pcm) <= 6144:
        print('note: at3tool refuses inputs of 6144 samples or fewer', file=sys.stderr)
    frames = encode_pcm(pcm)
    data = b''.join(frames)
    with open(at3_out, 'wb') as f:
        f.write(at3_header(pcm.shape[1], len(pcm), len(data)) + data)
    return len(frames)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    nf = encode_file(sys.argv[1], sys.argv[2])
    print('%d frames written' % nf)
