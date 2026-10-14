"""Non-looping RIFF/WAVE and call scheduling, ATRAC3.md §§3 and 6.4."""
import struct
import wave

import numpy as np


def read_pcm(path):
    # Parse chunks, rather than searching sample bytes for a chunk name.
    raw = path.read_bytes()
    if raw[:4] != b'RIFF' or raw[8:12] != b'WAVE':
        raise ValueError('Expected a RIFF/WAVE input')
    position = 12
    while position + 8 <= len(raw):
        chunk, length = raw[position:position + 4], int.from_bytes(raw[position + 4:position + 8], 'little')
        if chunk == b'smpl':
            raise ValueError('PSP loop framing is outside the implementation scope')
        position += 8 + length + (length % 2)
    with wave.open(str(path), 'rb') as wav:
        if wav.getnchannels() not in (1, 2) or wav.getsampwidth() != 2 or wav.getframerate() != 44100:
            raise ValueError('Expected mono/stereo 16-bit PCM at 44100 Hz')
        return np.frombuffer(wav.readframes(wav.getnframes()), '<i2').reshape(-1, wav.getnchannels())


def schedule(pcm):
    if pcm.ndim != 2 or len(pcm) < 6144:
        raise ValueError('Non-looping tool framing requires at least 6144 samples per channel')
    stream = np.concatenate([np.zeros((955, pcm.shape[1]), np.int16), pcm])
    for start in range(0, len(stream), 1024):
        block = np.zeros((1024, pcm.shape[1]), np.int16)
        part = stream[start:start + 1024]
        block[:len(part)] = part
        yield block
    remainder = len(stream) % 1024
    for _ in range(2 if 0 < remainder < 886 else 3):
        yield np.zeros((1024, pcm.shape[1]), np.int16)


def container(frames, samples, channels, align, joint=False):
    data = b''.join(frames)
    if len(data) % align:
        raise ValueError('Container payload is not a whole number of frames')
    average = (44100 * align + 512) // 1024
    fmt = struct.pack('<HHIIHHHHIHHHH', 0x270, channels, 44100, average, align, 0, 14,
                      1, channels * 2048, int(joint), int(joint), 1, 0)
    fact = struct.pack('<III', samples, 1024, 1024)
    body = b'WAVE' + b'fmt ' + struct.pack('<I', len(fmt)) + fmt
    body += b'fact' + struct.pack('<I', len(fact)) + fact
    body += b'data' + struct.pack('<I', len(data)) + data
    return b'RIFF' + struct.pack('<I', len(body)) + body
