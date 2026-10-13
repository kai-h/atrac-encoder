// ATRAC3 (MiniDisc LP2/LP4) encoder and decoder for the browser and Node, compiled from the
// C sources to WebAssembly.  Output is byte-identical to Sony's at3tool.
//
//   import { createAtrac3 } from './atrac3.mjs';
//   const atrac3 = await createAtrac3();
//   const at3 = atrac3.encodeWav(wavBytes, 132);   // Uint8Array .at3 file
//   const wav = atrac3.decodeAt3(at3);             // Uint8Array 16-bit WAV file
//
// encodeWav takes a 16-bit 44.1 kHz mono or stereo WAV file.  kbps: 132 (LP2), 105 or 66
// (LP4) for stereo, 66 or 52 for mono; 0 selects 132 for stereo and 66 for mono.  The
// ATRAC3 frames start at byte 80 of the .at3 file (a RIFF/WAVE header precedes them).
//
// The streaming encoder takes raw PCM and returns raw ATRAC3 frames, no header.  Its state
// carries over between write calls, so audio fed in pieces (such as gapless album tracks)
// encodes exactly as one continuous stream would; only finish() ends the stream:
//
//   const enc = atrac3.createEncoder(132, 2);   // kbps and channels as for encodeWav, no 0
//   const a = enc.write(pcm1);   // Int16Array of interleaved samples -> Uint8Array of frames
//   const b = enc.write(pcm2);   // continues pcm1 seamlessly
//   const c = enc.finish();      // the remaining frames
//   enc.free();
//
// Frames come out as soon as they are complete, with at3tool's framing: decoded sample i of
// the joined frames is sample i - 1093 of the joined input.  createAtrac3(options) passes
// options to the Emscripten module factory, e.g. { locateFile: () => url } when
// atrac3-core.wasm isn't beside atrac3-core.mjs.
import createAtrac3Module from './atrac3-core.mjs';

const CHUNK = 65536;            // samples per channel copied into the module at a time

export async function createAtrac3(options = {}) {
  const m = await createAtrac3Module(options);

  function run(fn, input, ...args) {
    const ptr = m._malloc(input.length);
    m.HEAPU8.set(input, ptr);
    try {
      const err = fn(ptr, input.length, ...args);
      if (err) throw new Error(m.UTF8ToString(err));
      const out = m._at3w_result();
      return m.HEAPU8.slice(out, out + m._at3w_result_len());
    } finally {
      m._free(ptr);
    }
  }

  function createEncoder(kbps, channels) {
    let enc = m._at3_encoder_create(kbps, channels);
    if (!enc) throw new Error(`unsupported mode: ${kbps} kbps with ${channels} channel(s)`);
    const frameBytes = m._at3_encoder_frame_bytes(enc);
    const pcmPtr = m._malloc(CHUNK * channels * 2);
    const outPtr = m._malloc((CHUNK / 1024 + 4) * frameBytes);  // a chunk's frames, or finish's 4

    function check() {
      if (!enc) throw new Error('encoder has been freed');
    }

    return {
      frameBytes,
      write(pcm) {
        check();
        if (!(pcm instanceof Int16Array) || pcm.length % channels)
          throw new Error('write takes an Int16Array of interleaved samples');
        const n = pcm.length / channels;
        const out = new Uint8Array(m._at3_encoder_max_output(enc, n));
        let pos = 0;
        for (let s = 0; s < n; s += CHUNK) {
          const piece = pcm.subarray(s * channels, Math.min(n, s + CHUNK) * channels);
          m.HEAPU8.set(new Uint8Array(piece.buffer, piece.byteOffset, piece.byteLength), pcmPtr);
          const w = m._at3_encoder_write(enc, pcmPtr, piece.length / channels, outPtr);
          out.set(m.HEAPU8.subarray(outPtr, outPtr + w), pos);
          pos += w;
        }
        return pos === out.length ? out : out.slice(0, pos);  // own buffer, exactly the frames
      },
      finish() {
        check();
        return m.HEAPU8.slice(outPtr, outPtr + m._at3_encoder_finish(enc, outPtr));
      },
      free() {
        if (!enc) return;
        m._at3_encoder_destroy(enc);
        m._free(pcmPtr);
        m._free(outPtr);
        enc = 0;
      },
    };
  }

  return {
    encodeWav: (wav, kbps = 0) => run(m._at3w_encode, wav, kbps),
    decodeAt3: (at3) => run(m._at3w_decode, at3),
    createEncoder,
  };
}
