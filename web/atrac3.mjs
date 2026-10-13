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
import createAtrac3Module from './atrac3-core.mjs';

export async function createAtrac3() {
  const m = await createAtrac3Module();

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

  return {
    encodeWav: (wav, kbps = 0) => run(m._at3w_encode, wav, kbps),
    decodeAt3: (at3) => run(m._at3w_decode, at3),
  };
}
