// The command-line tool's interface on top of the WebAssembly module, for tests/run_tests.py:
//   node tests/wasm_cli.mjs <dir with atrac3.mjs> -e [-split] [-br KBPS] in.wav out.at3
//   node tests/wasm_cli.mjs <dir with atrac3.mjs> -d in.at3 out.wav
// -split encodes with the streaming encoder instead, feeding the samples in pieces of awkward
// sizes; the frames must be identical to encodeWav's, whose header is reused.
import { readFileSync, writeFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const PIECES = [1, 0, 1023, 1025, 1024, 70001, 517, 4096];   // samples per channel, cycled

const [dir, mode, ...rest] = process.argv.slice(2);
const split = rest[0] === '-split';
if (split) rest.shift();
let kbps = 0;
if (rest[0] === '-br') {
  kbps = Number(rest[1]);
  rest.splice(0, 2);
}
const [input, output] = rest;
const { createAtrac3 } = await import(pathToFileURL(join(resolve(dir), 'atrac3.mjs')).href);
const atrac3 = await createAtrac3();
const data = readFileSync(input);

// channels and samples of a 16-bit PCM WAV (test vectors only: no validation)
function readPcm(wav) {
  const dv = new DataView(wav.buffer, wav.byteOffset, wav.byteLength);
  let channels = 0;
  for (let pos = 12; pos + 8 <= wav.length; ) {
    const id = wav.toString('latin1', pos, pos + 4), size = dv.getUint32(pos + 4, true);
    if (id === 'fmt ') channels = dv.getUint16(pos + 10, true);
    if (id === 'data') return [channels, new Int16Array(new Uint8Array(wav.subarray(pos + 8, pos + 8 + size)).buffer)];
    pos += 8 + size + (size & 1);
  }
}

function encodeSplit() {
  const [channels, pcm] = readPcm(data);
  const enc = atrac3.createEncoder(kbps || (channels === 2 ? 132 : 66), channels);
  const parts = [atrac3.encodeWav(data, kbps).subarray(0, 80)];
  for (let s = 0, i = 0; s < pcm.length; i++) {
    const n = Math.min(PIECES[i % PIECES.length] * channels, pcm.length - s);
    parts.push(enc.write(pcm.subarray(s, s + n)));
    s += n;
  }
  parts.push(enc.finish());
  enc.free();
  return Buffer.concat(parts);
}

writeFileSync(output, mode === '-d' ? atrac3.decodeAt3(data) : split ? encodeSplit() : atrac3.encodeWav(data, kbps));
