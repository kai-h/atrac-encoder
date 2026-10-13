// Packs the Web MiniDisc Pro plugin into a .wme file:
//   node webminidisc/package.mjs <WebAssembly build directory> <output .wme>
// A .wme is Web MiniDisc Pro's "SAR" archive: "SARFile", then the sub-magic and each file's
// name and contents, every one preceded by its big-endian length (u16 for the sub-magic and
// names, u32 for contents).  metadata.json has to be the first file.
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const [wasmDir, output] = process.argv.slice(2);
if (!output) {
  console.error('usage: node webminidisc/package.mjs <WebAssembly build directory> <output .wme>');
  process.exit(2);
}
const here = dirname(fileURLToPath(import.meta.url));
const files = [
  ['metadata.json', join(here, 'metadata.json')],
  ['index.js', join(here, 'index.js')],
  ['worker.js', join(here, 'worker.js')],
  ['atrac3.mjs', join(wasmDir, 'atrac3.mjs')],
  ['atrac3-core.mjs', join(wasmDir, 'atrac3-core.mjs')],
  ['atrac3-core.wasm', join(wasmDir, 'atrac3-core.wasm')],
  ['LICENSE', join(here, '..', 'LICENSE')],
];

function field(bytes, lengthSize) {
  const length = Buffer.alloc(lengthSize);
  length.writeUIntBE(bytes.length, 0, lengthSize);
  return [length, bytes];
}

const parts = [Buffer.from('SARFile'), ...field(Buffer.from('Web MiniDisc Dynamic Encoder'), 2)];
for (const [name, path] of files) parts.push(...field(Buffer.from(name), 2), ...field(readFileSync(path), 4));
writeFileSync(output, Buffer.concat(parts));
