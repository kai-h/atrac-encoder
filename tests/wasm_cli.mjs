// The command-line tool's interface on top of the WebAssembly module, for tests/run_tests.py:
//   node tests/wasm_cli.mjs <dir with atrac3.mjs> -e [-br KBPS] in.wav out.at3
//   node tests/wasm_cli.mjs <dir with atrac3.mjs> -d in.at3 out.wav
import { readFileSync, writeFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const [dir, mode, ...rest] = process.argv.slice(2);
let kbps = 0;
if (rest[0] === '-br') {
  kbps = Number(rest[1]);
  rest.splice(0, 2);
}
const [input, output] = rest;
const { createAtrac3 } = await import(pathToFileURL(join(resolve(dir), 'atrac3.mjs')).href);
const atrac3 = await createAtrac3();
const data = readFileSync(input);
writeFileSync(output, mode === '-e' ? atrac3.encodeWav(data, kbps) : atrac3.decodeAt3(data));
