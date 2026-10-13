// Module worker for index.js.  It owns the WebAssembly module and the open encoder, whose
// state carries over from one 'encode' message to the next until a message with last set ends
// the stream.  That is what makes consecutive tracks gapless.
let atrac3;
let encoder = null;

onmessage = async (ev) => {
  const { action, ...args } = ev.data;
  try {
    if (action === 'init') {
      const { createAtrac3 } = await import(args.wrapperURL);
      atrac3 = await createAtrac3({ locateFile: () => args.wasmURL });
      postMessage({});
    } else if (action === 'encode') {
      encoder ??= atrac3.createEncoder(args.kbps, 2);
      let result = encoder.write(new Int16Array(args.pcm));
      if (args.last) {
        const tail = encoder.finish();
        const joined = new Uint8Array(result.length + tail.length);
        joined.set(result);
        joined.set(tail, result.length);
        result = joined;
        encoder.free();
        encoder = null;
      }
      postMessage({ result: result.buffer }, [result.buffer]);
    }
  } catch (err) {
    postMessage({ error: String(err) });
  }
};
