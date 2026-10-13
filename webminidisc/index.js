// Web MiniDisc Pro encoder plugin (encoder API version 1): ATRAC3 LP2, 105 kbps and LP4 with
// the WebAssembly build, including gapless uploads.  package.mjs packs it into
// atrac-encoder.wme.  Web MiniDisc Pro imports this file from a blob URL and constructs the
// default export with its plugin API; the other files come from api.getOwnFile.

function blobURL(contents, type = 'text/javascript') {
  return URL.createObjectURL(new Blob([contents], { type }));
}

// worker.js, which owns the WebAssembly module and the open encoder
class EncoderWorker {
  constructor(worker) {
    this.worker = worker;
    worker.onmessage = (ev) => this.settle(ev.data.error && new Error(ev.data.error), ev.data);
    worker.onerror = (ev) => this.settle(new Error(ev.message || 'atrac-encoder worker failed to start'));
  }

  settle(err, data) {
    const pending = this.pending;
    this.pending = undefined;
    if (err) pending?.reject(err);
    else pending?.resolve(data);
  }

  request(message, transfer = []) {
    return new Promise((resolve, reject) => {
      this.pending = { resolve, reject };
      this.worker.postMessage(message, transfer);
    });
  }

  init(wrapperURL, wasmURL) {
    return this.request({ action: 'init', wrapperURL, wasmURL });
  }

  // pcm: 16-bit little-endian stereo.  Returns the finished frames; unless last is set, the
  // encoder keeps the remaining samples and its state for the next call.
  async encode(pcm, kbps, last) {
    return (await this.request({ action: 'encode', pcm, kbps, last }, [pcm])).result;
  }

  terminate() {
    this.worker.terminate();
  }
}

export default class AtracEncoder {
  static customConstructionParameters = [];

  constructor(api, _parameters) {
    this.api = api;
    this.ffmpeg = api.getBuiltinEncoder('ffmpeg');
    this.encoderWorker = undefined;
    this.urls = [];
  }

  async init() {
    // A new batch: drop any stream an aborted batch left open.
    this.stopWorker();
    await this.ffmpeg.init();
  }

  async deinit() {
    this.stopWorker();
    await this.ffmpeg.deinit();
  }

  stopWorker() {
    this.encoderWorker?.terminate();
    this.encoderWorker = undefined;
    this.urls.forEach((url) => URL.revokeObjectURL(url));
    this.urls = [];
  }

  async startWorker() {
    // Blob URLs can't import each other relatively, so point atrac3.mjs at the core's URL.
    const coreURL = blobURL(this.api.getOwnFile('atrac3-core.mjs'));
    const wrapper = new TextDecoder().decode(this.api.getOwnFile('atrac3.mjs'));
    if (!wrapper.includes("'./atrac3-core.mjs'")) throw new Error('Unexpected atrac3.mjs');
    const wrapperURL = blobURL(wrapper.replace("'./atrac3-core.mjs'", JSON.stringify(coreURL)));
    const wasmURL = blobURL(this.api.getOwnFile('atrac3-core.wasm'), 'application/wasm');
    const workerURL = blobURL(this.api.getOwnFile('worker.js'));
    this.urls.push(coreURL, wrapperURL, wasmURL, workerURL);
    const encoderWorker = new EncoderWorker(new Worker(workerURL, { type: 'module' }));
    await encoderWorker.init(wrapperURL, wasmURL);
    this.encoderWorker = encoderWorker;
  }

  async transcode(source, sourceFileName, exportParams, transcodeCallback) {
    if (['PCM', 'MP3'].includes(exportParams.format.codec)) {
      return this.ffmpeg.transcode(source, sourceFileName, exportParams, transcodeCallback);
    }
    if (exportParams.format.codec === 'A3+') throw new Error('ATRAC3plus is not supported');

    const ffmpegCommand = this.ffmpeg.createFfmpegParams(exportParams, 's16le');
    const dotIndex = sourceFileName.lastIndexOf('.') + 1;
    const fileExtension = dotIndex > 0 ? sourceFileName.substring(dotIndex) : 'unknown';
    const inFileName = `inAudioFile.${fileExtension}`;

    transcodeCallback?.({ stage: 'ffmpeg', progress: 0, total: 2 });
    await this.ffmpeg.ffmpegProcess.write(inFileName, source);
    await this.ffmpeg.ffmpegProcess.transcode(inFileName, 'outAudioFile.raw', ffmpegCommand);
    const { data } = await this.ffmpeg.ffmpegProcess.read('outAudioFile.raw');
    await this.ffmpeg.ffmpegProcess.remove('outAudioFile.raw');
    transcodeCallback?.({ stage: 'atrac', progress: 1, total: 2 });

    if (!this.encoderWorker) await this.startWorker();
    // The buffer is transferred to the worker, so it must hold exactly the samples.
    const pcm = data.byteLength === data.buffer.byteLength ? data.buffer : data.slice().buffer;
    // writeGapless: the next track continues this one, so leave the stream open.
    const result = await this.encoderWorker.encode(pcm, exportParams.format.bitrate, !exportParams.writeGapless);
    transcodeCallback?.({ stage: 'atrac', progress: 2, total: 2 });
    return new Uint8Array(result);
  }

  getUserFriendyStageName(stage) {
    return {
      ffmpeg: 'Transcoding to PCM...',
      atrac: 'Transcoding to ATRAC...',
    }[stage] ?? this.ffmpeg.getUserFriendyStageName(stage);
  }

  getSupportFor(codec) {
    switch (codec) {
      case 'AT3': return { state: 'perfect', gapless: true, bitrates: [132, 105, 66] };
      case 'PCM':
      case 'MP3': return { state: 'perfect', gapless: true };
      default: return { state: 'unsupported', gapless: false };
    }
  }
}
