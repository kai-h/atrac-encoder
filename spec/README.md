# ATRAC3 specification (MiniDisc LP2/LP4)

A specification of Sony's ATRAC3 audio codec, reconstructed from the reference implementation in Sony's `psp_at3tool.exe` for the purpose of preserving MiniDisc (MDLP) audio.

| Document | Contents |
|---|---|
| [ATRAC3.md](ATRAC3.md) | Main specification: modes, container, bitstream syntax, decoding process, encoders, tool behaviour, verification status |
| [ATRAC3-encoder-LP2.md](ATRAC3-encoder-LP2.md) | Chapter 7: Sony encoder "algorithm 1" (LP2 132 kbps, 66 kbps mono), complete and byte-exact |
| [ATRAC3-encoder-LP4.md](ATRAC3-encoder-LP4.md) | Chapter 8: Sony encoder "algorithm 0" (LP4 66 kbps joint stereo, 105 kbps, 52 kbps mono), complete and byte-exact |
| [ATRAC3-tables.md](ATRAC3-tables.md) | Appendix A: all tables (subbands, scale factors, code tables, window, QMF filter, gain and joint-stereo tables) |
| [reference/](reference/) | Executable reference models in Python: decoder (Sony-equivalent to ±1 LSB), LP2 and LP4 encoders (byte-exact) |

## Reference models

```bash
python spec/reference/atrac3_decoder.py input.at3 output.wav
python spec/reference/atrac3_encoder_lp2.py input.wav output.at3      # stereo → LP2, mono → 66 kbps
python spec/reference/atrac3_encoder_lp4.py input.wav output.at3 66   # stereo → LP4 (or 105)
```

The models need NumPy only. The decoder's output is the raw decoded stream, including the start-up delay; see ATRAC3.md §3.5 for trimming.

## How it was verified

The Sony binary runs under a CPU emulator (`tools/at3emu.py`), so its internal buffers can be compared stage by stage with the reference models. Re-running the checks needs the binary in `private/` and the Python virtual environment:

```bash
.venv/bin/python tools/gen_signals.py /tmp/sig
.venv/bin/python tools/at3emu.py -e /tmp/sig/mix.wav /tmp/sig/mix_66.at3 -br 66
.venv/bin/python tools/verify_decoder.py /tmp/sig/mix_66.at3
```

## What is not covered

- ATRAC1 (MiniDisc SP). It isn't present in this binary.
- ATRAC3plus and ATRAC Advanced Lossless. They are present in the binary but out of the current scope.
- Copy protection and the NetMD/OMA transfer containers. These wrap the same ATRAC3 frames described here.
