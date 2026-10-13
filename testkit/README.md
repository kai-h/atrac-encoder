# Hardware encoder test kit

These are 15 test tracks (7.0 minutes) for recording on MiniDisc hardware. You then rip the raw ATRAC data back with an MZ-RH1, so we can compare each encoder's decisions and quality: Type-R, Type-S, at3tool and atracdenc.

The WAV files aren't stored in git. Generate them (deterministically, the same every time) with:

```bash
python tools/make_testkit.py          # writes testkit/*.wav and testkit/manifest.json
```

Each track begins and ends with a short chirp (a "sync" sweep). That is how the analysis finds each track in a rip and lines it up with the source to the exact sample. So **the start of a recording doesn't need to be synchronised**: silence before the music, an extra track mark or a deleted lead-in track are all fine.

| # | Track | Length | What it probes |
|---|---|---|---|
| 01 | calibration | 22 s | **Capture path:** gain, channel mapping and clock (1 kHz tone, L-only and R-only segments, noise) |
| 02 | level_steps | 34 s | 1 kHz tone from −1 to −91 dBFS in 3 dB steps (scale factors, quiet threshold) |
| 03 | noise_steps | 20 s | Pink noise from −16 to −76 dBFS |
| 04 | transients | 26 s | Bursts over a quiet tone, impulses, hard-onset tones (pre-echo, gain control, block switching) |
| 05 | attacks | 19 s | 40 silence-to-noise onsets at random positions within a frame |
| 06 | tones | 35 s | Pure tones, including either side of the QMF band edges (5.5, 11.0 and 16.5 kHz) |
| 07 | sweep | 23 s | Log sweep from 20 Hz to 21 kHz (bandwidth limit of each mode) |
| 08 | harmonic | 16 s | Harmonic complexes, steady and with vibrato |
| 09 | masking | 80 s | Loud tone plus quieter probe tones at set frequency and level offsets (the psychoacoustic model) |
| 10 | noise_spectra | 31 s | White, pink and brown noise; band-limited noise; notched noise |
| 11 | tone_in_noise | 13 s | Tone at +20 to −20 dB SNR in noise (tonal component detection) |
| 12 | stereo | 45 s | L-only, R-only, centre, anti-phase, wide, a correlation sweep and pans (joint-stereo decisions) |
| 13 | quiet | 14 s | Digital silence, ±1 LSB dither, −80 dBFS noise, −90 dBFS tone |
| 14 | pressure | 15 s | Dense full-band material (rate control under pressure) |
| 15 | music | 27 s | Synthetic song: bass, piano-like chords, hi-hats and snare |

**Optional:** add 2 or 3 of your own CD tracks after track 15, and keep lossless copies of them. Real music is the final test, and I need the exact source to compare against.

## Step 0: check the pipeline first

Before recording everything, record **track 01 only** via each route you plan to use, rip it, and send me the rips. Running `python tools/kit_align.py <rips>` confirms:
- the rip format is readable,
- the path is bit-perfect: gain +0.000 dB, clock 0 ppm,
- the channels aren't swapped.

## Burning the CD

Use either app:
- **Music app:** create a playlist with tracks 01–15 in order, then choose **File → Burn Playlist to Disc** with these settings:
  - **Audio CD**,
  - **Gap Between Songs: None**,
  - **Use Sound Check: off**.
- **Burn app:** create an **Audio** disc, add the WAVs in order and set the pause between tracks to 0.

Don't normalise, apply EQ, crossfade or convert the files. The CD player's digital output is bit-perfect.

## Recording sessions

On every deck, set the digital recording level to **0 dB** if it has that control, and switch off any level or EQ processing. In each pass, start recording, let a second or two of silence go by, then press play on the CD player and record the whole CD.

### Disc 1: MDS-E12 (Type-R), CD player over coaxial, real time

| Pass | Deck record mode | Rip folder |
|---|---|---|
| E1 | SP stereo | `rips/E12_SP` |
| E2 | Mono | `rips/E12_MONO` |
| E3 | LP2 | `rips/E12_LP2` |
| E4 | LP4 | `rips/E12_LP4` |

Disc time used is about 16 minutes (SP equivalent).

### Disc 2: MDS-S500 (Type-S)

| Pass | Route | Mode | Rip folder |
|---|---|---|---|
| S1 | NetMD from Web MiniDisc Pro (upload the 15 WAVs) | SP | `rips/S500_NetMD_SP` |
| S2 | CD player into the S500's digital input, real time | LP2 | `rips/S500_LP2` |
| S3 | CD player into the S500's digital input, real time | LP4 | `rips/S500_LP4` |
| S4 | NetMD from Web MiniDisc Pro (computer-encoded, **control**) | LP2 | `rips/WMD_LP2` |
| S5 | NetMD from Web MiniDisc Pro (computer-encoded, **control**) | LP4 | `rips/WMD_LP4` |

For S4 and S5, note which encoder Web MiniDisc Pro was set to use (atracdenc by default, or a remote encoder).

Pass S1 also tests whether the S500 encodes NetMD SP transfers itself. If the rip matches atracdenc's SP output for the same WAVs, the computer encoded it; otherwise the deck did.

Disc time used is about 18 minutes (SP equivalent).

## Ripping with the MZ-RH1

Use Web MiniDisc Pro's homebrew-mode ripping to download every track of each pass, and keep the files exactly as delivered:
- ATRAC3 arrives in WAV,
- SP arrives as `.aea`,
- `.oma` also works if it's unencrypted.

Put each pass in its own folder as listed above, for example `testkit/rips/E12_LP2/` (gitignored). For each pass, also record in a text file: the deck, mode, input, any deck settings, the date and anything unusual.

## Analysis

```bash
python tools/kit_align.py testkit/rips/E12_LP2/*            # path check + track identification
python tools/kit_align.py testkit/rips/*/* --save /tmp/aligned  # aligned per-track excerpts for analysis
```

- **Track 01** gives gain and SNR measured on the 1 kHz tone, plus the channel map. A bit-perfect digital path shows **+0.000 dB and 0 ppm**.
- **Other tracks** show a whole-track waveform fit. Its SNR describes the codec, because noise is never coded as a waveform, so it says nothing about the capture path.
