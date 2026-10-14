# Gapless upload test set

These are two continuous signals, each cut into five tracks. They check two things about gapless uploads from Web MiniDisc Pro with the atrac-encoder plugin: that a MiniDisc deck plays the tracks back without a break, and that the disc holds exactly the frames the plugin sent. The second answers the open question for LP4, where a sound group holds two frames and a track with an odd number of frames may get padded.

The WAV files aren't stored in git. Generate them (deterministically, the same every time) with:

```bash
python tools/make_gapless_kit.py      # writes testkit/gapless/*.wav and manifest.json
```

| Set | Files | Length | Signal | Listen for |
|---|---|---|---|---|
| tone | `tone-1.wav` to `tone-5.wav` | 40 s | 440 Hz sine at −12 dBFS | any click, dropout or dip in level at a track change |
| sweep | `sweep-1.wav` to `sweep-5.wav` | 40 s | Log sweep from 20 Hz to 20 kHz at −12 dBFS | any jump or stutter in pitch at a track change |

Each signal runs continuously across its five tracks, with a fade only at the very start and end, so a click can only come from a join. The cuts fall between ATRAC3 frames, and the tracks mix odd and even frame counts: tone tracks 1, 3 and 5 and sweep tracks 2 and 3 are odd. The sweep fading out near the top of its range in LP4 is the mode's bandwidth limit, not a fault.

## Uploading

Build the plugin and install it in Web MiniDisc Pro as described under "Web MiniDisc Pro plugin" in the main [README](../../README.md), and select atrac-encoder as the encoder. For each pass below, add the five tracks of one set **in order** (check the order in the upload dialog), choose the mode, set Gapless as listed and leave ReplayGain off. ReplayGain changes the samples, and the check below would then fail.

| Pass | Set | Mode | Gapless | Rip folder |
|---|---|---|---|---|
| G1 | tone | LP2 | on | `rips/gapless/G1_tone_LP2` |
| G2 | sweep | LP2 | on | `rips/gapless/G2_sweep_LP2` |
| G3 | tone | LP4 | on | `rips/gapless/G3_tone_LP4` |
| G4 | sweep | LP4 | on | `rips/gapless/G4_sweep_LP4` |
| G5 | tone | LP2 | **off** (control) | `rips/gapless/G5_tone_LP2_off` |

Together the passes use under 2 minutes of disc time (SP equivalent).

## Listening

Play each pass from its first track straight through, on the MDS-S500 and, if you move the disc across, on the MDS-E12 too, because how the deck crosses a track mark is part of the test. In passes G1 to G4 the track changes should be inaudible. The track number changes 50 to 70 ms before the cut in the audio, which is expected. In the control pass G5, each track was encoded on its own, so there should be a short gap of 25 to 48 ms at every change. If you can't hear those gaps, the listening test can't tell either way.

Also try skipping straight to a middle track. Its first 50 to 70 ms is the end of the previous track, by design.

## Ripping and checking

Rip each pass with the MZ-RH1 (Web MiniDisc Pro's homebrew-mode download) into its folder under `testkit/rips/` (gitignored), keeping the files exactly as delivered. Then, with the native build in `build/`:

```bash
python tools/gapless_check.py testkit/rips/gapless/G3_tone_LP4/*
```

Pass the files in track order. The script works out the mode, the set and whether the upload was gapless, reproduces the frames the plugin sent and compares each track with them:

- **identical** for every track means the disc holds exactly what was sent, so the joins are continuous in the data;
- an LP4 track with an odd frame count that comes back **one frame longer** means the deck pads the last sound group of a track. The plugin would then need to hold back the last frame of odd-length tracks;
- any other difference is reported by position, for example missing frames or frames that differ from a given point.
