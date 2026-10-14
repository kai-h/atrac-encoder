"""Specification-artifact conformance tests; never import the existing codec models.

Run: python -m unittest cleanroom.test_conformance -v
"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .bitstream import (Bits, Tone, ToneGroup, pack_unit, codewords,
                        independent_frame, joint_frame)
from .encoder import encode_lp2, encode_lp4, LP4Encoder
from .framing import read_pcm, schedule, container
from .lp2_analysis import Analysis, Gain, profile, profile_points
from .lp2_quantisation import T as LP2, descending_order, scale_index
from .lp2_allocator import allocate as allocate_lp2
from .inspect_bitstream import unit as inspect_unit
from .lp4_analysis import QMF, MDCT
from .lp4_stereo import Stereo, apply_gain, split_budget
from .lp4_frontend import bandwidth_subbands
from .lp4_gain import gain_mode
from .lp4_quantisation import scale_index4, quantiser_scale, quantise as quantise_lp4
from .tables import SPEC, decode, load


class Conformance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.stages = json.loads((SPEC / 'vectors/stages.json').read_text())
        cls.cases = json.loads((SPEC / 'vectors/cases.json').read_text())['cases']
        cls.transients = read_pcm(SPEC / 'vectors/transients.wav')

    def exact_float(self, actual, expected):
        expected = decode(expected)
        self.assertEqual(actual.shape, expected.shape)
        np.testing.assert_array_equal(actual.view(np.uint32), expected.view(np.uint32))

    def test_lp2_analysis(self):
        analysis = Analysis()
        for k, expected in enumerate(self.stages['lp2']):
            with self.subTest(call=k):
                a, b, current, previous = analysis.process(self.transients[k * 1024:(k + 1) * 1024, 0])
                self.exact_float(analysis.history[:, 512:], expected['qmf_bands'])
                self.exact_float(a, expected['spectrum_A'])
                self.exact_float(b, expected['spectrum_B'])
                for actual, key in ((current, 'gain'), (previous, 'previous_gain')):
                    for g, e in zip(actual, expected[key]):
                        self.assertEqual([list(p) for p in g.points], e['points_level_location'])
                        self.exact_float(np.asarray(g.peak, np.float32), e['peak'])

    def test_lp4_analysis_kernels(self):
        qmf, transforms = [QMF(), QMF()], [MDCT(), MDCT()]
        stereo = Stereo()
        previous = [np.zeros((256, 4), np.float32) for _ in range(2)]
        for k, expected in enumerate(self.stages['lp4']):
            with self.subTest(call=k):
                raw = [qmf[c].process(self.transients[k * 1024:(k + 1) * 1024, c]) for c in range(2)]
                for c in range(2):
                    self.exact_float(raw[c], expected['qmf_before_stereo'][c])
                bands = stereo.process(*raw)
                self.assertEqual(stereo.sent, expected['sent_matrix'])
                self.assertEqual(stereo.sent_weight, expected['sent_weight'])
                # Supplied gain lists isolate the kernels; integrated detection is checked separately.
                for c in range(2):
                    gains = expected['gain'][4 * c:4 * c + 4]
                    modified = np.stack([apply_gain(previous[c][:, b], gains[b]) for b in range(4)], axis=1)
                    actual = transforms[c].process(modified, [g[0][0] if g else 4 for g in gains])
                    self.exact_float(actual, expected['spectrum'][c])
                previous = bands

    def test_lp4_integrated_stages(self):
        encoder = LP4Encoder()
        for k, expected in enumerate(self.stages['lp4']):
            with self.subTest(call=k):
                frame, budgets = encoder.process(self.transients[k * 1024:(k + 1) * 1024])
                self.assertEqual(list(budgets), expected['budgets'])
                gains = [g.points for control in encoder.analysis.gain for g in control.bands]
                self.assertEqual([[list(p) for p in g] for g in gains], expected['gain'])
                if expected['frame_hex'] is not None:
                    self.assertEqual(frame.hex(), expected['frame_hex'])

    def test_lp4_files(self):
        for case in self.cases:
            with self.subTest(signal=case['name']):
                actual = encode_lp4(read_pcm(SPEC / 'vectors' / case['pcm']))
                self.assertEqual(actual, (SPEC / 'vectors' / f'{case["name"]}.lp4.at3').read_bytes())

    def test_lp4_tone_replacement(self):
        actual = encode_lp4(read_pcm(SPEC / 'vectors/stereo.wav'))
        self.assertEqual(actual, (SPEC / 'vectors/stereo.lp4.at3').read_bytes())

    def test_lp4_configuration_and_gain_modes(self):
        self.assertEqual([bandwidth_subbands(hz) for hz in (12400, 2700, 13700)], [27, 12, 28])
        self.assertEqual(gain_mode(0, True, 3, -3), 5)
        self.assertEqual(gain_mode(0, True, 3, 3), -1)
        self.assertEqual(gain_mode(0, True, 0, 0), 0)
        self.assertEqual(gain_mode(0, True, 0, 1), 5)
        for band in range(4):
            self.assertEqual(gain_mode(band, False, 3, -3), band)

    def test_lp4_quantisation(self):
        self.assertEqual(quantiser_scale(15, 2), 2.5)
        self.assertEqual(scale_index4(np.zeros(4, np.float32)), 0)
        self.assertEqual(scale_index4(np.full(4, 65536, np.float32)), 0)
        x = np.array([1 / 3, 1, -1 / 3, -1], np.float32)
        np.testing.assert_array_equal(quantise_lp4(x, 15, 1), [0, 2, 0, -2])

    def test_lp2_files(self):
        for case in self.cases:
            with self.subTest(signal=case['name']):
                actual = encode_lp2(read_pcm(SPEC / 'vectors' / case['pcm']))
                expected = (SPEC / 'vectors' / f'{case["name"]}.lp2.at3').read_bytes()
                self.assertEqual(actual, expected)

    def test_lp2_single_tone(self):
        actual = encode_lp2(read_pcm(SPEC / 'vectors/tone.wav'))
        self.assertEqual(actual, (SPEC / 'vectors/tone.lp2.at3').read_bytes())

    def test_lp2_single_tone_band_extensions(self):
        # §10.3 decisions on crafted spectra, covering both tonal classes and band extensions.
        for group, value, bands, count in ((100, 10, 3, 1), (100, 1000, 3, 2), (200, 1000, 4, 2)):
            with self.subTest(group=group, value=value):
                spectrum = np.zeros(1024, np.float32)
                spectrum[4 * group] = value
                gain = [Gain() for _ in range(4)]
                encoded = allocate_lp2(spectrum, spectrum, gain, gain)
                parsed = inspect_unit(encoded)
                self.assertEqual(parsed['bands'], bands)
                self.assertEqual(len(parsed['groups']), count)
                self.assertLessEqual(parsed['bits'], 1536)

    def test_lp4_split_budget(self):
        expected = self.stages['lp4'][1]
        spectra = [decode(v) for v in expected['spectrum']]
        actual = split_budget(*spectra, expected['sent_matrix'][0], len(expected['gain'][4]))
        self.assertEqual(actual, expected['budgets'][0])

    def test_order_and_gain_profiles(self):
        for case in json.loads((SPEC / 'vectors/ordering.json').read_text())['cases']:
            self.assertEqual(descending_order(case['keys']), case['indices'])
        for case in json.loads((SPEC / 'vectors/gain-profiles.json').read_text())['cases']:
            values = profile(case['attacks'], case['releases'])
            self.assertEqual(values, case['profile'])
            self.assertEqual([list(p) for p in profile_points(values)], case['points_level_location'])

    def test_framing_and_container(self):
        for case in self.cases:
            pcm = read_pcm(SPEC / 'vectors' / case['pcm'])
            blocks = list(schedule(pcm))
            self.assertEqual(len(blocks) - 2, case['modes'][0]['frame_count'])
            np.testing.assert_array_equal(blocks[0][:955], 0)
            np.testing.assert_array_equal(blocks[0][955:], pcm[:69])
            np.testing.assert_array_equal(blocks[-1], 0)
            for mode in case['modes']:
                raw = (SPEC / 'vectors' / mode['raw_frames']).read_bytes()
                actual = container([raw], len(pcm), 2, mode['frame_bytes'], mode['mode'] == 'lp4')
                self.assertEqual(actual, (SPEC / 'vectors' / mode['at3']).read_bytes())
        # Original tool's accepted minimum and adjacent framing boundaries.
        for samples, count in ((6144, 8), (6213, 8), (7098, 8), (7099, 9)):
            self.assertEqual(len(list(schedule(np.zeros((samples, 2), np.int16)))) - 2, count)
        with self.assertRaises(ValueError):
            list(schedule(np.zeros((6143, 2), np.int16)))

    def test_bit_packing_and_codebooks(self):
        bits = Bits()
        bits.put(0b101, 3)
        bits.put(0b111100, 6)
        self.assertEqual(bytes(bits.bytes), bytes.fromhex('be00'))
        self.assertEqual(bits.count, 9)
        for value in (-1, 8):
            with self.assertRaises(ValueError):
                Bits().put(value, 3)
        # Every legal code from the shared decoder data must round-trip through each encoder's indexing.
        decoder = load('common')['CODES']
        for tables in (LP2, load('lp4')):
            for (mode, selector), codes in decoder.items():
                if (mode, selector) not in tables.get('ENC_CODES', tables.get('VLC', {})):
                    continue
                for text, value in codes.items():
                    values = list(value) if isinstance(value, np.ndarray) else [value]
                    words = list(codewords(values, selector, mode, tables))
                    self.assertEqual(words, [(int(text, 2), len(text))])
        self.assertEqual(joint_frame(b'ab', b'cd', 6), b'ab\0\0dc')
        with self.assertRaises(ValueError):
            independent_frame([b'abc'], 2)
        with self.assertRaises(ValueError):
            joint_frame(b'ab', b'cd', 3)

    def test_tonal_packing(self):
        # §4 fields checked by an independent literal bit sequence, including a negative CLC mantissa.
        group = ToneGroup(bands={0}, cells={0: [Tone(4, 15, [-31, 0, 1, 31])]})
        packed, used = pack_unit(1, [Gain()], [group], [0], [0], [[]], LP2)
        fields = ['101000', '00', '000', '00001', '01', '1', '011', '111', '001',
                  '001111', '000100', '100001', '000000', '000001', '011111',
                  '000', '000', '000', '00000', '0', '000']
        text = ''.join(fields)
        self.assertEqual(used, len(text))
        expected = int(text + '0' * ((-len(text)) % 8), 2).to_bytes((len(text) + 7) // 8, 'big')
        self.assertEqual(packed, expected)

    def test_scale_boundaries(self):
        for i, scale in enumerate(LP2['SF']):
            self.assertEqual(scale_index(scale), min(i + 1, 63))
            self.assertEqual(scale_index(np.nextafter(scale, np.float32(0))), i)

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            Analysis().process(np.zeros(1023))
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'invalid.wav'
            path.write_bytes(b'not a WAV')
            with self.assertRaises(ValueError):
                read_pcm(path)


if __name__ == '__main__':
    unittest.main()
