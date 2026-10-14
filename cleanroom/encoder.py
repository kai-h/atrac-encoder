"""LP2/LP4 non-looping file frontends, built only from the specification."""
import argparse
from pathlib import Path

from .bitstream import independent_frame, joint_frame
from .errors import SpecificationGap
from .framing import read_pcm, schedule, container
from .lp2_analysis import Analysis
from .lp2_allocator import allocate
from .lp4_frontend import Analysis as LP4Analysis
from .lp4_allocator import allocate as allocate_lp4


def encode_lp2(pcm):
    analysis = [Analysis() for _ in range(pcm.shape[1])]
    frames = []
    for call, block in enumerate(schedule(pcm)):
        # All calls advance analysis, including the two discarded call outputs.
        spectra = [analysis[c].process(block[:, c]) for c in range(pcm.shape[1])]
        units = [allocate(*stage) for stage in spectra]
        if call >= 2:
            frames.append(independent_frame(units, 192))
    return container(frames, len(pcm), pcm.shape[1], 192 * pcm.shape[1])


class LP4Encoder:
    def __init__(self):
        self.analysis = LP4Analysis()

    def process(self, block):
        analysis = self.analysis
        spectra, budget0, old_counts, _ = analysis.process(block)
        gains = [g.bands for g in analysis.gain]
        result0 = allocate_lp4(spectra[0], budget0, analysis.subbands[0], gains[0], old_counts[0])
        unit0, bits0 = result0.pack(gains[0])
        if bits0 > budget0 + 19:
            raise AssertionError(f'LP4 unit 0 exceeds allocation budget: {bits0} > {budget0 + 19}')
        budget1 = (192 - len(unit0)) * 8 - 27
        result1 = allocate_lp4(spectra[1], budget1, analysis.subbands[1], gains[1], old_counts[1],
                               joint=True, weight=analysis.stereo.sent_weight)
        header = (analysis.stereo.sent_weight, analysis.stereo.sent)
        unit1, bits1 = result1.pack(gains[1], header)
        if bits1 > budget1 + 27:
            raise AssertionError(f'LP4 unit 1 exceeds allocation budget: {bits1} > {budget1 + 27}')
        return joint_frame(unit0, unit1), (budget0, budget1)


def encode_lp4(pcm):
    if pcm.shape[1] != 2:
        raise ValueError('LP4 requires stereo input')
    encoder = LP4Encoder()
    frames = []
    for call, block in enumerate(schedule(pcm)):
        frame, _ = encoder.process(block)
        if call >= 2:
            frames.append(frame)
    return container(frames, len(pcm), 2, 192, joint=True)


def main():
    parser = argparse.ArgumentParser(description='Specification-driven LP2/LP4 encoder prototype')
    parser.add_argument('--mode', choices=('lp2', 'lp4'), default='lp2')
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    try:
        encode = encode_lp2 if args.mode == 'lp2' else encode_lp4
        output = encode(read_pcm(args.input))
    except (SpecificationGap, ValueError) as error:
        parser.exit(2, f'{error}\n')
    args.output.write_bytes(output)


if __name__ == '__main__':
    main()
