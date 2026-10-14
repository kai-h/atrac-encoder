"""LP4 signal analysis and mode configuration from Chapter 8 §§1–9."""
import numpy as np

from .lp4_analysis import QMF, MDCT
from .lp4_gain import GainControl
from .lp4_stereo import Stereo, split_budget, T


def bandwidth_subbands(hz):
    index = hz * 2048 // 44100
    return next(i + 1 for i in range(32) if index <= int(T['SB_START'][i + 1]))


class Analysis:
    def __init__(self):
        self.qmf = [QMF(), QMF()]
        self.mdct = [MDCT(), MDCT()]
        self.gain = [GainControl(), GainControl()]
        self.previous = [np.zeros((256, 4), np.float32) for _ in range(2)]
        self.stereo = Stereo()
        self.subbands = [bandwidth_subbands(hz) for hz in (12400, 2700)]

    def process(self, pcm):
        if np.shape(pcm) != (1024, 2):
            raise ValueError('LP4 requires 1024 stereo samples per analysis call')
        raw = [self.qmf[c].process(pcm[:, c]) for c in range(2)]
        converted = self.stereo.process(*raw)
        spectra, old_counts = [], []
        for c in range(2):
            modified, old_count = self.gain[c].process(converted[c], self.previous[c], c == 1,
                                                     self.stereo.sent, self.stereo.current)
            first_levels = [band.points[0][0] if band.points else 4 for band in self.gain[c].bands]
            spectra.append(self.mdct[c].process(modified, first_levels))
            old_counts.append(old_count)
        self.previous = converted
        budget = split_budget(*spectra, self.stereo.sent[0], len(self.gain[1].bands[0].points))
        return spectra, budget, old_counts, raw
