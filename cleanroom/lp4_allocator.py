"""Chapter 8 §10: LP4 spectrum and tonal-component allocation."""
import math
from dataclasses import dataclass

import numpy as np

from .bitstream import Tone, ToneGroup, pack_unit
from .lp4_quantisation import T, scale_index4, quantise, bit_cost, extract_tone, restore_tone
from .errors import SpecificationGap


@dataclass
class Allocation:
    bands: int
    selectors: list
    scales: list
    groups: list
    residual: np.ndarray

    def pack(self, gains, header=None):
        values = []
        edges = T['SB_START']
        for i, selector in enumerate(self.selectors):
            segment = self.residual[int(edges[i]):int(edges[i + 1])]
            values.append(quantise(segment, self.scales[i], selector) if selector else [])
        return pack_unit(self.bands, gains, self.groups, self.selectors, self.scales, values, T, header)


def clamp_selector(value):
    return min(7, max(1, value))


def all_tonal(x, nsub, nbands, left, gsf, sf, selectors, priorities):
    groups = [ToneGroup(selector=5), ToneGroup(selector=7)]
    cost, count = 22, 0
    edges = T['SB_START']
    done = False
    while not done:
        initial_count = count
        k = 0
        while k < int(edges[nsub]):
            if abs(float(x[k])) < 2:
                k += 1
                continue
            position = min(k, 1020)
            scale = scale_index4(x[position:position + 4])
            gi = int(scale > 30)
            cell = k >> 6
            reuse = None
            if len(groups[gi].cells.get(cell, [])) == 7:
                if len(groups[1 - gi].cells.get(cell, [])) == 7:
                    records = groups[gi].cells[cell]
                    eligible = [j for j, r in enumerate(records) if r.scale < scale]
                    if not eligible:
                        k += 1
                        continue
                    weakest = min(eligible, key=lambda j: records[j].scale)
                    reuse = records.pop(weakest)
                    restore_tone(x, reuse, groups[gi].selector)
                    position = reuse.position
                else:
                    gi = 1 - gi
            group = groups[gi]
            band = cell >> 2
            if band not in group.bands:
                cost += 12
                group.bands.add(band)
            record, _ = extract_tone(x, position, group.selector)
            group.cells.setdefault(cell, []).append(record)
            if reuse is None:
                count += 1
            cost += 28 + 8 * gi
            # Replacements reserve another component's bits but add no new record.
            if cost > left - 200 or count >= 64:
                done = True
                break
            k += 4
        if count == initial_count:
            break
    for i in range(nsub):
        start, end = int(edges[i]) // 4, int(edges[i + 1]) // 4
        for g in range(start, end):
            gsf[g] = scale_index4(x[4 * g:4 * g + 4])
        sf[i] = max(gsf[start:end])
        selectors[i] = int(sf[i] != 0)
        priorities[i] = 341 * sf[i]
    return groups, cost


def extract_normal_tones(x, nsub, nbands, left, sf, selectors, priorities):
    edges, widths = T['SB_START'], T['SB_WIDTH']
    sfvalues = T['SF_BITS'].astype(np.uint32).view(np.float32)
    group = ToneGroup(selector=3, mode=0)
    records = []
    header, bits, zeros = 8, 0, 0
    for i in range(8 if left < 1100 else 0, nsub):
        if len(records) > 64 - int(widths[i]) // 16 or 24 * len(records) + header > left - 600:
            break
        if not selectors[i]:
            continue
        cls = 0 if priorities[i] < 0x500 else selectors[i]
        band = 1 if int(widths[i]) == 32 else int(edges[i]) >> 8
        threshold = max(0, sf[i] - int(T['TONE_SF_OFFSET'][cls + 8 * band]))
        start, end = int(edges[i]), int(edges[i + 1])
        if records and records[-1].position + 4 > start:
            sf[i] = -sf[i]
            overlap = records[-1].position - start
            selectors[i] = overlap + 4
            selectors[i - 1] += overlap
            start = records[-1].position + 4
        peaks = []
        k = start
        while True:
            if abs(float(x[k])) >= float(sfvalues[threshold]):
                peaks.append(k)
                k += 3
            k += 1
            if k >= end:
                break
        if len(peaks) > (int(widths[i]) + 8) // 16:
            peaks = []
        if peaks and sf[i] >= 0:
            sf[i], selectors[i] = -sf[i], 0
        for p in peaks:
            position, covered = min(p, 1020), 4
            if records and records[-1].position + 4 == position:
                while covered > 1 and abs(float(x[position + 3])) <= abs(float(x[position - 1])):
                    position -= 1
                    covered -= 1
            cell = position >> 6
            if len(group.cells.get(cell, [])) == 7:
                continue
            record, cost = extract_tone(x, position, 3)
            if record.values[0] == 0:
                skip = next((j for j in range(1, 4) if record.values[j]), 4)
                if position + skip >= end:
                    restore_tone(x, record, 3)
                    break
                if position + skip <= 1020:
                    record.position += skip
                    record.values = record.values[skip:] + [0] * skip
            bits += cost
            selectors[i] += covered
            cell = record.position >> 6
            if cell >> 2 not in group.bands:
                group.bands.add(cell >> 2)
                header += 12
            group.cells.setdefault(cell, []).append(record)
            records.append(record)
            zeros += record.values[3] == 0
    count = len(records)
    if count == 0:
        return [], 0
    if count == 1 and (records[0].position < 128 or left < 1100):
        restore_tone(x, records[0], 3)
        return [], 0
    header += nbands
    vlc, clc = bits, 24 * count
    shrink = zeros == count
    if not shrink:
        nz, fields = 2 * count - 3 * zeros, 12 * (count - zeros)
        ve, ce = nz + fields + bits, 3 * nz + fields + 24 * count
        if min(bits, 24 * count) > min(ve, ce):
            additions = []
            cell_counts = {cell: len(values) for cell, values in group.cells.items()}
            for record in records:
                if record.values[3]:
                    position = record.position + 2
                    cell = position >> 6
                    if (position > 1021 or cell >> 2 not in group.bands
                            or cell_counts.get(cell, 0) >= 7 or count + len(additions) == 64):
                        additions = None
                        break
                    cell_counts[cell] = cell_counts.get(cell, 0) + 1
                    additions.append((record, Tone(position, record.scale, record.values[2:] + [0, 0])))
            if additions is not None:
                for old, new in additions:
                    old.values[2:] = [0, 0]
                    records.append(new)
                    group.cells.setdefault(new.position >> 6, []).append(new)
                count = len(records)
                vlc, clc, shrink = count + ve, 3 * count + ce, True
    if shrink:
        v = 3
        while True:
            stop = False
            for record in reversed(records):
                if record.values[v]:
                    record.position += 1
                    record.values[:v] = record.values[1:v + 1]
                if record.values[v - 1] and (record.values[0] or record.position % 64 == 63):
                    stop = True
            v -= 1
            if stop:
                break
            if v < 1:
                raise SpecificationGap('LP4 §10.4: shrink termination for an all-zero component group is undefined')
        removed = (3 - v) * count
        vlc -= removed
        clc -= 3 * removed
        group.value_count = v + 1
        for record in records:
            record.values = record.values[:v + 1]
    group.mode = int(vlc > clc)
    return [group], header + min(vlc, clc)


def allocate(spectrum, budget, nsub_max, gains, previous_count0, joint=False, weight=7):
    x = spectrum.reshape(1024).copy()
    if budget < 40:
        gains[0].points = []
        return Allocation(1, [0], [0], [], x)
    edges, widths = T['SB_START'], T['SB_WIDTH']
    gsf = [scale_index4(x[4 * g:4 * g + 4]) for g in range(256)]
    sf, sums, cumulative = [0] * 32, [0] * 32, [0] * 32
    selectors, priorities = [0] * 32, [0] * 32
    nsub, count, total, top = nsub_max, 0, 0, 0
    for i in range(32):
        part = gsf[int(edges[i]) // 4:int(edges[i + 1]) // 4]
        level_sum, maximum = sum(part), max(part)
        count += sum(v > 7 for v in part)
        if top < maximum:
            top = maximum
            if not joint:
                multiple = 16 if int(edges[i]) & 0x80 and i < 22 else 8
                if count * multiple < budget:
                    nsub = max(nsub, 28 if i <= 26 else min(i + 2, 32))
        elif level_sum < int(T['MIN_SFSUM'][i]) or maximum < 3:
            level_sum = maximum = 0
        sf[i], sums[i], cumulative[i] = maximum, level_sum, count
        total += level_sum
    shift = 11 if 2 * budget < 12 * cumulative[nsub - 1] else 10
    nbands = (int(edges[nsub]) + 255) // 256
    left = budget - 3 * nsub - 3 * (3 * sum(len(g.points) for g in gains[:nbands]) + nbands)
    if left > 32 * cumulative[nsub - 1]:
        groups, cost = all_tonal(x, nsub, nbands, left, gsf, sf, selectors, priorities)
        left -= cost
    else:
        for i in range(nsub):
            m = sf[i]
            if m >= 3:
                m = 256 * m - total
                selectors[i] = clamp_selector(m >> shift)
            priorities[i] = m
        groups, cost = ([], 0) if joint else extract_normal_tones(x, nsub, nbands, left, sf, selectors, priorities)
        left -= cost
        for i in range(nsub):
            m = sf[i]
            if m < 0:
                for g in range(int(edges[i]) // 4, int(edges[i + 1]) // 4):
                    gsf[g] = scale_index4(x[4 * g:4 * g + 4])
                part = gsf[int(edges[i]) // 4:int(edges[i + 1]) // 4]
                maximum, level_sum = max(part), sum(part)
                sf[i] = maximum
                if selectors[i] == 0:
                    priorities[i] += 256
                else:
                    k = min((selectors[i] - 1) >> 1, 4)
                    penalty = int(T['TONE_PENALTY'][k]) * (maximum + m)
                    if budget < 1100:
                        penalty *= 2
                    elif int(widths[i]) <= 16 and priorities[i] > 0x1800:
                        penalty = min(2 * penalty, priorities[i] - 0x400)
                    if i > 1:
                        priorities[i - 1] -= penalty >> 3
                    if i < 31:
                        priorities[i + 1] -= penalty >> 2
                    priorities[i] -= penalty
                selectors[i] = (clamp_selector(priorities[i] >> shift)
                                if maximum > 7 and level_sum > int(T['MIN_SFSUM'][i])
                                and (priorities[i] >= 0x400 or budget > 0x708) else 0)
            elif shift == 11:
                if int(widths[i]) * m < 5 * sums[i]:
                    if m >= 62:
                        continue
                    sf[i] = m + 1
                if sf[i] < 62:
                    sf[i] += 1
            elif int(widths[i]) * m > 8 * sums[i]:
                priorities[i] += 256

    def estimate(i):
        selector = selectors[i]
        quiet = sum(gsf[g] < sf[i] - int(T['EST_THR'][selector])
                    for g in range(int(edges[i]) // 4, int(edges[i + 1]) // 4))
        return int(T['EST_BITS'][selector]) * int(widths[i]) + 60 - int(T['EST_SAVE'][selector]) * quiet

    width, full_width, estimate_total = 0, 0, 0
    for i in range(nsub):
        if not selectors[i]:
            continue
        if ((i > 1 and priorities[i] < priorities[i - 1] - 0xf00)
                or (i < nsub - 1 and priorities[i] < priorities[i + 1] - 0x1400)):
            selectors[i] = 0
        else:
            width += int(widths[i])
            if selectors[i] == 7:
                full_width += int(widths[i])
            estimate_total += estimate(i)
    while nsub > 1 and selectors[nsub - 1] == 0:
        nsub -= 1
        left += 3
    if not joint or weight & 7 == 7:
        old_bands = nbands
        while nbands > 1 and not gains[nbands - 1].points and not any(nbands - 1 in g.bands for g in groups):
            nbands -= 1
        left += (len(groups) + 3) * (old_bands - nbands)
    l10 = 10 * left
    denominator = width if l10 <= estimate_total else width - full_width
    offset = 0
    if denominator:
        diff = (l10 if left < 600 else 9 * left) - estimate_total
        offset = math.trunc(1024 * diff / (10 * denominator))
        if diff < 0 and shift == 10:
            offset += offset >> 3
            priorities[0] -= offset
            priorities[1] -= (188 * offset) >> 8
            for k in range(2, 8):
                priorities[k] -= (97 * offset) >> 8
            for k in range(8, 18):
                priorities[k] -= (74 * offset) >> 8
            for k in range(18):
                if priorities[k] + offset < 0 and sf[k] < 63:
                    sf[k] += 1
            for k in range(18, nsub):
                if priorities[k] + offset < 1024:
                    sf[k] = min(63, sf[k] - ((priorities[k] + offset - 1024) >> 10))
    for i in range(nsub - 1, -1, -1):
        if selectors[i]:
            selectors[i] = clamp_selector((priorities[i] + offset) >> shift)
            l10 -= estimate(i)
            priorities[i] = l10
    used, order = 0, [-1] * 32

    def bits(i, selector=None):
        return bit_cost(x[int(edges[i]):int(edges[i + 1])], sf[i], selectors[i] if selector is None else selector)

    for i in range(nsub):
        if not selectors[i]:
            continue
        if int(widths[i]) > 2 * (left - used) - 16:
            selectors[i] = 0
            continue
        diff = priorities[i] - 10 * used
        if selectors[i] < 7 and diff > 4000:
            selectors[i] += 1
        elif i != 0 and diff < 300:
            if selectors[i] > 1:
                selectors[i] -= 1
            elif sf[i] != 63:
                sf[i] += 1
        b = bits(i)
        priorities[i] = b
        zero_cost = int(T['ZERO_LEN2'][selectors[i]]) * (int(widths[i]) >> 1) + 6
        if b == zero_cost:
            selectors[i] = 0
        else:
            used += b
            order[i] = sf[i]
    while nsub > 1 and selectors[nsub - 1] == 0:
        nsub -= 1
        left += 3
    slack = left - used
    if not joint and previous_count0 + len(gains[0].points):
        for k in range(5):
            if 0 <= order[k] <= 64:
                order[k] = 64
    spent = 0
    while max(order) >= 0:
        i = order.index(max(order))
        order[i] = -1
        b = priorities[i]
        if spent + b > left:
            # The clearing loop starts at the last subband and stops at a processed one;
            # an unprocessed band somewhere below that is not sufficient to enter it.
            if nsub - 1 > i and order[nsub - 1] >= 0:
                raise SpecificationGap('LP4 §10.9: unreachable clearing branch encountered')
            zero_cost = int(T['ZERO_LEN2'][selectors[i]]) * (int(widths[i]) >> 1) + 6
            while spent + b > left:
                sf[i] += 1
                if sf[i] >= 64 or b <= zero_cost:
                    break
                b = bits(i)
        elif selectors[i] < 7 and int(T['STEP_COST'][selectors[i]]) * (int(widths[i]) >> 1) < left - spent - 7:
            new_selector = selectors[i] + 1
            if slack > 4 * int(widths[i]) and new_selector < 7:
                new_selector += 1
            new_cost = bits(i, new_selector)
            if new_cost - b <= max(slack, 0):
                selectors[i] = new_selector
                slack += b - new_cost
                b = new_cost
        spent += b
        if spent > left:
            spent -= b
            slack += b
            selectors[i] = 0
    return Allocation(nbands, selectors[:nsub], sf[:nsub], groups, x)
