"""Chapter 7 §10: LP2 tonal extraction, allocation, fitting, and refinement."""
import numpy as np

from .bitstream import ToneGroup, pack_unit
from .lp2_quantisation import (T, group_scales, subband_scales, scale_index,
                              extract_component, descending_order, word_lengths, quantise)


def tonal_class(groups):
    peak = int(np.argmax(groups))
    band, center = divmod(peak, 64)
    neighbor = band - 1 if center < 32 else band + 1
    other = []
    for i, value in enumerate(groups):
        b, r = divmod(i, 64)
        if b == band and abs(r - center) <= 2:
            continue
        if b == neighbor and abs(r - (63 - center)) <= 2:
            continue
        other.append(value)
    value = groups[peak]
    if max(other) < max(value - 17, 9):
        return 2 if value >= 36 else 1 if value >= 21 else 0
    return 0


def allocate(spectrum_a, spectrum_b, gains, previous):
    a = spectrum_a.copy()
    ga, gb = group_scales(a), group_scales(spectrum_b)
    tone_class = tonal_class(ga)
    peak = int(np.argmax(gb))
    nsub, nbands = 29, 3
    if tone_class:
        if peak >= 160:
            nsub, nbands = 32, 4
        elif peak >= 96:
            nsub = 30
    gain_bits = 3 * (nbands + 3 * sum(len(g.points) for g in gains[:nbands]))
    budget = 1536 - 3 * nsub - 8 - gain_bits - 13
    groups = []
    single_bits = 0
    if tone_class:
        band, center = divmod(peak, 64)
        neighbor = band - 1 if center < 32 else band + 1
        single_bits = (nbands + 6) * tone_class
        go_on = True
        for _ in range(tone_class):
            group = ToneGroup()
            groups.append(group)
            if not go_on:
                continue
            for b, r in ((band, center), (neighbor, 63 - center)):
                if not 0 <= b < nbands:
                    continue
                group.bands.add(b)
                single_bits += 12
                if single_bits >= budget:
                    group.bands.remove(b)
                    single_bits -= 12
                    continue
                c = 64 * b + r
                for g in range(max(c - 2, 0), min(c + 2, 255) + 1):
                    cell = g // 16
                    if len(group.cells.get(cell, [])) == 7:
                        continue
                    if single_bits + 36 > budget:
                        go_on = False
                        continue
                    component = extract_component(a, 4 * g)
                    group.cells.setdefault(cell, []).append(component)
                    ga[g] = scale_index(np.max(np.abs(a[4 * g:4 * g + 4])))
                    single_bits += 36
        budget -= single_bits
    slope = 60 if single_bits else 30
    group_limit = int(T['SB_START'][nsub]) // 4
    mean = float(np.float32(sum(gb[:group_limit]) / group_limit))

    def attacking(points):
        levels = [4] + [level for level, _ in points]
        return max(levels) - min(levels) > 2

    attack = (0 if any(attacking(g.points) for g in gains[:nbands]) else
              1 if any(attacking(g.points) for g in previous[:nbands]) else -1)
    tone_bits = 0
    tone_count = 0
    threshold = float(np.float32(mean + 32))
    for _ in range(0 if single_bits else 2):
        candidates = []
        cell_counts = [0] * 16
        for g in range(group_limit):
            cell = g // 16
            if gb[g] > threshold and cell_counts[cell] < 7 and len(candidates) + tone_count < 64:
                candidates.append(g)
                cell_counts[cell] += 1
        if not candidates:
            continue
        if tone_bits + nbands + 6 > budget:
            break
        tone_bits += nbands + 6
        group = ToneGroup()
        groups.append(group)
        stop = False
        for g in candidates:
            band, cell = g // 64, g // 16
            if band not in group.bands:
                if tone_bits + 12 <= budget:
                    group.bands.add(band)
                    tone_bits += 12
                else:
                    stop = True
            if band in group.bands:
                if tone_bits + 36 > budget:
                    stop = True
                else:
                    component = extract_component(a, 4 * g)
                    group.cells.setdefault(cell, []).append(component)
                    tone_count += 1
                    tone_bits += 36
                    residual = scale_index(np.max(np.abs(a[4 * g:4 * g + 4])))
                    gb[g] += residual - ga[g]
                    ga[g] = residual
        if stop:
            break
    budget -= tone_bits
    if not groups:
        budget += 2
    sa, sb = subband_scales(ga, nsub), subband_scales(gb, nsub)
    offsets = np.array([(sb[i] - mean) * float(T['OFFSET_FACTORS'][i]) for i in range(nsub)], np.float32)
    if attack != -1:
        offsets[:8] = (offsets[:8].astype(np.float64) + 0.7).astype(np.float32)
        offsets[8:18] = (offsets[8:18].astype(np.float64) + 0.5).astype(np.float32)
        offsets[0] = max(float(offsets[0]), 6)
    else:
        offsets[0] = max(float(offsets[0]), 6)
        offsets[1:4] = np.maximum(offsets[1:4], 3)
    selectors, quiet = word_lengths(offsets, sb, slope)
    fixed = [False] * nsub
    for i in range(nsub):
        if sa[i] < quiet:
            selectors[i], fixed[i] = 0, True
        if tone_bits == 0 and (selectors[i] == 0 or sa[i] < mean):
            fixed[i] = True
    original = selectors.copy()
    edges = T['SB_START']
    for i in range(nsub):
        lo, hi = int(edges[i]), int(edges[i + 1])
        inverse = float(np.float32(1 / float(T['SF'][sa[i]])))
        a[lo:hi] = (a[lo:hi].astype(np.float64) * inverse).astype(np.float32)
    density = budget / (int(edges[nsub]) + 256)
    if density > 1.15:
        limit = 0
    else:
        rounded = float(np.float32(density))
        limit = next((v for threshold, v in ((1.1, 3), (1.05, 5), (1.0, 7), (0.95, 9), (0.9, 11))
                      if rounded > threshold), 12)
    dead = [0] * nsub
    values, costs = [], []

    def trial(i, selector, dz):
        return quantise(a[int(edges[i]):int(edges[i + 1])], selector, dz)

    for i in range(nsub):
        q, cost = trial(i, selectors[i], 0)
        values.append(q)
        costs.append(cost)
    total = sum(costs)
    offset, step = 2.0, 4.0
    adjusted = offsets.copy()
    for iteration in range(15):
        weights = T['NEGATIVE_OFFSET_WEIGHTS'][:nsub]
        adjusted = (offsets.astype(np.float64) + offset * (1 if offset > 0 else weights)).astype(np.float32)
        fresh, _ = word_lengths(adjusted, sb, slope)
        fresh = [v if original[i] else 0 for i, v in enumerate(fresh)]
        for i in range(nsub):
            if fresh[i] != selectors[i]:
                q, cost = trial(i, fresh[i], dead[i])
                total += cost - costs[i]
                values[i], costs[i] = q, cost
        selectors = fresh
        if total > budget:
            offset = float(np.float32(offset - step))
        elif iteration >= 6:
            break
        else:
            offset = float(np.float32(offset + step if total < budget else offset - step))
        step = float(np.float32(step * (0.5 if iteration < 7 else 1.5)))
    for _ in range(4):
        for i in range(nsub - 1, -1, -1):
            if total <= budget:
                break
            if selectors[i] and dead[i] < limit:
                dead[i] += 1
                q, cost = trial(i, selectors[i], dead[i])
                total += cost - costs[i]
                values[i], costs[i] = q, cost
    if total > budget:
        keys = ([(sb[i] + 1) * 32 - i for i in range(nsub)] if single_bits else
                [sb[i] - i // 2 for i in range(nsub)])
        order = descending_order(keys)
        for i in reversed(order):
            while selectors[i] and total > budget:
                selectors[i] -= 1
                q, cost = trial(i, selectors[i], dead[i])
                total += cost - costs[i]
                values[i], costs[i] = q, cost
            if total <= budget:
                break
    nused = max([1] + [i + 1 for i, v in enumerate(selectors) if v])
    used = 1536 - (budget + 3 * (nsub - nused) - total)
    margins = [float(np.float32(float(adjusted[i]) - selectors[i])) if selectors[i]
               else float(np.float32(-9.99)) for i in range(nused)]
    low, high = min(margins), max(margins)
    order = []
    for k in range(1, 11):
        boundary = float(np.float32(high - k * (high - low) * 0.1))
        for i in range(nused - 1, -1, -1):
            if i not in order and margins[i] >= boundary:
                order.append(i)

    def improve(i, selector, dz):
        nonlocal used
        q, cost = trial(i, selector, dz)
        if used + cost - costs[i] > 1536:
            return False
        used += cost - costs[i]
        costs[i], values[i] = cost, q
        selectors[i], dead[i] = selector, dz
        return True

    for i in order:
        if fixed[i] or not 0 < selectors[i] < 7:
            continue
        if selectors[i] < original[i]:
            accepted = improve(i, selectors[i] + 1, dead[i])
        elif dead[i] > 0:
            accepted = improve(i, selectors[i], dead[i] - 1)
        else:
            if used > 1528:
                break
            continue
        if accepted and used > 1528:
            break
    for i in range(nused):
        if selectors[i] and dead[i]:
            improve(i, selectors[i], 0)
    for _ in range(7):
        for i in range(nused):
            if used > 1534:
                break
            if not fixed[i] and 0 < selectors[i] < original[i]:
                if not improve(i, selectors[i] + 1, dead[i]):
                    fixed[i] = True
    if nused > 1:
        order = descending_order([(sb[i] + 1) * 32 - i for i in range(nused)])
        for p in range(4):
            for i in order:
                if used > 1534:
                    break
                if 0 < selectors[i] < p + 3:
                    improve(i, selectors[i] + 1, dead[i])
    for _ in range(7):
        for i in range(nused):
            if used > 1534:
                break
            if 0 < selectors[i] < 7:
                improve(i, selectors[i] + 1, dead[i])
    for i in range(nused):
        if not selectors[i]:
            continue
        energy_original = 0.0
        energy_quantised = 0.0
        scale = float(np.float32(1 / (int(T['MAXQ'][selectors[i]]) + 0.5)))
        for x, q in zip(a[int(edges[i]):int(edges[i + 1])], values[i]):
            energy_original += float(x) * float(x)
            decoded = int(q) * scale
            energy_quantised += decoded * decoded
        if 1.25 * energy_original < energy_quantised and sa[i] > 0:
            sa[i] -= 1
    packed, actual_bits = pack_unit(nbands, gains, groups, selectors[:nused], sa[:nused], values[:nused], T)
    if actual_bits != used:
        raise AssertionError(f'Spec allocation and packing disagree: {used} versus {actual_bits} bits')
    return packed
