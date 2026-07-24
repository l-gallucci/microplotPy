"""Ordering + nested color palette shared by barplot/bubbleplot.

Mirrors R's mp_nested_order_palette() in microplotr -- keep both in sync.
Visual result matches (bold group header above the first taxon of each
group in the legend), but since a single matplotlib Text cannot mix bold
and plain font weight, the Python side returns `group_of` so the caller
builds the bold header as a separate legend entry rather than embedding
markdown in one label (as the R/ggtext side does).
"""
from __future__ import annotations

import colorsys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


def _qualitative_hues(n: int) -> list[str]:
    if n <= 0:
        return []
    base = list(plt.get_cmap("Dark2").colors)
    if n <= len(base):
        return [mcolors.to_hex(c) for c in base[:n]]
    base2 = list(plt.get_cmap("tab20").colors)
    return [mcolors.to_hex(base2[i % len(base2)]) for i in range(n)]


def _lighten(hex_color: str, amount: float) -> str:
    r, g, b = mcolors.to_rgb(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = min(l + amount * (1 - l), 1.0)
    return mcolors.to_hex(colorsys.hls_to_rgb(h, l, s))


def nested_order_palette(taxa, groups, totals, other_label="Other", other_color="#b3b3b3"):
    """Order taxa grouped by `groups` (each group contiguous, groups and
    within-group taxa both ordered by descending `totals`), assign each
    group a base hue shaded across its members, and report which taxon is
    first in its group (for the caller to render a bold header there).

    Returns dict(order, palette, group_of) -- `group_of` maps taxon -> group
    for real taxa (excludes `other_label`, which is always last / grey).
    """
    taxa = list(taxa)
    groups = list(groups)
    totals = list(totals)
    is_other = [t == other_label for t in taxa]

    main = [(t, g, tot) for t, g, tot, o in zip(taxa, groups, totals, is_other) if not o]
    other_present = any(is_other)

    group_totals: dict[str, float] = {}
    for _, g, tot in main:
        group_totals[g] = group_totals.get(g, 0.0) + tot
    group_order = sorted(group_totals, key=lambda g: -group_totals[g])

    main.sort(key=lambda x: (group_order.index(x[1]), -x[2]))

    order = [t for t, _, _ in main]
    if other_present:
        order.append(other_label)

    base_hues = dict(zip(group_order, _qualitative_hues(len(group_order))))

    palette: dict[str, str] = {}
    group_of: dict[str, str] = {}
    for g in group_order:
        members = [t for t, gg, _ in main if gg == g]
        n = len(members)
        amounts = [0.0] if n == 1 else list(np.linspace(0, 0.55, n))
        for taxon, amount in zip(members, amounts):
            palette[taxon] = _lighten(base_hues[g], amount)
            group_of[taxon] = g

    if other_present:
        palette[other_label] = other_color

    return {"order": order, "palette": palette, "group_of": group_of}
