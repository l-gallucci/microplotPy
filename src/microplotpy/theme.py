"""Minimal publication-figure matplotlib styling.

No title/subtitle chrome, no top/right spines -- just axis titles, tick
labels, and the legend. Plot functions never set a figure title (callers add
one externally if they want one at all).
"""
from __future__ import annotations

from matplotlib.axes import Axes


def apply_pub_style(ax: Axes) -> Axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("")
    return ax
