"""MAG (metagenome-assembled genome) quality plots.

One dot per MAG/bin, x = Contamination (%), y = Completeness (%) -- the
standard first figure in any MAG-recovery paper. Accepts CheckM2's
`quality_report.tsv` or CheckM (v1)'s `qa`/`bin_stats` table directly
(column-name differences normalized via `normalize_mag_table()`). Draws
MIMAG quality-tier reference lines by default (Bowers et al. 2017, *Nature
Biotechnology*): high quality >=90% completeness & <5% contamination;
medium quality >=50% completeness & <10% contamination.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from ..theme import apply_pub_style
from ..validate_mag import normalize_mag_table

_MIMAG_COMPLETENESS = (50, 90)
_MIMAG_CONTAMINATION = (5, 10)


def mag_quality_plot(
    mag_table: pd.DataFrame,
    size_col: str | None = "Genome_Size",
    color_col: str | None = None,
    show_mimag_thresholds: bool = True,
    figsize: tuple[float, float] = (7, 5),
    ax: plt.Axes | None = None,
):
    mag_table = normalize_mag_table(mag_table)
    completeness = pd.to_numeric(mag_table["Completeness"])
    contamination = pd.to_numeric(mag_table["Contamination"])

    if size_col is not None and size_col not in mag_table.columns:
        size_col = None
    if color_col is not None and color_col not in mag_table.columns:
        raise ValueError(f"color_col '{color_col}' not found in mag_table.")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    sizes = None
    if size_col is not None:
        vals = pd.to_numeric(mag_table[size_col])
        sizes = (vals / vals.max()) * 300 + 20

    if color_col is not None:
        categories = mag_table[color_col].astype("category")
        cmap = plt.get_cmap("Dark2")
        colors = [cmap(i % 8) for i in categories.cat.codes]
        ax.scatter(contamination, completeness, s=sizes if sizes is not None else 60,
                   c=colors, alpha=0.85, edgecolors="none", zorder=2)
        handles = [
            plt.Line2D([], [], marker="o", linestyle="", color=cmap(i % 8), label=cat)
            for i, cat in enumerate(categories.cat.categories)
        ]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, title=None)
    else:
        ax.scatter(contamination, completeness, s=sizes if sizes is not None else 60,
                   color="steelblue", alpha=0.85, edgecolors="none", zorder=2)

    if show_mimag_thresholds:
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((xlim[0], 90), 5 - xlim[0], ylim[1] - 90,
                                facecolor="grey", alpha=0.15, zorder=0, edgecolor="none"))
        for y in _MIMAG_COMPLETENESS:
            ax.axhline(y, linestyle="--", color="grey", linewidth=1, zorder=1)
        for x in _MIMAG_CONTAMINATION:
            ax.axvline(x, linestyle="--", color="grey", linewidth=1, zorder=1)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    ax.set_xlabel("Contamination (%)")
    ax.set_ylabel("Completeness (%)")
    apply_pub_style(ax)
    fig.tight_layout()
    return fig


def mag_quality_distribution(
    mag_table: pd.DataFrame,
    bins: int = 20,
    show_mimag_thresholds: bool = True,
    figsize: tuple[float, float] = (9, 4),
):
    mag_table = normalize_mag_table(mag_table)
    completeness = pd.to_numeric(mag_table["Completeness"])
    contamination = pd.to_numeric(mag_table["Contamination"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    ax1.hist(completeness, bins=bins, color="steelblue", edgecolor="white")
    if show_mimag_thresholds:
        for x in _MIMAG_COMPLETENESS:
            ax1.axvline(x, linestyle="--", color="grey", linewidth=1)

    ax2.hist(contamination, bins=bins, color="steelblue", edgecolor="white")
    if show_mimag_thresholds:
        for x in _MIMAG_CONTAMINATION:
            ax2.axvline(x, linestyle="--", color="grey", linewidth=1)

    ax1.set_ylabel("Number of MAGs")
    apply_pub_style(ax1)
    apply_pub_style(ax2)
    # set titles *after* apply_pub_style, which blanks any existing title --
    # these act as facet-style panel labels, not a plot title, so they stay.
    ax1.set_title("Completeness (%)", fontweight="bold", fontsize=10)
    ax2.set_title("Contamination (%)", fontweight="bold", fontsize=10)
    fig.tight_layout()
    return fig
