"""Assembly QC plots: Nx curve and summary barplot.

Mirrors R's mp_assembly_nx_plot()/mp_assembly_summary_barplot() in
microplotr -- keep both in sync.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from ..theme import apply_pub_style
from ..validate_assembly import normalize_assembly_summary


def assembly_nx_plot(contig_lengths: pd.DataFrame, figsize: tuple[float, float] = (7, 5), ax: plt.Axes | None = None):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    for assembly_id, group in contig_lengths.groupby("Assembly_ID"):
        lengths = group["Length"].astype(float).sort_values(ascending=False).to_numpy()
        total = lengths.sum()
        cum_pct = lengths.cumsum() / total * 100
        ax.step(cum_pct, lengths, where="post", label=str(assembly_id))

    ax.set_yscale("log")
    ax.set_xlabel("Nx (% of total assembly length)")
    ax.set_ylabel("Contig length (bp)")
    apply_pub_style(ax)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    return fig


def assembly_summary_barplot(
    assembly_summary: pd.DataFrame,
    stat_col: str = "N50",
    figsize: tuple[float, float] = (6, 4),
    ax: plt.Axes | None = None,
):
    assembly_summary = normalize_assembly_summary(assembly_summary)
    if stat_col not in assembly_summary.columns:
        raise ValueError(f"stat_col '{stat_col}' not found in assembly_summary.")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    values = pd.to_numeric(assembly_summary[stat_col])
    ax.bar(assembly_summary["Assembly_ID"].astype(str), values, color="steelblue", width=0.6)
    ax.set_ylabel(stat_col.replace("_", " "))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    apply_pub_style(ax)
    fig.tight_layout()
    return fig
