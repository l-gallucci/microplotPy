"""Taxa bubble (dot) plot.

One row per taxon at `rank` (top-N by total relative abundance), one column
per sample; dot size = relative abundance (%), dot color = the upper
taxonomy level `group_rank` -- the classic microbiome dot-plot convention
(psmelt + ggplot bubble / MicrobiomeAnalyst). Every sample x taxon
combination is drawn (size 0 where absent) so presence/absence patterns read
directly off the grid. Taxa (rows) are ordered grouped by `group_rank` block
(blocks and within-block taxa both ordered by descending total abundance),
matching the ordering convention used by `taxa_barplot()`.

By default (`fix_taxonomy=True`), unknown/missing values at `rank` are first
resolved via `tax_fix()`.
"""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .._long import long_abundance
from ..filter_taxa import filter_taxa
from ..io import MicrobiomeData
from ..palette import _qualitative_hues, nested_order_palette
from ..theme import apply_pub_style

_OTHER_SENTINEL = "__none__"


def taxa_bubbleplot(
    data: MicrobiomeData,
    rank: str = "Genus",
    group_rank: str | None = "Phylum",
    top_n: int | None = 25,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    facet_var: str | None = None,
    sample_order: str | list[str] | None = None,
    fix_taxonomy: bool = True,
    tax_fix_ranks: list[str] | None = None,
    max_size: float = 300.0,
    figsize: tuple[float, float] = (7, 6),
    ax: plt.Axes | None = None,
):
    """
    top_n: number of `rank` taxa to show (ranked by total relative abundance
        among taxa passing min_rel_abund/min_prevalence). None shows every
        taxon that passes the filters, uncapped.
    min_rel_abund, min_prevalence, detection: see `filter_taxa()`.
    tax_fix_ranks: passed through to tax_fix()'s `ranks` argument. None
        (default) auto-detects standard 16S rank columns; pass explicitly
        when data.taxonomy is a non-taxonomic hierarchy.
    """
    long = long_abundance(data, fix_taxonomy=fix_taxonomy, tax_fix_ranks=tax_fix_ranks)
    long["rel_abund"] = long.groupby("Sample_ID")["Count"].transform(lambda x: 100 * x / x.sum())

    agg_cols = list(dict.fromkeys(
        ["Sample_ID", rank] + ([group_rank] if group_rank else []) + ([facet_var] if facet_var else [])
    ))
    agg = long.groupby(agg_cols, as_index=False)["rel_abund"].sum()

    if min_rel_abund > 0 or min_prevalence > 0:
        survivors = set(filter_taxa(agg[["Sample_ID", rank, "rel_abund"]], taxon_col=rank,
                                     min_rel_abund=min_rel_abund, min_prevalence=min_prevalence,
                                     detection=detection))
    else:
        survivors = set(agg[rank].unique())

    taxon_total_cols = [rank] + ([group_rank] if group_rank else [])
    taxon_totals = (
        agg[agg[rank].isin(survivors)].groupby(taxon_total_cols, as_index=False)["rel_abund"].sum()
        .rename(columns={"rel_abund": "total"})
        .sort_values("total", ascending=False)
    )
    keep = set(taxon_totals[rank].iloc[:top_n])
    taxon_totals = taxon_totals[taxon_totals[rank].isin(keep)]

    if group_rank and taxon_totals[rank].duplicated().any():
        # A taxon should map to exactly one group; if the data assigns it
        # inconsistent groups across rows (e.g. messy/partial annotation),
        # collapse to the group with the highest partial abundance, but keep
        # the taxon's full cross-group total so ranking still reflects it.
        full_total = taxon_totals.groupby(rank)["total"].transform("sum")
        taxon_totals = (
            taxon_totals.assign(_full_total=full_total)
            .sort_values("total", ascending=False)
            .drop_duplicates(subset=rank, keep="first")
            .assign(total=lambda d: d["_full_total"])
            .drop(columns="_full_total")
        )

    if group_rank:
        spec = nested_order_palette(
            taxon_totals[rank], taxon_totals[group_rank], taxon_totals["total"], other_label=_OTHER_SENTINEL
        )
        taxon_order = spec["order"]
        group_of = spec["group_of"]
        groups_present = list(dict.fromkeys(group_of.values()))
        group_palette = dict(zip(groups_present, _qualitative_hues(len(groups_present))))
    else:
        taxon_order = list(taxon_totals.sort_values("total", ascending=False)[rank])
        group_of = {}
        group_palette = {}

    samples = list(dict.fromkeys(agg["Sample_ID"]))
    if sample_order is None:
        sample_levels = samples
    elif isinstance(sample_order, str):
        sample_levels = list(data.metadata.sort_values(sample_order)["Sample_ID"].astype(str))
    else:
        sample_levels = list(sample_order)

    facet_levels = [None] if not facet_var else list(dict.fromkeys(agg[facet_var]))
    n_facets = len(facet_levels)
    if ax is not None and n_facets > 1:
        raise ValueError("Cannot pass a single `ax` with facet_var set (multiple panels needed).")

    if ax is not None:
        fig = ax.figure
        axes = [ax]
    else:
        fig, axes_arr = plt.subplots(1, n_facets, figsize=(figsize[0] * n_facets, figsize[1]), sharey=True, squeeze=False)
        axes = list(axes_arr[0])

    taxon_order_top_to_bottom = list(reversed(taxon_order))
    y_pos = {t: i for i, t in enumerate(taxon_order_top_to_bottom)}
    global_max = float(agg["rel_abund"].max()) if len(agg) else 1.0

    for panel_ax, facet_level in zip(axes, facet_levels):
        sub = agg if facet_level is None else agg[agg[facet_var] == facet_level]
        panel_samples = [s for s in sample_levels if s in set(sub["Sample_ID"])]

        grid = pd.MultiIndex.from_product([panel_samples, taxon_order], names=["Sample_ID", rank]).to_frame(index=False)
        merged = grid.merge(sub[["Sample_ID", rank, "rel_abund"]], on=["Sample_ID", rank], how="left")
        merged["rel_abund"] = merged["rel_abund"].fillna(0)
        merged["x"] = merged["Sample_ID"].map({s: i for i, s in enumerate(panel_samples)})
        merged["y"] = merged[rank].map(y_pos)
        if group_rank:
            merged["_color"] = merged[rank].map(group_of).map(group_palette)
        else:
            merged["_color"] = "#1f77b4"

        sizes = merged["rel_abund"] / global_max * max_size
        panel_ax.scatter(merged["x"], merged["y"], s=sizes, c=merged["_color"])
        panel_ax.set_xticks(range(len(panel_samples)))
        panel_ax.set_xticklabels(panel_samples, rotation=45, ha="right")
        panel_ax.set_yticks(range(len(taxon_order_top_to_bottom)))
        panel_ax.set_yticklabels(taxon_order_top_to_bottom)
        apply_pub_style(panel_ax)
        if facet_level is not None:
            panel_ax.set_title(str(facet_level), fontweight="bold", fontsize=10)

    size_values = [v for v in (global_max, global_max / 2, global_max / 4) if v > 0.05]
    size_handles = [
        Line2D([], [], marker="o", linestyle="",
               markersize=math.sqrt(4 * (v / global_max * max_size) / math.pi), color="grey")
        for v in size_values
    ]
    size_labels = [f"{v:.1f}" for v in size_values]

    handles = list(size_handles)
    labels = list(size_labels)
    if group_rank:
        handles += [Patch(facecolor=group_palette[g], edgecolor="none") for g in groups_present]
        labels += list(groups_present)

    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9,
               title="Rel. abundance (%)" if not group_rank else None)

    fig.tight_layout()
    return fig
