"""Stacked relative-abundance taxa barplot, optionally with a nested legend.

Standard 16S/marker-gene stacked barplot: one bar per sample, segments are
relative abundance (%) of the top-N taxa at `rank` (everything else pooled
into `other_label`). With `nested_legend=True` (default), taxa are ordered
and colored grouped by `group_rank` -- each group gets a base hue, its member
taxa get shades of that hue, and the legend shows a bold group header above
the first taxon of each group -- matching the grouped-legend convention used
e.g. in the R `ggnested` package for microbiome barplots. Bars/legend are
ordered by abundance at `rank` (default Genus), grouped by `group_rank`
(default Phylum) abundance.

By default (`fix_taxonomy=True`), unknown/missing values at `rank` are first
resolved via `tax_fix()` so no bar segment is left unlabeled.

The plot carries no title -- only axis titles and the legend.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from .._long import long_abundance
from ..filter_taxa import filter_taxa
from ..io import MicrobiomeData
from ..palette import _qualitative_hues, nested_order_palette
from ..theme import apply_pub_style


def taxa_barplot(
    data: MicrobiomeData,
    rank: str = "Genus",
    group_rank: str | None = "Phylum",
    top_n: int | None = 10,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    facet_var: str | None = None,
    sample_order: str | list[str] | None = None,
    fix_taxonomy: bool = True,
    tax_fix_ranks: list[str] | None = None,
    nested_legend: bool = True,
    other_label: str = "Other",
    figsize: tuple[float, float] = (7, 5),
    ax: plt.Axes | None = None,
):
    """
    top_n: number of `rank` taxa to show (ranked by total relative abundance
        among taxa passing min_rel_abund/min_prevalence). None shows every
        taxon that passes the filters, uncapped.
    min_rel_abund: minimum mean relative abundance (%) across all samples for
        a taxon to be eligible; taxa below this are pooled into other_label
        regardless of top_n. 0 (default) disables this filter.
    min_prevalence: minimum fraction of samples (0-1) -- or, if > 1, an
        absolute sample count -- in which a taxon must be detected (relative
        abundance > detection) to be eligible. 0 (default) disables this filter.
    detection: relative-abundance (%) threshold above which a taxon counts as
        "detected" in a sample, for min_prevalence. Default 0.
    tax_fix_ranks: passed through to tax_fix()'s `ranks` argument. None
        (default) auto-detects standard 16S rank columns; pass explicitly
        (e.g. [group_rank, rank]) when data.taxonomy is a non-taxonomic
        hierarchy (e.g. a functional annotation table).
    """
    if nested_legend and group_rank is None:
        raise ValueError(
            "nested_legend=True requires group_rank (e.g. 'Phylum'); pass group_rank=None "
            "together with nested_legend=False for a flat legend."
        )

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

    totals = agg[agg[rank].isin(survivors)].groupby(rank)["rel_abund"].sum().sort_values(ascending=False)
    keep = set(totals.index[:top_n])

    agg["_taxon"] = np.where(agg[rank].isin(keep), agg[rank], other_label)
    if group_rank:
        agg["_group"] = np.where(agg["_taxon"] == other_label, other_label, agg[group_rank])

    bucket_cols = list(dict.fromkeys(
        ["Sample_ID", "_taxon"] + (["_group"] if group_rank else []) + ([facet_var] if facet_var else [])
    ))
    agg = agg.groupby(bucket_cols, as_index=False)["rel_abund"].sum()

    taxon_total_cols = ["_taxon"] + (["_group"] if group_rank else [])
    taxon_totals = (
        agg.groupby(taxon_total_cols, as_index=False)["rel_abund"].sum().rename(columns={"rel_abund": "total"})
    )

    if group_rank and taxon_totals["_taxon"].duplicated().any():
        # A taxon should map to exactly one group; if the data assigns it
        # inconsistent groups across rows (e.g. messy/partial annotation),
        # collapse to the group with the highest partial abundance, but keep
        # the taxon's full cross-group total so ranking still reflects it.
        full_total = taxon_totals.groupby("_taxon")["total"].transform("sum")
        taxon_totals = (
            taxon_totals.assign(_full_total=full_total)
            .sort_values("total", ascending=False)
            .drop_duplicates(subset="_taxon", keep="first")
            .assign(total=lambda d: d["_full_total"])
            .drop(columns="_full_total")
        )

    if nested_legend:
        spec = nested_order_palette(
            taxon_totals["_taxon"], taxon_totals["_group"], taxon_totals["total"], other_label=other_label
        )
        order, palette, group_of = spec["order"], spec["palette"], spec["group_of"]
    else:
        ordered = taxon_totals.sort_values("total", ascending=False)
        order = [t for t in ordered["_taxon"] if t != other_label]
        has_other = other_label in set(ordered["_taxon"])
        hues = _qualitative_hues(len(order))
        palette = dict(zip(order, hues))
        if has_other:
            order = order + [other_label]
            palette[other_label] = "#b3b3b3"
        group_of = {}

    if sample_order is None:
        sample_levels = list(dict.fromkeys(agg["Sample_ID"]))
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

    for panel_ax, facet_level in zip(axes, facet_levels):
        sub = agg if facet_level is None else agg[agg[facet_var] == facet_level]
        wide = sub.pivot_table(index="Sample_ID", columns="_taxon", values="rel_abund", fill_value=0)
        present_samples = [s for s in sample_levels if s in wide.index]
        wide = wide.reindex(index=present_samples, columns=order, fill_value=0)

        x = np.arange(len(wide.index))
        bottom = np.zeros(len(wide.index))
        for taxon in order:
            heights = wide[taxon].to_numpy()
            panel_ax.bar(x, heights, bottom=bottom, width=0.9, color=palette[taxon])
            bottom = bottom + heights

        panel_ax.set_xticks(x)
        panel_ax.set_xticklabels(wide.index, rotation=45, ha="right")
        apply_pub_style(panel_ax)
        if facet_level is not None:
            panel_ax.set_title(str(facet_level), fontweight="bold", fontsize=10)

    axes[0].set_ylabel("Relative abundance (%)")

    handles, labels = [], []
    seen_groups = set()
    for taxon in order:
        g = group_of.get(taxon)
        if g is not None and g not in seen_groups:
            handles.append(Patch(facecolor="none", edgecolor="none"))
            labels.append(g)
            seen_groups.add(g)
        handles.append(Patch(facecolor=palette[taxon], edgecolor="none"))
        labels.append(("    " if g is not None else "") + taxon)

    legend = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
    group_names = set(group_of.values())
    for text in legend.get_texts():
        if text.get_text() in group_names:
            text.set_fontweight("bold")

    fig.tight_layout()
    return fig
