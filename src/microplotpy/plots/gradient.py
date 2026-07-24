"""ASV/taxon abundance over a continuous gradient.

Relative abundance (%) of the top-N taxa at `rank` plotted against a
continuous metadata gradient (depth, pH, time...) -- the standard way to
show which taxa respond to an environmental gradient. Same taxa filtering
and nested-legend grouping conventions as `taxa_barplot()`, so a barplot and
gradient plot of the same data select/color/order taxa consistently.
Mirrors R's mp_asv_gradient_plot() in microplotr -- keep both in sync.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from .._long import long_abundance
from ..filter_taxa import filter_taxa
from ..io import MicrobiomeData
from ..palette import _qualitative_hues, nested_order_palette
from ..theme import apply_pub_style


def _lowess_smooth(x: np.ndarray, y: np.ndarray, frac: float = 0.6, n_out: int = 100):
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        order = np.argsort(x)
        result = lowess(y[order], x[order], frac=frac, return_sorted=True)
        return result[:, 0], result[:, 1]
    except ImportError:
        # Fallback: simple polynomial fit if statsmodels isn't installed.
        order = np.argsort(x)
        xs, ys = x[order], y[order]
        deg = min(2, len(xs) - 1)
        coeffs = np.polyfit(xs, ys, deg)
        xx = np.linspace(xs.min(), xs.max(), n_out)
        return xx, np.polyval(coeffs, xx)


def taxa_gradient_plot(
    data: MicrobiomeData,
    gradient_var: str,
    rank: str = "Genus",
    group_rank: str | None = "Phylum",
    top_n: int | None = 10,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    facet: bool = False,
    smooth: bool = True,
    fix_taxonomy: bool = True,
    tax_fix_ranks: list[str] | None = None,
    other_label: str = "Other",
    figsize: tuple[float, float] = (8, 5),
):
    long = long_abundance(data, fix_taxonomy=fix_taxonomy, tax_fix_ranks=tax_fix_ranks)
    long["rel_abund"] = long.groupby("Sample_ID")["Count"].transform(lambda x: 100 * x / x.sum())

    agg_cols = list(dict.fromkeys(["Sample_ID", rank] + ([group_rank] if group_rank else [])))
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

    bucket_cols = list(dict.fromkeys(["Sample_ID", "_taxon"] + (["_group"] if group_rank else [])))
    agg = agg.groupby(bucket_cols, as_index=False)["rel_abund"].sum()

    taxon_total_cols = ["_taxon"] + (["_group"] if group_rank else [])
    taxon_totals = agg.groupby(taxon_total_cols, as_index=False)["rel_abund"].sum().rename(columns={"rel_abund": "total"})

    if group_rank and taxon_totals["_taxon"].duplicated().any():
        full_total = taxon_totals.groupby("_taxon")["total"].transform("sum")
        taxon_totals = (
            taxon_totals.assign(_full_total=full_total)
            .sort_values("total", ascending=False)
            .drop_duplicates(subset="_taxon", keep="first")
            .assign(total=lambda d: d["_full_total"])
            .drop(columns="_full_total")
        )

    if group_rank:
        spec = nested_order_palette(taxon_totals["_taxon"], taxon_totals["_group"], taxon_totals["total"],
                                     other_label=other_label)
        order, palette, group_of = spec["order"], spec["palette"], spec["group_of"]
    else:
        ordered = taxon_totals.sort_values("total", ascending=False)
        order = [t for t in ordered["_taxon"] if t != other_label]
        has_other = other_label in set(ordered["_taxon"])
        palette = dict(zip(order, _qualitative_hues(len(order))))
        if has_other:
            order = order + [other_label]
            palette[other_label] = "#999999"
        group_of = {}

    agg = agg.merge(data.metadata[["Sample_ID", gradient_var]], on="Sample_ID", how="left")
    agg[gradient_var] = pd.to_numeric(agg[gradient_var])

    if facet:
        n = len(order)
        ncols = min(3, n)
        nrows = -(-n // ncols)
        fig, axes_arr = plt.subplots(nrows, ncols, figsize=(figsize[0] * ncols / 2, figsize[1] * nrows / 2),
                                      squeeze=False, constrained_layout=True)
        axes = axes_arr.flatten()
        for ax, taxon in zip(axes, order):
            sub = agg[agg["_taxon"] == taxon]
            ax.scatter(sub[gradient_var], sub["rel_abund"], color=palette[taxon], alpha=0.8, s=30)
            if smooth and len(sub) > 3:
                xs, ys = _lowess_smooth(sub[gradient_var].to_numpy(), sub["rel_abund"].to_numpy())
                ax.plot(xs, ys, color=palette[taxon], linewidth=1.5)
            # apply_pub_style blanks any existing title, so set the panel
            # label (a facet strip, not a plot title) afterward.
            apply_pub_style(ax)
            ax.set_title(taxon, fontweight="bold", fontsize=10)
        for ax in axes[n:]:
            ax.axis("off")
        fig.supxlabel(gradient_var)
        fig.supylabel("Relative abundance (%)")
        return fig

    fig, ax = plt.subplots(figsize=figsize)
    for taxon in order:
        sub = agg[agg["_taxon"] == taxon]
        ax.scatter(sub[gradient_var], sub["rel_abund"], color=palette[taxon], alpha=0.7, s=30)
        if smooth and len(sub) > 3:
            xs, ys = _lowess_smooth(sub[gradient_var].to_numpy(), sub["rel_abund"].to_numpy())
            ax.plot(xs, ys, color=palette[taxon], linewidth=1.5)

    ax.set_xlabel(gradient_var)
    ax.set_ylabel("Relative abundance (%)")
    apply_pub_style(ax)

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
