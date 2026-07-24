"""Clustered taxa heatmap.

Top-N taxa at `rank` x samples, CLR- or log10-transformed relative
abundance, with hierarchical clustering (UPGMA) on both axes -- samples
clustered by Bray-Curtis dissimilarity on relative abundance (the
community-ecology standard), taxa clustered by Euclidean distance on the
transformed values shown -- matching the `ampvis2::amp_heatmap` / `pheatmap`
/ `ComplexHeatmap` convention for marker-gene heatmaps. Dendrograms are drawn
on both margins by default (`show_dendrogram=True`).

By default (`fix_taxonomy=True`), unknown/missing values at `rank` are first
resolved via `tax_fix()`.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform

from .._long import long_abundance
from ..filter_taxa import filter_taxa
from ..io import MicrobiomeData


def taxa_heatmap(
    data: MicrobiomeData,
    rank: str = "Genus",
    top_n: int | None = 25,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    transform: str = "clr",
    pseudocount: float | None = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    hclust_method: str = "average",
    show_dendrogram: bool = True,
    fix_taxonomy: bool = True,
    tax_fix_ranks: list[str] | None = None,
    figsize: tuple[float, float] = (9, 7),
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
    if transform not in ("clr", "log10"):
        raise ValueError("transform must be 'clr' or 'log10'")

    long = long_abundance(data, fix_taxonomy=fix_taxonomy, tax_fix_ranks=tax_fix_ranks)
    long["rel_abund"] = long.groupby("Sample_ID")["Count"].transform(lambda x: 100 * x / x.sum())

    agg = long.groupby(["Sample_ID", rank], as_index=False)["rel_abund"].sum()

    if min_rel_abund > 0 or min_prevalence > 0:
        survivors = set(filter_taxa(agg, taxon_col=rank, min_rel_abund=min_rel_abund,
                                     min_prevalence=min_prevalence, detection=detection))
    else:
        survivors = set(agg[rank].unique())

    totals = agg[agg[rank].isin(survivors)].groupby(rank)["rel_abund"].sum().sort_values(ascending=False)
    keep = list(totals.index[:top_n])
    agg = agg[agg[rank].isin(keep)]

    rel_mat = agg.pivot_table(index=rank, columns="Sample_ID", values="rel_abund", fill_value=0).reindex(keep)

    pc = pseudocount if pseudocount is not None else (rel_mat.to_numpy()[rel_mat.to_numpy() > 0].min() / 2)
    if transform == "log10":
        val_mat = np.log10(rel_mat + pc)
        legend_title = "log10(rel. abundance %)"
    else:
        log_mat = np.log(rel_mat + pc)
        val_mat = log_mat - log_mat.mean(axis=0)
        legend_title = "CLR abundance"

    n_row, n_col = val_mat.shape
    row_order = list(range(n_row))
    col_order = list(range(n_col))
    row_link = col_link = None

    if cluster_rows and n_row > 2:
        row_link = linkage(val_mat.to_numpy(), method=hclust_method, metric="euclidean")
        row_order = dendrogram(row_link, no_plot=True)["leaves"]
    if cluster_cols and n_col > 2:
        bray = pdist(rel_mat.to_numpy().T, metric="braycurtis")
        col_link = linkage(bray, method=hclust_method)
        col_order = dendrogram(col_link, no_plot=True)["leaves"]

    val_mat = val_mat.iloc[row_order, col_order]
    taxa_ord = list(val_mat.index)
    sample_ord = list(val_mat.columns)
    n_row, n_col = val_mat.shape

    show_dend = show_dendrogram and (row_link is not None or col_link is not None)
    if show_dend:
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[1, 4], wspace=0.02, hspace=0.02)
        ax_top = fig.add_subplot(gs[0, 1])
        ax_left = fig.add_subplot(gs[1, 0])
        ax_heat = fig.add_subplot(gs[1, 1])

        if col_link is not None:
            dendrogram(col_link, ax=ax_top, no_labels=True, color_threshold=0, above_threshold_color="black")
        ax_top.axis("off")

        if row_link is not None:
            dendrogram(row_link, ax=ax_left, orientation="left", no_labels=True,
                       color_threshold=0, above_threshold_color="black")
            ax_left.invert_yaxis()
        ax_left.axis("off")
    else:
        fig, ax_heat = plt.subplots(figsize=figsize)

    im = ax_heat.imshow(val_mat.to_numpy(), aspect="auto", cmap="viridis", origin="upper")
    ax_heat.set_xticks(range(n_col))
    ax_heat.set_xticklabels(sample_ord, rotation=45, ha="right")
    ax_heat.set_yticks(range(n_row))
    ax_heat.set_yticklabels(taxa_ord)
    ax_heat.yaxis.tick_right()
    ax_heat.set_title("")

    # Reserve exact room for the (right-side) taxon labels, a colorbar beyond
    # them, and the (bottom, rotated) sample labels -- computed from actual
    # rendered text extents so nothing overlaps or gets clipped regardless of
    # how long the taxon/sample names are.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_w_px, fig_h_px = fig.get_size_inches() * fig.dpi

    label_w = max((t.get_window_extent(renderer=renderer).width for t in ax_heat.get_yticklabels()), default=0)
    label_h = max((t.get_window_extent(renderer=renderer).height for t in ax_heat.get_xticklabels()), default=0)
    label_w_frac = label_w / fig_w_px
    label_h_frac = label_h / fig_h_px

    cbar_width_frac = 0.03
    gap_frac = 0.02
    right_edge = max(0.3, 1 - label_w_frac - 2 * gap_frac - cbar_width_frac)
    bottom_edge = min(0.6, label_h_frac + 0.12)

    fig.subplots_adjust(right=right_edge, bottom=bottom_edge)

    heat_pos = ax_heat.get_position()
    cbar_ax = fig.add_axes([right_edge + label_w_frac + gap_frac, heat_pos.y0, cbar_width_frac, heat_pos.height])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(legend_title)

    return fig
