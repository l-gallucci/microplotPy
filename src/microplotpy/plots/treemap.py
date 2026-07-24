"""Taxa (or functional) treemap.

Hierarchical composition as nested rectangles -- a static, publication-ready
alternative to an interactive Krona chart: area = relative abundance (mean
across samples, so areas sum to ~100%), grouped/colored by `group_rank` with
a subgroup border, each rectangle labeled directly (no external legend
needed). Works generically on taxonomic or functional data. Mirrors R's
mp_taxa_treemap()/mp_function_treemap() in microplotr -- keep both in sync.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import squarify
from matplotlib.patches import Rectangle
from matplotlib.transforms import TransformedPatchPath

from .._long import long_abundance
from ..filter_taxa import filter_taxa
from ..io import MicrobiomeData
from ..palette import _qualitative_hues

_CANVAS = 100.0


def taxa_treemap(
    data: MicrobiomeData,
    rank: str = "Genus",
    group_rank: str | None = "Phylum",
    top_n: int | None = None,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    fix_taxonomy: bool = True,
    tax_fix_ranks: list[str] | None = None,
    other_label: str = "Other",
    figsize: tuple[float, float] = (9, 7),
):
    long = long_abundance(data, fix_taxonomy=fix_taxonomy, tax_fix_ranks=tax_fix_ranks)
    long["rel_abund"] = long.groupby("Sample_ID")["Count"].transform(lambda x: 100 * x / x.sum())

    group_cols = list(dict.fromkeys(["Sample_ID", rank] + ([group_rank] if group_rank else [])))
    agg = long.groupby(group_cols, as_index=False)["rel_abund"].sum()

    if min_rel_abund > 0 or min_prevalence > 0:
        survivors = set(filter_taxa(agg[["Sample_ID", rank, "rel_abund"]], taxon_col=rank,
                                     min_rel_abund=min_rel_abund, min_prevalence=min_prevalence,
                                     detection=detection))
    else:
        survivors = set(agg[rank].unique())

    # Mean per-sample relative abundance per taxon -- sums to ~100% across
    # taxa (linearity of the mean), so treemap area is directly comparable
    # to the barplot's per-sample relative abundance (%).
    taxon_cols = list(dict.fromkeys([rank] + ([group_rank] if group_rank else [])))
    taxon_mean = (
        agg[agg[rank].isin(survivors)]
        .groupby(taxon_cols, as_index=False)["rel_abund"].mean()
        .sort_values("rel_abund", ascending=False)
    )

    keep = set(taxon_mean[rank].iloc[:top_n])
    taxon_mean["_taxon"] = np.where(taxon_mean[rank].isin(keep), taxon_mean[rank], other_label)
    if group_rank:
        taxon_mean["_group"] = np.where(taxon_mean["_taxon"] == other_label, other_label, taxon_mean[group_rank])

    bucket_cols = ["_taxon"] + (["_group"] if group_rank else [])
    final_df = taxon_mean.groupby(bucket_cols, as_index=False)["rel_abund"].sum()

    fig, ax = plt.subplots(figsize=figsize)

    if group_rank:
        _draw_grouped(ax, final_df, other_label)
    else:
        _draw_flat(ax, final_df, other_label)

    ax.set_xlim(0, _CANVAS)
    ax.set_ylim(0, _CANVAS)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    return fig


def _label_rect(ax, patch, x, y, dx, dy, text, fontsize=8):
    # Scale down for narrow/short rectangles, and clip to the rectangle's own
    # patch so text that still doesn't fit is cropped at its own cell's
    # boundary instead of bleeding into neighboring cells or off the canvas.
    if dx <= 6 or dy <= 4:
        return
    size = max(5, min(fontsize, dx / (len(str(text)) * 0.55)))
    t = ax.text(x + dx / 2, y + dy / 2, text, ha="center", va="center", color="white", fontsize=size)
    t.set_clip_path(TransformedPatchPath(patch))


def _draw_grouped(ax, final_df, other_label):
    group_totals = final_df.groupby("_group")["rel_abund"].sum().sort_values(ascending=False)
    groups = list(group_totals.index)
    real_groups = [g for g in groups if g != other_label]
    hues = _qualitative_hues(len(real_groups))
    group_colors = dict(zip(real_groups, hues))
    if other_label in groups:
        group_colors[other_label] = "#b3b3b3"

    outer_sizes = squarify.normalize_sizes(group_totals.to_numpy().tolist(), _CANVAS, _CANVAS)
    outer_rects = squarify.squarify(outer_sizes, 0, 0, _CANVAS, _CANVAS)

    for group, orect in zip(groups, outer_rects):
        sub = final_df[final_df["_group"] == group].sort_values("rel_abund", ascending=False)
        sizes = sub["rel_abund"].tolist()
        if sum(sizes) <= 0:
            continue
        norm = squarify.normalize_sizes(sizes, orect["dx"], orect["dy"])
        inner_rects = squarify.squarify(norm, orect["x"], orect["y"], orect["dx"], orect["dy"])
        color = group_colors[group]
        for taxon, irect in zip(sub["_taxon"], inner_rects):
            rect_patch = Rectangle((irect["x"], irect["y"]), irect["dx"], irect["dy"],
                                    facecolor=color, edgecolor="white", linewidth=1)
            ax.add_patch(rect_patch)
            _label_rect(ax, rect_patch, irect["x"], irect["y"], irect["dx"], irect["dy"], taxon)

        outer_patch = Rectangle((orect["x"], orect["y"]), orect["dx"], orect["dy"],
                                 facecolor="none", edgecolor="white", linewidth=3)
        ax.add_patch(outer_patch)
        group_fontsize = max(6, min(13, orect["dx"] / (len(str(group)) * 0.6)))
        label = ax.text(orect["x"] + 1, orect["y"] + orect["dy"] - 2, str(group),
                         ha="left", va="top", color="white", fontsize=group_fontsize,
                         fontweight="bold", alpha=0.85)
        label.set_clip_path(TransformedPatchPath(outer_patch))


def _draw_flat(ax, final_df, other_label):
    final_df = final_df.sort_values("rel_abund", ascending=False)
    real_taxa = [t for t in final_df["_taxon"] if t != other_label]
    hues = iter(_qualitative_hues(len(real_taxa)))
    color_map = {t: ("#b3b3b3" if t == other_label else next(hues)) for t in final_df["_taxon"]}

    sizes = final_df["rel_abund"].tolist()
    norm = squarify.normalize_sizes(sizes, _CANVAS, _CANVAS)
    rects = squarify.squarify(norm, 0, 0, _CANVAS, _CANVAS)
    for taxon, rect in zip(final_df["_taxon"], rects):
        rect_patch = Rectangle((rect["x"], rect["y"]), rect["dx"], rect["dy"],
                                facecolor=color_map[taxon], edgecolor="white", linewidth=1)
        ax.add_patch(rect_patch)
        _label_rect(ax, rect_patch, rect["x"], rect["y"], rect["dx"], rect["dy"], taxon)


def function_treemap(
    feature_table,
    function_annotation,
    rank: str = "KEGG_ko",
    group_rank: str | None = "COG_category",
    top_n: int | None = None,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    fix_taxonomy: bool = True,
    other_label: str = "Other",
    figsize: tuple[float, float] = (9, 7),
):
    from .function_profile import _metadata_stub

    data = MicrobiomeData(
        feature_table=feature_table,
        taxonomy=function_annotation,
        metadata=_metadata_stub(feature_table, None),
    )
    tax_fix_ranks = [r for r in (group_rank, rank) if r is not None]
    return taxa_treemap(
        data, rank=rank, group_rank=group_rank, top_n=top_n,
        min_rel_abund=min_rel_abund, min_prevalence=min_prevalence, detection=detection,
        fix_taxonomy=fix_taxonomy, tax_fix_ranks=tax_fix_ranks, other_label=other_label,
        figsize=figsize,
    )
