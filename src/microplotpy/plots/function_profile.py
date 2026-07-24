"""Functional profile barplot/heatmap.

Thin wrappers around taxa_barplot()/taxa_heatmap() -- identical engine,
different input names and defaults suited to functional-annotation tools.

Realistic input shape: a `feature_table` of per-sample gene/KO counts
(from any read quantifier, e.g. featureCounts/Salmon on predicted ORFs)
joined by Feature_ID to a `function_annotation` table from a gene
functional-annotation tool:

- **eggNOG-mapper** (`*.emapper.annotations`): use rank="KEGG_ko",
  group_rank="COG_category" (defaults) -- eggNOG-mapper's own column names,
  unchanged.
- **KofamScan** (KO-only, no broader category): pass a function_annotation
  with just Feature_ID/KO/KO_definition, set rank="KO", group_rank=None,
  nested_legend=False (flat legend, no grouping tier to nest under).

Missing/blank category values (e.g. a gene with no COG hit) are resolved by
tax_fix() the same way an unclassified genus is on the taxonomic side -- see
`fix_taxonomy`. Mirrors R's mp_function_barplot()/mp_function_heatmap() in
microplotr -- keep both in sync.
"""
from __future__ import annotations

import pandas as pd

from ..io import MicrobiomeData
from .barplot import taxa_barplot
from .heatmap import taxa_heatmap


def _metadata_stub(feature_table: pd.DataFrame, metadata: pd.DataFrame | None) -> pd.DataFrame:
    if metadata is not None:
        return metadata
    sample_cols = [c for c in feature_table.columns if c != "Feature_ID"]
    return pd.DataFrame({"Sample_ID": sample_cols})


def function_barplot(
    feature_table: pd.DataFrame,
    function_annotation: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
    rank: str = "KEGG_ko",
    group_rank: str | None = "COG_category",
    top_n: int | None = 10,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
    facet_var: str | None = None,
    sample_order: str | list[str] | None = None,
    fix_taxonomy: bool = True,
    nested_legend: bool = True,
    other_label: str = "Other",
    figsize: tuple[float, float] = (7, 5),
):
    data = MicrobiomeData(
        feature_table=feature_table,
        taxonomy=function_annotation,
        metadata=_metadata_stub(feature_table, metadata),
    )
    tax_fix_ranks = [r for r in (group_rank, rank) if r is not None]
    return taxa_barplot(
        data, rank=rank, group_rank=group_rank, top_n=top_n,
        min_rel_abund=min_rel_abund, min_prevalence=min_prevalence, detection=detection,
        facet_var=facet_var, sample_order=sample_order, fix_taxonomy=fix_taxonomy,
        tax_fix_ranks=tax_fix_ranks, nested_legend=nested_legend, other_label=other_label,
        figsize=figsize,
    )


def function_heatmap(
    feature_table: pd.DataFrame,
    function_annotation: pd.DataFrame,
    rank: str = "KEGG_ko",
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
    figsize: tuple[float, float] = (9, 7),
):
    data = MicrobiomeData(
        feature_table=feature_table,
        taxonomy=function_annotation,
        metadata=_metadata_stub(feature_table, None),
    )
    return taxa_heatmap(
        data, rank=rank, top_n=top_n, min_rel_abund=min_rel_abund, min_prevalence=min_prevalence,
        detection=detection, transform=transform, pseudocount=pseudocount, cluster_rows=cluster_rows,
        cluster_cols=cluster_cols, hclust_method=hclust_method, show_dendrogram=show_dendrogram,
        fix_taxonomy=fix_taxonomy, tax_fix_ranks=[rank], figsize=figsize,
    )
