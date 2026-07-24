from pathlib import Path

import pandas as pd
from matplotlib.figure import Figure

from microplotpy import MicrobiomeData
from microplotpy.plots import function_barplot, function_heatmap, taxa_barplot, taxa_bubbleplot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str, file: str) -> pd.DataFrame:
    return pd.read_csv(EXAMPLES / name / file, sep="\t")


def test_function_barplot_builds_with_cog_grouping():
    ft = _load("example_function_profile", "gene_count_table.tsv")
    fa = _load("example_function_profile", "function_annotation.tsv")
    fig = function_barplot(ft, fa, rank="KEGG_ko", group_rank="COG_category", top_n=8)
    assert isinstance(fig, Figure)


def test_function_barplot_resolves_blank_cog_category():
    ft = _load("example_function_profile", "gene_count_table.tsv")
    fa = _load("example_function_profile", "function_annotation.tsv")
    assert fa["COG_category"].isna().any()  # sanity: raw data has blanks (read as NaN)
    fig = function_barplot(ft, fa, rank="KEGG_ko", group_rank="COG_category", top_n=20)
    assert isinstance(fig, Figure)


def test_function_heatmap_builds():
    ft = _load("example_function_profile", "gene_count_table.tsv")
    fa = _load("example_function_profile", "function_annotation.tsv")
    fig = function_heatmap(ft, fa, rank="KEGG_ko", top_n=20, show_dendrogram=False)
    assert isinstance(fig, Figure)


def test_flat_kofamscan_style_annotation_works():
    ft = _load("example_function_profile", "gene_count_table.tsv")
    fa = _load("example_function_profile", "function_annotation.tsv")
    fa_kofam = fa[["Feature_ID", "KEGG_ko", "Description"]].rename(
        columns={"KEGG_ko": "KO", "Description": "KO_definition"}
    )
    fig = function_barplot(ft, fa_kofam, rank="KO", group_rank=None, nested_legend=False, top_n=8)
    assert isinstance(fig, Figure)


def test_regression_rank_mapping_to_two_groups_does_not_crash_barplot():
    ft = pd.DataFrame({"Feature_ID": ["g1", "g2", "g3"], "S1": [10, 20, 5], "S2": [15, 5, 10]})
    fa = pd.DataFrame({
        "Feature_ID": ["g1", "g2", "g3"],
        "Category": ["A", "A", "B"],
        "Item": ["X", "X", "X"],  # "X" maps to both "A" and "B" across rows
    })
    data = MicrobiomeData(feature_table=ft, taxonomy=fa, metadata=pd.DataFrame({"Sample_ID": ["S1", "S2"]}))
    fig = taxa_barplot(data, rank="Item", group_rank="Category", top_n=None, fix_taxonomy=False)
    assert isinstance(fig, Figure)


def test_regression_rank_mapping_to_two_groups_does_not_crash_bubbleplot():
    ft = pd.DataFrame({"Feature_ID": ["g1", "g2", "g3"], "S1": [10, 20, 5], "S2": [15, 5, 10]})
    fa = pd.DataFrame({
        "Feature_ID": ["g1", "g2", "g3"],
        "Category": ["A", "A", "B"],
        "Item": ["X", "X", "X"],
    })
    data = MicrobiomeData(feature_table=ft, taxonomy=fa, metadata=pd.DataFrame({"Sample_ID": ["S1", "S2"]}))
    fig = taxa_bubbleplot(data, rank="Item", group_rank="Category", top_n=None, fix_taxonomy=False)
    assert isinstance(fig, Figure)
