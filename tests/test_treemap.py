from pathlib import Path

import pandas as pd
from matplotlib.figure import Figure

from microplotpy import load
from microplotpy.plots import function_treemap, taxa_treemap

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load_taxa():
    d = EXAMPLES / "example_valid"
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


def _load_function():
    ft = pd.read_csv(EXAMPLES / "example_function_profile" / "gene_count_table.tsv", sep="\t")
    fa = pd.read_csv(EXAMPLES / "example_function_profile" / "function_annotation.tsv", sep="\t")
    return ft, fa


def test_taxa_treemap_builds_with_grouping():
    data = _load_taxa()
    fig = taxa_treemap(data, rank="Genus", group_rank="Phylum")
    assert isinstance(fig, Figure)


def test_taxa_treemap_areas_sum_to_100():
    data = _load_taxa()
    fig = taxa_treemap(data, rank="Genus", group_rank="Phylum")
    rects = [p for p in fig.axes[0].patches if p.get_facecolor()[3] > 0]  # filled rects only (skip outer borders)
    total_area = sum(p.get_width() * p.get_height() for p in rects)
    assert abs(total_area - 100 * 100) < 1e-3


def test_taxa_treemap_works_without_group_rank():
    data = _load_taxa()
    fig = taxa_treemap(data, rank="Genus", group_rank=None, top_n=8)
    assert isinstance(fig, Figure)


def test_top_n_pools_remainder_into_other_label():
    data = _load_taxa()
    fig = taxa_treemap(data, rank="Genus", group_rank="Phylum", top_n=3)
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert "Other" in texts


def test_fix_taxonomy_resolves_missing_genus():
    data = _load_taxa()
    fig = taxa_treemap(data, rank="Genus", group_rank="Phylum", top_n=None, fix_taxonomy=True)
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert any(t.startswith("Unclassified_") for t in texts)


def test_function_treemap_builds_on_eggnog_shaped_data():
    ft, fa = _load_function()
    fig = function_treemap(ft, fa, rank="KEGG_ko", group_rank="COG_category")
    assert isinstance(fig, Figure)


def test_all_labels_clipped_to_own_rectangle():
    data = _load_taxa()
    fig = taxa_treemap(data, rank="Genus", group_rank="Phylum")
    for text in fig.axes[0].texts:
        assert text.get_clip_path() is not None
