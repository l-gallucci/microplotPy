from pathlib import Path

import pytest
from matplotlib.figure import Figure

from microplotpy import load
from microplotpy.plots import taxa_barplot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str):
    d = EXAMPLES / name
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


def test_nested_legend_barplot_builds_without_error():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=5)
    assert isinstance(fig, Figure)


def test_top_n_caps_number_of_taxa_shown_plus_other():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=3)
    ax = fig.axes[0]
    n_samples = 8
    assert len(ax.patches) == (3 + 1) * n_samples  # 3 kept taxa + Other, one bar segment per sample


def test_nested_legend_true_without_group_rank_raises():
    data = _load("example_valid")
    with pytest.raises(ValueError, match="requires group_rank"):
        taxa_barplot(data, rank="Genus", group_rank=None, nested_legend=True)


def test_nested_legend_false_gives_flat_legend_no_bold_headers():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank=None, nested_legend=False, top_n=5)
    legend = fig.legends[0]
    weights = {t.get_fontweight() for t in legend.get_texts()}
    assert weights == {"normal"} or weights == set()


def test_fix_taxonomy_resolves_missing_genus():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=20, fix_taxonomy=True)
    legend = fig.legends[0]
    labels = [t.get_text() for t in legend.get_texts()]
    assert any(lbl.strip().startswith("Unclassified_") for lbl in labels)


def test_facet_var_creates_one_axes_per_level():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=5, facet_var="Group")
    assert len(fig.axes) == 2  # Group has 2 levels (A, B)


def test_explicit_sample_order_sets_xtick_labels():
    data = _load("example_valid")
    custom = ["S8", "S7", "S6", "S5", "S4", "S3", "S2", "S1"]
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=5, sample_order=custom)
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels == custom


def test_min_rel_abund_pools_low_abundance_taxa_into_other_even_with_top_n_none():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=None, min_rel_abund=8)
    legend = fig.legends[0]
    labels = [t.get_text().strip() for t in legend.get_texts()]
    assert "Other" in labels
    assert len(labels) < 20  # fewer entries than the full 10 taxa + headers


def test_top_n_none_with_no_thresholds_shows_every_taxon_no_other():
    data = _load("example_valid")
    fig = taxa_barplot(data, rank="Genus", group_rank="Phylum", top_n=None)
    legend = fig.legends[0]
    labels = [t.get_text().strip() for t in legend.get_texts()]
    assert "Other" not in labels
