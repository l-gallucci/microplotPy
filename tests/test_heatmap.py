from pathlib import Path

from matplotlib.figure import Figure

from microplotpy import load
from microplotpy.plots import taxa_heatmap

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str):
    d = EXAMPLES / name
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


def test_default_heatmap_with_dendrograms_has_three_axes():
    data = _load("example_valid")
    fig = taxa_heatmap(data, rank="Genus", top_n=10)
    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 3  # top dendro, left dendro, heatmap (+ colorbar axes)


def test_show_dendrogram_false_has_single_axes_plus_colorbar():
    data = _load("example_valid")
    fig = taxa_heatmap(data, rank="Genus", top_n=10, show_dendrogram=False)
    assert len(fig.axes) == 2  # heatmap + colorbar


def test_top_n_limits_rows_shown():
    data = _load("example_valid")
    fig = taxa_heatmap(data, rank="Genus", top_n=5, show_dendrogram=False)
    ax = fig.axes[0]
    assert len(ax.get_yticklabels()) == 5


def test_log10_and_clr_give_different_legend_titles():
    data = _load("example_valid")
    fig_clr = taxa_heatmap(data, rank="Genus", top_n=8, transform="clr", show_dendrogram=False)
    fig_log = taxa_heatmap(data, rank="Genus", top_n=8, transform="log10", show_dendrogram=False)
    assert fig_clr.axes[-1].get_ylabel() != fig_log.axes[-1].get_ylabel()


def test_cluster_flags_false_skip_clustering_without_error():
    data = _load("example_valid")
    fig = taxa_heatmap(data, rank="Genus", top_n=8, cluster_rows=False, cluster_cols=False, show_dendrogram=False)
    assert isinstance(fig, Figure)


def test_fix_taxonomy_resolves_missing_genus():
    data = _load("example_valid")
    fig = taxa_heatmap(data, rank="Genus", top_n=20, show_dendrogram=False, fix_taxonomy=True)
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert any(lbl.startswith("Unclassified_") for lbl in labels)


def test_min_prevalence_excludes_taxa_detected_in_too_few_samples():
    data = _load("example_valid")
    fig_all = taxa_heatmap(data, rank="Genus", top_n=None, show_dendrogram=False)
    fig_filtered = taxa_heatmap(data, rank="Genus", top_n=None, min_prevalence=0.99, show_dendrogram=False)
    n_all = len(fig_all.axes[0].get_yticklabels())
    n_filtered = len(fig_filtered.axes[0].get_yticklabels())
    assert n_filtered <= n_all
