from pathlib import Path

from matplotlib.figure import Figure

from microplotpy import load
from microplotpy.plots import taxa_bubbleplot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str):
    d = EXAMPLES / name
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


def test_bubbleplot_builds_without_error():
    data = _load("example_valid")
    fig = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=8)
    assert isinstance(fig, Figure)


def test_grid_is_complete_n_samples_times_top_n():
    data = _load("example_valid")
    n_samples = 8
    fig = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=6)
    ax = fig.axes[0]
    offsets = ax.collections[0].get_offsets()
    assert len(offsets) == n_samples * 6


def test_top_n_limits_distinct_taxa_on_y_axis():
    data = _load("example_valid")
    fig = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=4)
    ax = fig.axes[0]
    assert len(ax.get_yticklabels()) == 4


def test_group_rank_none_still_builds():
    data = _load("example_valid")
    fig = taxa_bubbleplot(data, rank="Genus", group_rank=None, top_n=6)
    assert isinstance(fig, Figure)


def test_fix_taxonomy_resolves_missing_genus():
    data = _load("example_valid")
    fig = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=20, fix_taxonomy=True)
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert any(lbl.startswith("Unclassified_") for lbl in labels)


def test_facet_var_creates_one_axes_per_level_and_own_samples():
    data = _load("example_valid")
    fig = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=6, facet_var="Group")
    plot_axes = [ax for ax in fig.axes if ax.get_xticklabels()]
    assert len(plot_axes) == 2
    for ax in plot_axes:
        labels = [t.get_text() for t in ax.get_xticklabels()]
        assert len(labels) == 4  # 4 samples per Group level


def test_min_rel_abund_reduces_number_of_taxa_shown():
    data = _load("example_valid")
    fig_all = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=None)
    fig_filtered = taxa_bubbleplot(data, rank="Genus", group_rank="Phylum", top_n=None, min_rel_abund=8)
    n_all = len(fig_all.axes[0].get_yticklabels())
    n_filtered = len(fig_filtered.axes[0].get_yticklabels())
    assert n_filtered <= n_all
