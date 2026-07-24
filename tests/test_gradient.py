from pathlib import Path

from matplotlib.figure import Figure

from microplotpy import load
from microplotpy.plots import taxa_gradient_plot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load():
    d = EXAMPLES / "example_valid"
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


def test_builds_with_nested_grouping_and_smoothing():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum", top_n=8)
    assert isinstance(fig, Figure)
    ax = fig.axes[0]
    assert len(ax.collections) > 0  # scatter points
    assert len(ax.lines) > 0  # smoothing lines


def test_gradient_column_used_as_numeric_x():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum", top_n=8)
    ax = fig.axes[0]
    xdata = ax.collections[0].get_offsets()[:, 0]
    assert xdata.dtype.kind in "fi"


def test_top_n_caps_taxa_legend_shows_other():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum", top_n=3)
    legend = fig.legends[0]
    labels = [t.get_text().strip() for t in legend.get_texts()]
    assert "Other" in labels


def test_group_rank_none_gives_flat_legend_no_bold_headers():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank=None, top_n=5)
    legend = fig.legends[0]
    weights = {t.get_fontweight() for t in legend.get_texts()}
    assert weights == {"normal"} or weights == set()


def test_smooth_false_omits_smoothing_lines():
    data = _load()
    fig_smooth = taxa_gradient_plot(data, gradient_var="Depth_m", top_n=5, smooth=True)
    fig_flat = taxa_gradient_plot(data, gradient_var="Depth_m", top_n=5, smooth=False)
    assert len(fig_smooth.axes[0].lines) > len(fig_flat.axes[0].lines)


def test_facet_true_creates_one_axes_per_taxon_with_visible_title():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum",
                              top_n=5, facet=True)
    used_axes = [ax for ax in fig.axes if ax.get_title()]
    assert len(used_axes) == 6  # 5 taxa + Other
    for ax in used_axes:
        assert ax.get_title() != ""


def test_fix_taxonomy_resolves_missing_genus():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum", top_n=20)
    legend = fig.legends[0]
    labels = [t.get_text().strip() for t in legend.get_texts()]
    assert any(lbl.startswith("Unclassified_") for lbl in labels)


def test_min_rel_abund_filters_taxa_into_other():
    data = _load()
    fig = taxa_gradient_plot(data, gradient_var="Depth_m", rank="Genus", group_rank="Phylum",
                              top_n=None, min_rel_abund=8)
    legend = fig.legends[0]
    labels = [t.get_text().strip() for t in legend.get_texts()]
    assert "Other" in labels
