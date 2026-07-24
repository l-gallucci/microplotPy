from pathlib import Path

from matplotlib.figure import Figure

from microplotpy import load
from microplotpy.plots import alpha_diversity_plot, beta_diversity_plot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load():
    d = EXAMPLES / "example_valid"
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


# --- alpha diversity ---

def test_alpha_diversity_builds_one_axes_per_metric():
    data = _load()
    fig = alpha_diversity_plot(data, group_var="Group")
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 3


def test_observed_richness_no_fake_jitter_stays_constant():
    data = _load()
    fig = alpha_diversity_plot(data, group_var="Group", metrics=["Observed"])
    ax = fig.axes[0]
    ydata = [pt[1] for coll in ax.collections for pt in coll.get_offsets()]
    assert len(set(ydata)) == 1


def test_degenerate_metric_omits_p_value_text_instead_of_nan():
    data = _load()
    fig = alpha_diversity_plot(data, group_var="Group", metrics=["Observed"])
    ax = fig.axes[0]
    texts = [t.get_text() for t in ax.texts]
    assert not any("nan" in t for t in texts)


def test_metrics_subset_respected():
    data = _load()
    fig = alpha_diversity_plot(data, group_var="Group", metrics=["Shannon", "Simpson"])
    assert len(fig.axes) == 2
    titles = [ax.get_title() for ax in fig.axes]
    assert titles == ["Shannon", "Simpson"]


def test_test_none_omits_p_value_annotations():
    data = _load()
    fig = alpha_diversity_plot(data, group_var="Group", metrics=["Shannon"], test=None)
    ax = fig.axes[0]
    assert len(ax.texts) == 0


# --- beta diversity ---

def test_beta_diversity_pcoa_builds_with_permanova_title():
    data = _load()
    fig = beta_diversity_plot(data, group_var="Group", method="bray", ordination="pcoa")
    assert isinstance(fig, Figure)
    assert "PERMANOVA" in fig.axes[0].get_title(loc="left")


def test_beta_diversity_nmds_builds_with_stress_annotated():
    data = _load()
    fig = beta_diversity_plot(data, group_var="Group", method="bray", ordination="nmds")
    assert "stress" in fig.axes[0].get_title(loc="left")


def test_jaccard_method_runs_without_error():
    data = _load()
    fig = beta_diversity_plot(data, group_var="Group", method="jaccard", ordination="pcoa")
    assert isinstance(fig, Figure)


def test_show_ellipse_false_omits_ellipse_patches():
    data = _load()
    fig_ellipse = beta_diversity_plot(data, group_var="Group", show_ellipse=True, permanova=False)
    fig_flat = beta_diversity_plot(data, group_var="Group", show_ellipse=False, permanova=False)
    assert len(fig_ellipse.axes[0].patches) > len(fig_flat.axes[0].patches)


def test_permanova_false_and_ellipse_false_gives_no_title():
    data = _load()
    fig = beta_diversity_plot(data, group_var="Group", permanova=False, show_ellipse=False)
    assert fig.axes[0].get_title() == ""


def test_invalid_method_raises():
    data = _load()
    try:
        beta_diversity_plot(data, group_var="Group", method="euclidean")
        assert False, "expected ValueError"
    except ValueError:
        pass
