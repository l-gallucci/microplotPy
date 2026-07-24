"""Alpha diversity boxplot.

Observed richness, Shannon, and Simpson diversity per sample, boxplot +
jitter by a grouping metadata variable, faceted one panel per metric, with a
Mann-Whitney U (2 groups) or Kruskal-Wallis (>2 groups) test annotated.
Mirrors R's mp_alpha_diversity_plot() in microplotr -- keep both in sync.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from ..io import MicrobiomeData
from ..theme import apply_pub_style


def _diversity_table(feature_table: pd.DataFrame, metrics) -> pd.DataFrame:
    sample_cols = [c for c in feature_table.columns if c != "Feature_ID"]
    mat = feature_table[sample_cols].apply(pd.to_numeric).to_numpy().T  # samples x features

    out = {"Sample_ID": sample_cols}
    if "Observed" in metrics:
        out["Observed"] = (mat > 0).sum(axis=1)

    row_sums = mat.sum(axis=1, keepdims=True)
    props = np.divide(mat, row_sums, out=np.zeros_like(mat, dtype=float), where=row_sums > 0)

    if "Shannon" in metrics:
        with np.errstate(divide="ignore", invalid="ignore"):
            plogp = np.where(props > 0, props * np.log(props), 0.0)
        out["Shannon"] = -plogp.sum(axis=1)
    if "Simpson" in metrics:
        out["Simpson"] = 1 - (props ** 2).sum(axis=1)

    return pd.DataFrame(out)


def alpha_diversity_plot(
    data: MicrobiomeData,
    group_var: str,
    metrics=("Observed", "Shannon", "Simpson"),
    test: str | None = "auto",
    figsize: tuple[float, float] = (9, 4),
):
    metrics_present = [m for m in metrics if m in ("Observed", "Shannon", "Simpson")]
    df = _diversity_table(data.feature_table, metrics_present)
    long = df.melt(id_vars="Sample_ID", value_vars=metrics_present, var_name="metric", value_name="value")
    long = long.merge(data.metadata[["Sample_ID", group_var]], on="Sample_ID", how="left")

    groups = list(dict.fromkeys(long[group_var]))
    n_groups = len(groups)
    method = test
    if test == "auto":
        method = "wilcoxon" if n_groups == 2 else "kruskal"

    fig, axes_arr = plt.subplots(1, len(metrics_present), figsize=figsize, squeeze=False)
    axes = axes_arr[0]
    rng = np.random.default_rng(0)

    for ax, metric in zip(axes, metrics_present):
        sub = long[long["metric"] == metric]
        group_vals = [sub.loc[sub[group_var] == g, "value"].to_numpy() for g in groups]

        ax.boxplot(group_vals, tick_labels=groups, showfliers=False)
        for i, vals in enumerate(group_vals, start=1):
            x = rng.normal(i, 0.05, size=len(vals))
            ax.scatter(x, vals, alpha=0.7, s=15, color="black", zorder=3)

        apply_pub_style(ax)
        ax.set_title(metric, fontweight="bold", fontsize=10)

        if method is not None:
            try:
                if method == "wilcoxon" and n_groups == 2:
                    _, p = stats.mannwhitneyu(group_vals[0], group_vals[1])
                else:
                    _, p = stats.kruskal(*group_vals)
            except ValueError:
                p = np.nan  # e.g. all values identical/tied -- nothing meaningful to test
            if not np.isnan(p):
                p_str = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
                ax.text(0.5, 1.12, p_str, transform=ax.transAxes, ha="center", fontsize=9)

    axes[0].set_ylabel("Diversity value")
    fig.tight_layout()
    return fig
