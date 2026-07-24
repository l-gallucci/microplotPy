"""Beta diversity ordination plot.

Bray-Curtis or Jaccard dissimilarity (on relative abundance), PCoA (classical
MDS, implemented directly -- no extra dependency) or NMDS (`sklearn.manifold.MDS`,
metric=False) ordination, colored by a grouping metadata variable with 95%
confidence ellipses, PERMANOVA (implemented directly, Anderson 2001 pseudo-F
permutation test) annotated as a plot title. No phylogenetic (UniFrac)
distances -- Bray-Curtis/Jaccard only. Mirrors R's mp_beta_diversity_plot()
in microplotr -- keep both in sync.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2

from ..io import MicrobiomeData
from ..palette import _qualitative_hues
from ..theme import apply_pub_style


def _pcoa(dist: np.ndarray):
    n = dist.shape[0]
    d2 = dist ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    positive_sum = eigvals[eigvals > 0].sum()
    coords = eigvecs[:, :2] * np.sqrt(np.clip(eigvals[:2], 0, None))
    if positive_sum > 0:
        var_explained = eigvals[:2] / positive_sum * 100
    else:
        var_explained = np.zeros(2)  # degenerate (e.g. all-identical) distance matrix
    return coords, var_explained


def _permanova(dist: np.ndarray, groups: np.ndarray, permutations: int = 999, seed: int = 42):
    d2 = dist ** 2
    n = d2.shape[0]
    groups = np.asarray(groups)
    uniq = np.unique(groups)
    a = len(uniq)

    def pseudo_f_and_r2(labels):
        ss_total = d2.sum() / (2 * n)
        ss_within = 0.0
        for g in uniq:
            idx = np.where(labels == g)[0]
            nk = len(idx)
            if nk > 1:
                ss_within += d2[np.ix_(idx, idx)].sum() / (2 * nk)
        ss_between = ss_total - ss_within
        if ss_total <= 0 or ss_within <= 0:
            return 0.0, 0.0  # degenerate (e.g. all-identical) distance matrix
        f_stat = (ss_between / (a - 1)) / (ss_within / (n - a))
        r2 = ss_between / ss_total
        return f_stat, r2

    obs_f, obs_r2 = pseudo_f_and_r2(groups)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(permutations):
        perm_labels = rng.permutation(groups)
        perm_f, _ = pseudo_f_and_r2(perm_labels)
        if perm_f >= obs_f:
            count += 1
    p_value = (count + 1) / (permutations + 1)
    return {"F": obs_f, "R2": obs_r2, "p_value": p_value}


def _draw_confidence_ellipse(ax, x, y, color, level=0.95):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    mean_x, mean_y = np.mean(x), np.mean(y)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    scale = chi2.ppf(level, df=2)
    width, height = 2 * np.sqrt(np.clip(eigvals, 0, None) * scale)
    ellipse = Ellipse((mean_x, mean_y), width, height, angle=angle,
                       facecolor="none", edgecolor=color, linewidth=1.5)
    ax.add_patch(ellipse)


def beta_diversity_plot(
    data: MicrobiomeData,
    group_var: str,
    method: str = "bray",
    ordination: str = "pcoa",
    show_ellipse: bool = True,
    permanova: bool = True,
    permutations: int = 999,
    figsize: tuple[float, float] = (7, 6),
):
    if method not in ("bray", "jaccard"):
        raise ValueError("method must be 'bray' or 'jaccard'")
    if ordination not in ("pcoa", "nmds"):
        raise ValueError("ordination must be 'pcoa' or 'nmds'")

    ft = data.feature_table
    sample_cols = [c for c in ft.columns if c != "Feature_ID"]
    mat = ft[sample_cols].apply(pd.to_numeric).to_numpy().T
    rel = mat / mat.sum(axis=1, keepdims=True) * 100

    metric = "braycurtis" if method == "bray" else "jaccard"
    dist = squareform(pdist(rel, metric=metric))

    meta = data.metadata.set_index("Sample_ID").loc[sample_cols].reset_index()
    groups = meta[group_var].to_numpy()

    stress_label = None
    if ordination == "pcoa":
        coords, var_explained = _pcoa(dist)
        xlab = f"PCoA1 ({var_explained[0]:.1f}%)"
        ylab = f"PCoA2 ({var_explained[1]:.1f}%)"
    else:
        from sklearn.manifold import MDS

        mds = MDS(n_components=2, metric_mds=False, metric="precomputed", random_state=42,
                  normalized_stress=True, n_init=4, init="random")
        coords = mds.fit_transform(dist)
        xlab, ylab = "NMDS1", "NMDS2"
        stress_label = f"stress = {mds.stress_:.3f}"

    fig, ax = plt.subplots(figsize=figsize)
    unique_groups = list(dict.fromkeys(groups))
    palette = dict(zip(unique_groups, _qualitative_hues(len(unique_groups))))
    for g in unique_groups:
        idx = groups == g
        ax.scatter(coords[idx, 0], coords[idx, 1], color=palette[g], label=str(g), alpha=0.85, s=40)
        if show_ellipse:
            _draw_confidence_ellipse(ax, coords[idx, 0], coords[idx, 1], color=palette[g])

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    apply_pub_style(ax)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    caption_parts = []
    if permanova:
        result = _permanova(dist, groups, permutations=permutations)
        p_str = "p < 0.001" if result["p_value"] < 0.001 else f"p = {result['p_value']:.3f}"
        caption_parts.append(f"PERMANOVA R² = {result['R2']:.3f}, {p_str}")
    if stress_label:
        caption_parts.append(stress_label)
    if caption_parts:
        ax.set_title("; ".join(caption_parts), fontsize=10, loc="left")

    fig.tight_layout()
    return fig
