"""Filter taxa by minimum abundance and/or prevalence.

Mirrors R's mp_filter_taxa() in microplotr -- keep both in sync. Shared by
taxa_barplot(), taxa_heatmap(), taxa_bubbleplot(): narrows the candidate
taxa *before* top_n is applied, so a taxon is excluded from the "top" set
either because it failed a threshold here or because it wasn't in the top N
of what passed -- both cases are treated identically by the calling plot
function (pooled into other_label for the barplot, simply dropped for the
heatmap/bubbleplot).
"""
from __future__ import annotations

import pandas as pd


def filter_taxa(
    agg: pd.DataFrame,
    taxon_col: str,
    min_rel_abund: float = 0,
    min_prevalence: float = 0,
    detection: float = 0,
) -> list:
    """Return the taxa in `agg[taxon_col]` passing both filters.

    `agg` must have one row per (Sample_ID, taxon) with a `rel_abund` column
    (0 where absent) so prevalence/mean are computed correctly.

    `min_rel_abund`: minimum mean relative abundance (%) across all samples.
    `min_prevalence`: minimum fraction of samples (0-1), or if > 1 an
    absolute sample count, in which the taxon must be "detected" (rel_abund
    > `detection`).
    """
    n_samples = agg["Sample_ID"].nunique()
    prevalence_threshold = min_prevalence if min_prevalence > 1 else min_prevalence * n_samples

    summary = agg.groupby(taxon_col).agg(
        mean_rel_abund=("rel_abund", "mean"),
        prevalence=("rel_abund", lambda s: (s > detection).sum()),
    )

    keep = summary[(summary["mean_rel_abund"] >= min_rel_abund) & (summary["prevalence"] >= prevalence_threshold)]
    return list(keep.index)
