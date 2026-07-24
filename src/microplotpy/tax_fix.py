"""Fix unknown/missing taxonomy labels by falling back to the last known rank.

Mirrors R's mp_tax_fix() in microplotr -- keep both in sync.
"""
from __future__ import annotations

import pandas as pd

ALL_POSSIBLE_RANKS = ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
DEFAULT_UNKNOWN_STRINGS = ("", "na", "unknown", "unclassified", "uncultured", "unidentified", "metagenome")


def tax_fix(
    taxonomy: pd.DataFrame,
    ranks: list[str] | None = None,
    unknown_strings: tuple[str, ...] = DEFAULT_UNKNOWN_STRINGS,
    label_format: str = "Unclassified_{last}",
) -> pd.DataFrame:
    """Replace unknown/missing taxonomy values with the last known ancestor rank.

    Genus-level (or any rank) plots otherwise show blank/NA/"uncultured" bars,
    which is common in 16S taxonomy tables and misleading in a stacked
    barplot. This fills an unknown value at a given rank with the last known
    ancestor rank so every feature still gets a meaningful label -- the same
    idea as `tax_fix()` in the microViz R package (Barnett et al. 2021,
    Bioinformatics). Cascading unknowns (e.g. both Genus and Species missing)
    all fall back to the same true ancestor rather than chaining onto each
    other's fabricated labels.
    """
    if ranks is None:
        ranks = [r for r in ALL_POSSIBLE_RANKS if r in taxonomy.columns]
    if not ranks:
        raise ValueError(
            "No recognized taxonomy rank columns found in `taxonomy`. "
            f"Expected one or more of: {', '.join(ALL_POSSIBLE_RANKS)}"
        )

    unknown_set = {s.lower() for s in unknown_strings}

    def is_unknown(series: pd.Series) -> pd.Series:
        return series.isna() | series.astype(str).str.strip().str.lower().isin(unknown_set)

    out = taxonomy.copy()
    n = len(taxonomy)
    last_good = pd.array([None] * n, dtype="object")

    for rank in ranks:
        values = taxonomy[rank]
        unk = is_unknown(values).to_numpy()
        fixed = values.astype(object).to_numpy(copy=True)

        has_ancestor = unk & pd.notna(last_good)
        no_ancestor = unk & pd.isna(last_good)

        for i in range(n):
            if has_ancestor[i]:
                fixed[i] = label_format.replace("{last}", str(last_good[i]))
            elif no_ancestor[i]:
                fixed[i] = "Unclassified"

        out[rank] = fixed
        # Advance the ancestor pointer only where the real value was known --
        # a fixed/fabricated label never becomes the next rank's ancestor.
        real_values = values.to_numpy(copy=True)
        last_good = pd.array(
            [real_values[i] if not unk[i] else last_good[i] for i in range(n)], dtype="object"
        )

    return out
