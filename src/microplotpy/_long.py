"""Melt feature table into long format, joined with taxonomy and metadata.

Internal helper shared by the taxa-level plot functions.
"""
from __future__ import annotations

import pandas as pd

from .io import MicrobiomeData
from .tax_fix import tax_fix


def long_abundance(data: MicrobiomeData, fix_taxonomy: bool = True, tax_fix_ranks: list[str] | None = None) -> pd.DataFrame:
    """`tax_fix_ranks` is passed through to `tax_fix()`'s `ranks` argument.
    `None` (default) auto-detects standard 16S rank columns; pass explicitly
    when `data.taxonomy` holds a non-taxonomic hierarchy (e.g. a functional
    annotation table) whose column names aren't in the fixed 16S rank set.
    """
    taxonomy = tax_fix(data.taxonomy, ranks=tax_fix_ranks) if fix_taxonomy else data.taxonomy
    ft = data.feature_table
    sample_cols = [c for c in ft.columns if c != "Feature_ID"]

    long = ft.melt(id_vars="Feature_ID", value_vars=sample_cols, var_name="Sample_ID", value_name="Count")
    long["Count"] = pd.to_numeric(long["Count"])
    long = long.merge(taxonomy, on="Feature_ID", how="left")
    long = long.merge(data.metadata, on="Sample_ID", how="left")
    return long
