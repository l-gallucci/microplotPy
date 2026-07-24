"""Loading tidy microbiome input files (feature table, taxonomy, metadata)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class MicrobiomeData:
    feature_table: pd.DataFrame  # index = Feature_ID, columns = Sample_ID
    taxonomy: pd.DataFrame       # index = Feature_ID, columns = ranks
    metadata: pd.DataFrame       # index = Sample_ID, columns = variables


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    return pd.read_csv(path, sep=sep, dtype=str)


def load(feature_table_path, taxonomy_path, metadata_path) -> MicrobiomeData:
    """Read the three raw input files as string-typed DataFrames (no validation, no coercion).

    Values are kept as strings here so the validator can distinguish
    "not numeric" from "numeric" without pandas silently coercing bad
    values to NaN beforehand.
    """
    feature_table = _read_table(feature_table_path)
    taxonomy = _read_table(taxonomy_path)
    metadata = _read_table(metadata_path)
    return MicrobiomeData(feature_table=feature_table, taxonomy=taxonomy, metadata=metadata)
