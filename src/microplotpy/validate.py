"""Validate tidy microbiome input (feature table + taxonomy + metadata).

Mirrors the checks documented in data-format.md and implemented on the R
side by mp_validate() in microplotr — keep both in sync.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import pandas as pd

from .io import MicrobiomeData

REQUIRED_TAXONOMY_RANKS = ("Phylum", "Genus")

Level = Literal["error", "warning"]


@dataclass
class Finding:
    level: Level
    file: str
    field: str
    message: str


@dataclass
class ValidationReport:
    findings: list[Finding] = field(default_factory=list)

    @property
    def errors(self) -> list[Finding]:
        return [f for f in self.findings if f.level == "error"]

    @property
    def warnings(self) -> list[Finding]:
        return [f for f in self.findings if f.level == "warning"]

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0

    def __bool__(self) -> bool:
        return self.is_valid


def _to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def validate(
    data: MicrobiomeData,
    gradient_column: str | None = None,
    group_column: str | None = None,
    required_ranks: tuple[str, ...] = REQUIRED_TAXONOMY_RANKS,
) -> ValidationReport:
    """`required_ranks`: column names that must be present in `data.taxonomy`.
    Defaults to the standard 16S taxonomy shape (`Phylum`, `Genus`). Pass
    `()` when `data.taxonomy` is actually a non-taxonomic hierarchy (e.g. a
    functional annotation table from eggNOG-mapper/KofamScan, see
    `function_barplot()`) whose column names aren't ranks.
    """
    findings: list[Finding] = []
    ft, tax, meta = data.feature_table, data.taxonomy, data.metadata

    # --- required columns present ---
    if "Feature_ID" not in ft.columns:
        findings.append(Finding("error", "feature_table.tsv", "Feature_ID",
                                 "Required column 'Feature_ID' not found."))
    if "Feature_ID" not in tax.columns:
        findings.append(Finding("error", "taxonomy.tsv", "Feature_ID",
                                 "Required column 'Feature_ID' not found."))
    if "Sample_ID" not in meta.columns:
        findings.append(Finding("error", "metadata.tsv", "Sample_ID",
                                 "Required column 'Sample_ID' not found."))
    for rank in required_ranks:
        if rank not in tax.columns:
            findings.append(Finding("error", "taxonomy.tsv", rank,
                                     f"Required taxonomy rank column '{rank}' not found."))

    # bail out early if the basic skeleton isn't there — everything below assumes it exists
    if any(f.level == "error" for f in findings):
        return ValidationReport(findings)

    # --- uniqueness ---
    dup_features_ft = ft["Feature_ID"][ft["Feature_ID"].duplicated()].unique().tolist()
    if dup_features_ft:
        findings.append(Finding("error", "feature_table.tsv", "Feature_ID",
                                 f"Duplicate Feature_ID values found: {dup_features_ft}."))
    dup_features_tax = tax["Feature_ID"][tax["Feature_ID"].duplicated()].unique().tolist()
    if dup_features_tax:
        findings.append(Finding("error", "taxonomy.tsv", "Feature_ID",
                                 f"Duplicate Feature_ID values found: {dup_features_tax}."))
    dup_samples = meta["Sample_ID"][meta["Sample_ID"].duplicated()].unique().tolist()
    if dup_samples:
        findings.append(Finding("error", "metadata.tsv", "Sample_ID",
                                 f"Duplicate Sample_ID values found: {dup_samples}."))

    # --- sample ID match: feature_table columns vs metadata Sample_ID ---
    ft_samples = set(ft.columns) - {"Feature_ID"}
    meta_samples = set(meta["Sample_ID"].astype(str))
    missing_in_ft = sorted(meta_samples - ft_samples)
    missing_in_meta = sorted(ft_samples - meta_samples)
    if missing_in_ft:
        findings.append(Finding("error", "metadata.tsv", "Sample_ID",
                                 f"Sample(s) {missing_in_ft} found in metadata.tsv but missing from "
                                 f"feature_table.tsv columns."))
    if missing_in_meta:
        findings.append(Finding("error", "feature_table.tsv", "Sample_ID",
                                 f"Sample column(s) {missing_in_meta} found in feature_table.tsv but "
                                 f"missing from metadata.tsv Sample_ID."))

    # --- feature ID match: feature_table vs taxonomy ---
    ft_features = set(ft["Feature_ID"].astype(str))
    tax_features = set(tax["Feature_ID"].astype(str))
    missing_in_tax = sorted(ft_features - tax_features)
    missing_in_ft_feat = sorted(tax_features - ft_features)
    if missing_in_tax:
        findings.append(Finding("error", "taxonomy.tsv", "Feature_ID",
                                 f"Feature(s) {missing_in_tax} found in feature_table.tsv but missing "
                                 f"from taxonomy.tsv."))
    if missing_in_ft_feat:
        findings.append(Finding("error", "feature_table.tsv", "Feature_ID",
                                 f"Feature(s) {missing_in_ft_feat} found in taxonomy.tsv but missing "
                                 f"from feature_table.tsv."))

    # --- numeric, non-negative abundance values ---
    sample_cols = [c for c in ft.columns if c != "Feature_ID"]
    numeric_ft = ft[sample_cols].apply(_to_numeric)
    non_numeric_mask = numeric_ft.isna() & ft[sample_cols].notna()
    if non_numeric_mask.to_numpy().any():
        bad_cells = [f"{ft.loc[r, 'Feature_ID']}/{c}" for r in non_numeric_mask.index
                     for c in sample_cols if non_numeric_mask.loc[r, c]]
        findings.append(Finding("error", "feature_table.tsv", "values",
                                 f"Non-numeric abundance value(s) at: {bad_cells[:10]}"
                                 + (" (truncated)" if len(bad_cells) > 10 else "") + "."))
    else:
        negative_mask = numeric_ft < 0
        if negative_mask.to_numpy().any():
            bad_cells = [f"{ft.loc[r, 'Feature_ID']}/{c}" for r in negative_mask.index
                         for c in sample_cols if negative_mask.loc[r, c]]
            findings.append(Finding("error", "feature_table.tsv", "values",
                                     f"Negative abundance value(s) at: {bad_cells[:10]}"
                                     + (" (truncated)" if len(bad_cells) > 10 else "") + "."))

        # --- all-zero rows / columns (warning only) ---
        zero_features = ft.loc[(numeric_ft.fillna(0).sum(axis=1) == 0), "Feature_ID"].tolist()
        if zero_features:
            findings.append(Finding("warning", "feature_table.tsv", "Feature_ID",
                                     f"Feature(s) {zero_features} are all-zero across every sample."))
        zero_samples = [c for c in sample_cols if numeric_ft[c].fillna(0).sum() == 0]
        if zero_samples:
            findings.append(Finding("warning", "feature_table.tsv", "Sample_ID",
                                     f"Sample(s) {zero_samples} are all-zero across every feature."))

    # --- missing taxonomy values (warning, treated as Unclassified) ---
    for rank in tax.columns:
        if rank == "Feature_ID":
            continue
        missing_mask = tax[rank].isna() | (tax[rank].astype(str).str.strip() == "")
        if missing_mask.any():
            feats = tax.loc[missing_mask, "Feature_ID"].tolist()
            findings.append(Finding("warning", "taxonomy.tsv", rank,
                                     f"Missing '{rank}' for feature(s) {feats}; will be labeled 'Unclassified'."))

    # --- gradient column check ---
    if gradient_column is not None:
        if gradient_column not in meta.columns:
            findings.append(Finding("error", "metadata.tsv", gradient_column,
                                     f"Gradient column '{gradient_column}' not found in metadata.tsv."))
        else:
            coerced = _to_numeric(meta[gradient_column])
            if coerced.isna().any() and meta[gradient_column].notna().any():
                findings.append(Finding("error", "metadata.tsv", gradient_column,
                                         f"Gradient column '{gradient_column}' contains non-numeric value(s)."))

    # --- group column check ---
    if group_column is not None:
        if group_column not in meta.columns:
            findings.append(Finding("error", "metadata.tsv", group_column,
                                     f"Group column '{group_column}' not found in metadata.tsv."))
        else:
            levels = meta[group_column].dropna().unique().tolist()
            if len(levels) < 2:
                findings.append(Finding("error", "metadata.tsv", group_column,
                                         f"Group column '{group_column}' has {len(levels)} level(s); "
                                         f"need at least 2 for group comparison."))

    return ValidationReport(findings)
