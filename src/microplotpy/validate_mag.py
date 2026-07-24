"""Validate a MAG (metagenome-assembled genome) quality table.

Checks a single flat table shaped like CheckM2's `quality_report.tsv` or
CheckM (v1)'s `qa`/`bin_stats` table (one row per bin/MAG) -- not the 3-file
taxonomy/metadata schema used by `validate()`. Mirrors R's
mp_validate_mag()/mp_mag_quality_plot() in microplotr -- keep both in sync.
"""
from __future__ import annotations

import pandas as pd

from .validate import Finding, ValidationReport

# CheckM1's qa/bin_stats table uses "Bin Id"; CheckM2's quality_report.tsv
# uses "Name". Both tools agree on "Completeness"/"Contamination" already.
_MAG_ALIAS_MAP = {
    "Name": ["Name", "Bin Id", "Bin_Id", "bin_id"],
    "Genome_Size": ["Genome_Size", "Genome size (bp)", "Genome_size"],
    "Contig_N50": ["Contig_N50", "N50 (contigs)", "N50"],
    "GC_Content": ["GC_Content", "GC", "GC (%)"],
}


def normalize_mag_table(mag_table: pd.DataFrame) -> pd.DataFrame:
    """Rename whichever alias of each known MAG-table column is present to
    its canonical name (CheckM1 vs CheckM2 column-naming differences)."""
    out = mag_table.copy()
    for canonical, aliases in _MAG_ALIAS_MAP.items():
        if canonical in out.columns:
            continue
        for alias in aliases:
            if alias in out.columns:
                out = out.rename(columns={alias: canonical})
                break
    return out


def validate_mag(mag_table: pd.DataFrame) -> ValidationReport:
    mag_table = normalize_mag_table(mag_table)
    findings: list[Finding] = []

    required = ["Name", "Completeness", "Contamination", "Genome_Size"]
    for col in required:
        if col not in mag_table.columns:
            findings.append(Finding("error", "mag_quality.tsv", col, f"Required column '{col}' not found."))
    if any(f.level == "error" for f in findings):
        return ValidationReport(findings)

    dup_names = mag_table["Name"][mag_table["Name"].duplicated()].unique().tolist()
    if dup_names:
        findings.append(Finding("error", "mag_quality.tsv", "Name",
                                 f"Duplicate MAG Name value(s) found: {dup_names}."))

    completeness = pd.to_numeric(mag_table["Completeness"], errors="coerce")
    non_numeric = completeness.isna() & mag_table["Completeness"].notna()
    if non_numeric.any():
        findings.append(Finding("error", "mag_quality.tsv", "Completeness",
                                 "Non-numeric value(s) in 'Completeness'."))
    elif ((completeness < 0) | (completeness > 100)).any():
        findings.append(Finding("error", "mag_quality.tsv", "Completeness",
                                 "Value(s) outside the valid [0, 100] range."))

    contamination = pd.to_numeric(mag_table["Contamination"], errors="coerce")
    non_numeric = contamination.isna() & mag_table["Contamination"].notna()
    if non_numeric.any():
        findings.append(Finding("error", "mag_quality.tsv", "Contamination",
                                 "Non-numeric value(s) in 'Contamination'."))
    elif (contamination < 0).any():
        findings.append(Finding("error", "mag_quality.tsv", "Contamination",
                                 "Negative value(s) found; contamination cannot be negative."))

    genome_size = pd.to_numeric(mag_table["Genome_Size"], errors="coerce")
    non_numeric = genome_size.isna() & mag_table["Genome_Size"].notna()
    if non_numeric.any():
        findings.append(Finding("error", "mag_quality.tsv", "Genome_Size",
                                 "Non-numeric value(s) in 'Genome_Size'."))
    elif (genome_size <= 0).any():
        findings.append(Finding("error", "mag_quality.tsv", "Genome_Size",
                                 "Value(s) must be positive (genome size in bp)."))

    return ValidationReport(findings)
