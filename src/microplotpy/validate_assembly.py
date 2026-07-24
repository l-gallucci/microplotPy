"""Validate assembly QC tables: per-contig lengths and QUAST-shaped summary.

Mirrors R's mp_validate_contig_lengths()/mp_validate_assembly_summary() in
microplotr -- keep both in sync.
"""
from __future__ import annotations

import pandas as pd

from .validate import Finding, ValidationReport

# QUAST's report.tsv uses these column names; map to our canonical names.
_ASSEMBLY_SUMMARY_ALIAS_MAP = {
    "Assembly_ID": ["Assembly_ID", "Assembly"],
    "N_contigs": ["N_contigs", "# contigs"],
    "Largest_contig": ["Largest_contig", "Largest contig"],
    "Total_length": ["Total_length", "Total length"],
    "N50": ["N50"],
    "L50": ["L50"],
    "GC_percent": ["GC_percent", "GC (%)", "GC"],
}


def _normalize(df: pd.DataFrame, alias_map: dict) -> pd.DataFrame:
    out = df.copy()
    for canonical, aliases in alias_map.items():
        if canonical in out.columns:
            continue
        for alias in aliases:
            if alias in out.columns:
                out = out.rename(columns={alias: canonical})
                break
    return out


def normalize_assembly_summary(assembly_summary: pd.DataFrame) -> pd.DataFrame:
    return _normalize(assembly_summary, _ASSEMBLY_SUMMARY_ALIAS_MAP)


def validate_contig_lengths(contig_lengths: pd.DataFrame) -> ValidationReport:
    findings: list[Finding] = []
    required = ["Assembly_ID", "Contig_ID", "Length"]
    for col in required:
        if col not in contig_lengths.columns:
            findings.append(Finding("error", "contig_lengths.tsv", col, f"Required column '{col}' not found."))
    if any(f.level == "error" for f in findings):
        return ValidationReport(findings)

    dup_mask = contig_lengths.duplicated(subset=["Assembly_ID", "Contig_ID"])
    if dup_mask.any():
        dup = contig_lengths.loc[dup_mask, ["Assembly_ID", "Contig_ID"]]
        pairs = [f"{a}/{c}" for a, c in zip(dup["Assembly_ID"], dup["Contig_ID"])]
        findings.append(Finding("error", "contig_lengths.tsv", "Contig_ID",
                                 f"Duplicate Contig_ID within an Assembly_ID: {pairs}."))

    length = pd.to_numeric(contig_lengths["Length"], errors="coerce")
    non_numeric = length.isna() & contig_lengths["Length"].notna()
    if non_numeric.any():
        findings.append(Finding("error", "contig_lengths.tsv", "Length", "Non-numeric value(s) in 'Length'."))
    elif (length <= 0).any():
        findings.append(Finding("error", "contig_lengths.tsv", "Length",
                                 "Value(s) must be positive (contig length in bp)."))

    return ValidationReport(findings)


def validate_assembly_summary(assembly_summary: pd.DataFrame) -> ValidationReport:
    assembly_summary = normalize_assembly_summary(assembly_summary)
    findings: list[Finding] = []
    required = ["Assembly_ID", "N50", "Total_length"]
    for col in required:
        if col not in assembly_summary.columns:
            findings.append(Finding("error", "assembly_summary.tsv", col, f"Required column '{col}' not found."))
    if any(f.level == "error" for f in findings):
        return ValidationReport(findings)

    dup_names = assembly_summary["Assembly_ID"][assembly_summary["Assembly_ID"].duplicated()].unique().tolist()
    if dup_names:
        findings.append(Finding("error", "assembly_summary.tsv", "Assembly_ID",
                                 f"Duplicate Assembly_ID value(s) found: {dup_names}."))

    for col in ("N50", "Total_length"):
        vals = pd.to_numeric(assembly_summary[col], errors="coerce")
        non_numeric = vals.isna() & assembly_summary[col].notna()
        if non_numeric.any():
            findings.append(Finding("error", "assembly_summary.tsv", col, f"Non-numeric value(s) in '{col}'."))
        elif (vals <= 0).any():
            findings.append(Finding("error", "assembly_summary.tsv", col, f"Value(s) in '{col}' must be positive."))

    return ValidationReport(findings)
