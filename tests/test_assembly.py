from pathlib import Path

import pandas as pd
import pytest
from matplotlib.figure import Figure

from microplotpy import validate_assembly_summary, validate_contig_lengths
from microplotpy.plots import assembly_nx_plot, assembly_summary_barplot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str, file: str) -> pd.DataFrame:
    return pd.read_csv(EXAMPLES / name / file, sep="\t")


def test_valid_contig_lengths_passes():
    contigs = _load("example_assembly", "contig_lengths.tsv")
    assert validate_contig_lengths(contigs).is_valid


def test_valid_assembly_summary_passes():
    summary_tbl = _load("example_assembly", "assembly_summary.tsv")
    assert validate_assembly_summary(summary_tbl).is_valid


def test_negative_contig_length_detected():
    bad = _load("example_assembly_broken_negative_length", "contig_lengths.tsv")
    report = validate_contig_lengths(bad)
    assert not report.is_valid
    assert any("must be positive" in f.message for f in report.errors)


def test_duplicate_contig_id_within_assembly_detected():
    bad = _load("example_assembly_broken_duplicate_contig_ids", "contig_lengths.tsv")
    report = validate_contig_lengths(bad)
    assert not report.is_valid
    assert any("Duplicate Contig_ID" in f.message for f in report.errors)


def test_missing_required_summary_column_detected():
    bad = _load("example_assembly_broken_missing_columns", "assembly_summary.tsv")
    report = validate_assembly_summary(bad)
    assert not report.is_valid
    assert any(f.field == "N50" for f in report.errors)


def test_duplicate_assembly_id_detected():
    bad = _load("example_assembly_broken_duplicate_assembly_ids", "assembly_summary.tsv")
    report = validate_assembly_summary(bad)
    assert not report.is_valid
    assert any("Duplicate Assembly_ID" in f.message for f in report.errors)


def test_quast_style_column_names_normalized_and_accepted():
    summary_tbl = _load("example_assembly", "assembly_summary.tsv")
    quast_style = summary_tbl.rename(columns={
        "Assembly_ID": "Assembly", "Total_length": "Total length", "GC_percent": "GC (%)",
    })
    report = validate_assembly_summary(quast_style)
    assert report.is_valid


def test_assembly_nx_plot_builds_one_line_per_assembly():
    contigs = _load("example_assembly", "contig_lengths.tsv")
    fig = assembly_nx_plot(contigs)
    assert isinstance(fig, Figure)
    ax = fig.axes[0]
    assert len(ax.get_lines()) == contigs["Assembly_ID"].nunique()


def test_assembly_summary_barplot_uses_n50_by_default():
    summary_tbl = _load("example_assembly", "assembly_summary.tsv")
    fig = assembly_summary_barplot(summary_tbl)
    assert isinstance(fig, Figure)
    assert fig.axes[0].get_ylabel() == "N50"


def test_assembly_summary_barplot_errors_on_unknown_stat_col():
    summary_tbl = _load("example_assembly", "assembly_summary.tsv")
    with pytest.raises(ValueError, match="not found"):
        assembly_summary_barplot(summary_tbl, stat_col="NotAColumn")
