from pathlib import Path

import pandas as pd
import pytest
from matplotlib.figure import Figure

from microplotpy import validate_mag
from microplotpy.plots import mag_quality_distribution, mag_quality_plot

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str) -> pd.DataFrame:
    return pd.read_csv(EXAMPLES / name / "mag_quality.tsv", sep="\t")


def test_valid_mag_table_passes():
    mag = _load("example_mag_quality")
    report = validate_mag(mag)
    assert report.is_valid


def test_completeness_out_of_range_is_error():
    mag = _load("example_mag_broken_out_of_range")
    report = validate_mag(mag)
    assert not report.is_valid
    assert any("outside the valid" in f.message for f in report.errors)


def test_negative_contamination_is_error():
    mag = _load("example_mag_broken_negative")
    report = validate_mag(mag)
    assert not report.is_valid
    assert any("cannot be negative" in f.message for f in report.errors)


def test_duplicate_names_detected():
    mag = _load("example_mag_broken_duplicate_ids")
    report = validate_mag(mag)
    assert not report.is_valid
    assert any("Duplicate MAG Name" in f.message for f in report.errors)


def test_missing_required_column_detected():
    mag = _load("example_mag_broken_missing_columns")
    report = validate_mag(mag)
    assert not report.is_valid
    assert any(f.field == "Contamination" for f in report.errors)


def test_checkm1_style_columns_normalized_and_accepted():
    mag = _load("example_mag_quality")
    checkm1_style = mag.rename(columns={"Name": "Bin Id", "Genome_Size": "Genome size (bp)"})
    report = validate_mag(checkm1_style)
    assert report.is_valid


def test_mag_quality_plot_builds_without_error():
    mag = _load("example_mag_quality")
    fig = mag_quality_plot(mag, size_col="Genome_Size", color_col="Phylum")
    assert isinstance(fig, Figure)


def test_mag_quality_plot_with_no_size_or_color():
    mag = _load("example_mag_quality")
    fig = mag_quality_plot(mag, size_col=None, color_col=None, show_mimag_thresholds=False)
    assert isinstance(fig, Figure)


def test_mag_quality_plot_errors_on_unknown_color_col():
    mag = _load("example_mag_quality")
    with pytest.raises(ValueError, match="not found"):
        mag_quality_plot(mag, color_col="NotAColumn")


def test_mag_quality_distribution_builds_two_panels():
    mag = _load("example_mag_quality")
    fig = mag_quality_distribution(mag)
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    titles = [ax.get_title() for ax in fig.axes]
    assert "Completeness (%)" in titles
    assert "Contamination (%)" in titles
