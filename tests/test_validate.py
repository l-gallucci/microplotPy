from pathlib import Path

import pandas as pd
import pytest

from microplotpy import MicrobiomeData, load, validate

EXAMPLES = Path(__file__).parent.parent / "examples"


def _load(name: str):
    d = EXAMPLES / name
    return load(d / "feature_table.tsv", d / "taxonomy.tsv", d / "metadata.tsv")


def test_valid_dataset_has_no_errors():
    data = _load("example_valid")
    report = validate(data)
    assert report.is_valid
    assert report.errors == []


def test_valid_dataset_gradient_and_group_checks_pass():
    data = _load("example_valid")
    report = validate(data, gradient_column="Depth_m", group_column="Group")
    assert report.is_valid


def test_id_mismatch_detected_both_directions():
    data = _load("example_broken_id_mismatch")
    report = validate(data)
    assert not report.is_valid
    messages = " ".join(f.message for f in report.errors)
    assert "S9" in messages  # metadata has S9, missing from feature_table
    assert "S8" in messages  # feature_table has S8, missing from metadata


def test_negative_counts_detected():
    data = _load("example_broken_negative_counts")
    report = validate(data)
    assert not report.is_valid
    assert any("Negative abundance" in f.message for f in report.errors)
    assert any(f.field == "values" for f in report.errors)


def test_missing_taxonomy_column_detected():
    data = _load("example_broken_missing_columns")
    report = validate(data)
    assert not report.is_valid
    assert any(f.field == "Phylum" for f in report.errors)


def test_duplicate_feature_ids_detected():
    data = _load("example_broken_duplicate_ids")
    report = validate(data)
    assert not report.is_valid
    assert any("Duplicate Feature_ID" in f.message for f in report.errors)


def test_missing_group_column_is_error():
    data = _load("example_valid")
    report = validate(data, group_column="NotAColumn")
    assert not report.is_valid


def test_single_level_group_column_is_error():
    data = _load("example_valid")
    # Group has 2 levels (A/B) in the valid set; use Sample_ID itself as a fake
    # "group" with 8 levels vs a constant column to trigger the <2-levels case.
    data.metadata["Constant"] = "only_one_value"
    report = validate(data, group_column="Constant")
    assert not report.is_valid


def test_required_ranks_can_be_relaxed_for_functional_annotation():
    ft = pd.DataFrame({"Feature_ID": ["g1", "g2"], "S1": [10, 20], "S2": [5, 15]})
    fa = pd.DataFrame({"Feature_ID": ["g1", "g2"], "COG_category": ["C", "J"], "KEGG_ko": ["K00001", "K00002"]})
    meta = pd.DataFrame({"Sample_ID": ["S1", "S2"]})
    data = MicrobiomeData(feature_table=ft, taxonomy=fa, metadata=meta)

    report_default = validate(data)
    assert not report_default.is_valid  # Phylum/Genus missing by default

    report_relaxed = validate(data, required_ranks=())
    assert report_relaxed.is_valid
