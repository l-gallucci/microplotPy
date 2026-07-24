import pandas as pd

from microplotpy import filter_taxa


def make_agg():
    # 4 samples, 3 taxa: A abundant+prevalent, B abundant but only in 1 sample,
    # C low abundance but prevalent everywhere.
    return pd.DataFrame({
        "Sample_ID": ["S1", "S2", "S3", "S4"] * 3,
        "taxon": ["A"] * 4 + ["B"] * 4 + ["C"] * 4,
        "rel_abund": [20, 25, 22, 18, 40, 0, 0, 0, 1, 1, 1, 1],
    })


def test_no_filters_returns_every_taxon():
    out = filter_taxa(make_agg(), taxon_col="taxon")
    assert set(out) == {"A", "B", "C"}


def test_min_rel_abund_excludes_low_mean_abundance_taxa():
    out = filter_taxa(make_agg(), taxon_col="taxon", min_rel_abund=5)
    assert set(out) == {"A", "B"}


def test_min_prevalence_fraction_excludes_taxa_in_too_few_samples():
    out = filter_taxa(make_agg(), taxon_col="taxon", min_prevalence=0.5)
    assert set(out) == {"A", "C"}


def test_min_prevalence_above_one_is_absolute_sample_count():
    out = filter_taxa(make_agg(), taxon_col="taxon", min_prevalence=3)
    assert set(out) == {"A", "C"}


def test_both_filters_combine_with_and_logic():
    out = filter_taxa(make_agg(), taxon_col="taxon", min_rel_abund=5, min_prevalence=0.5)
    assert set(out) == {"A"}


def test_detection_threshold_changes_what_counts_as_present():
    out = filter_taxa(make_agg(), taxon_col="taxon", min_prevalence=0.5, detection=1)
    assert "C" not in out
