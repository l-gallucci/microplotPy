import pandas as pd

from microplotpy import tax_fix


def test_known_values_pass_through_unchanged():
    tax = pd.DataFrame({
        "Feature_ID": ["F1"], "Phylum": ["Firmicutes"],
        "Family": ["Lactobacillaceae"], "Genus": ["Lactobacillus"],
    })
    fixed = tax_fix(tax)
    assert fixed["Genus"].iloc[0] == "Lactobacillus"


def test_single_unknown_falls_back_to_last_known_ancestor():
    tax = pd.DataFrame({
        "Feature_ID": ["F1"], "Phylum": ["Firmicutes"],
        "Family": ["Lactobacillaceae"], "Genus": [None],
    })
    fixed = tax_fix(tax)
    assert fixed["Genus"].iloc[0] == "Unclassified_Lactobacillaceae"


def test_cascading_unknowns_anchor_to_same_true_ancestor():
    tax = pd.DataFrame({
        "Feature_ID": ["F1"], "Phylum": ["Firmicutes"], "Family": ["Lactobacillaceae"],
        "Genus": ["unclassified"], "Species": ["uncultured"],
    })
    fixed = tax_fix(tax)
    assert fixed["Genus"].iloc[0] == "Unclassified_Lactobacillaceae"
    assert fixed["Species"].iloc[0] == "Unclassified_Lactobacillaceae"


def test_unknown_with_no_known_ancestor_becomes_plain_unclassified():
    tax = pd.DataFrame({"Feature_ID": ["F1"], "Domain": [None], "Phylum": [None]})
    fixed = tax_fix(tax)
    assert fixed["Domain"].iloc[0] == "Unclassified"
    assert fixed["Phylum"].iloc[0] == "Unclassified"


def test_unknown_string_matching_case_insensitive_and_trimmed():
    tax = pd.DataFrame({
        "Feature_ID": ["F1", "F2", "F3"], "Phylum": ["Bacteroidota"] * 3,
        "Genus": ["UNCULTURED", "  ", "Unknown"],
    })
    fixed = tax_fix(tax)
    assert (fixed["Genus"] == "Unclassified_Bacteroidota").all()


def test_errors_informatively_when_no_recognized_rank_columns():
    tax = pd.DataFrame({"Feature_ID": ["F1"], "Notes": ["x"]})
    try:
        tax_fix(tax)
        assert False, "expected ValueError"
    except ValueError as e:
        assert "No recognized taxonomy rank columns" in str(e)
