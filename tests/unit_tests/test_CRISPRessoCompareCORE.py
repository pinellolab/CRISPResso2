"""Unit tests for CRISPRessoCompareCORE."""

from CRISPResso2 import CRISPRessoCompareCORE

from copy import deepcopy
import pytest


@pytest.fixture(scope="function")
def run_info():
    """Fixture providing run_info dictionary for tests."""
    return {
        "results": {
            "refs": {
                "Reference": {
                    "sequence": (
                        "CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG"
                    ),
                    "sgRNA_orig_sequences": ["GGCCCTTAAAA"],
                    "sgRNA_cut_points": [50],
                    "allele_frequency_files": ["Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt"],
                },
            },
        },
    }


@pytest.fixture(scope="function")
def run_info_1(run_info):
    """Fixture providing a copy of run_info for test 1."""
    return deepcopy(run_info)


@pytest.fixture(scope="function")
def run_info_2(run_info):
    """Fixture providing a copy of run_info for test 2."""
    return deepcopy(run_info)


def test_get_matching_allele_files(run_info):
    """Test get_matching_allele_files with identical run_info."""
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info, run_info)
    assert matching_allele_files == [("Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt", "Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt")]


def test_get_matching_allele_files_different_cut_points(run_info_1, run_info_2):
    """Test get_matching_allele_files with different cut points."""
    run_info_2["results"]["refs"]["Reference"]["sgRNA_cut_points"] = [50, 51]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_different_guides(run_info_1, run_info_2):
    """Test get_matching_allele_files with different guides."""
    run_info_2["results"]["refs"]["Reference"]["sgRNA_orig_sequences"] = ["GGCCCTTAAAC"]
    run_info_2["results"]["refs"]["Reference"]["allele_frequency_files"] = ["Alleles_frequency_table_around_sgRNA_GGCCCTTAAAC.txt"]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_multiple_alleles(run_info_1, run_info_2):
    """Test get_matching_allele_files with multiple alleles."""
    run_info_1["results"]["refs"]["Other_Amplicon"] = deepcopy(run_info_1["results"]["refs"]["Reference"])
    run_info_1["results"]["refs"]["Other_Amplicon"]["sequence"] = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    run_info_1["results"]["refs"]["Other_Amplicon"]["allele_frequency_files"] = [
        "Other_Amplicon.Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt"
    ]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [("Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt", "Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt")]


def test_get_matching_allele_files_different_amplicon_names_same_sequence(run_info_1, run_info_2):
    """Test get_matching_allele_files with different amplicon names but same sequence."""
    run_info_2["results"]["refs"]["Other_Amplicon"] = deepcopy(run_info_1["results"]["refs"]["Reference"])
    del run_info_2["results"]["refs"]["Reference"]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [("Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt", "Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt")]


def test_get_matching_allele_files_some_different_guides(run_info_1, run_info_2):
    """Test get_matching_allele_files with some different guides."""
    run_info_1["results"]["refs"]["Reference"]["sgRNA_orig_sequences"] += ["AAAAAAAAAAAAAAAAAAA"]
    run_info_1["results"]["refs"]["Reference"]["allele_frequency_files"] += ["Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt"]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [("Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt", "Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt")]


def test_get_matching_allele_files_multiple_guides(run_info_1, run_info_2):
    """Test get_matching_allele_files with multiple guides."""
    run_info_1["results"]["refs"]["Reference"]["sgRNA_orig_sequences"] += ["AAAAAAAAAAAAAAAAAAA"]
    run_info_1["results"]["refs"]["Reference"]["allele_frequency_files"] += ["Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt"]
    run_info_2["results"]["refs"]["Reference"]["sgRNA_orig_sequences"] += ["AAAAAAAAAAAAAAAAAAA"]
    run_info_2["results"]["refs"]["Reference"]["allele_frequency_files"] += ["Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt"]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [
        ("Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt", "Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt"),
        ("Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt", "Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt"),
    ]


# =============================================================================
# Additional edge case tests
# =============================================================================


def test_get_matching_allele_files_empty_refs(run_info_1, run_info_2):
    """Test get_matching_allele_files with empty refs in one run."""
    run_info_1["results"]["refs"] = {}
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_both_empty_refs():
    """Test get_matching_allele_files with empty refs in both runs."""
    run_info_1 = {"results": {"refs": {}}}
    run_info_2 = {"results": {"refs": {}}}
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_different_sequences(run_info_1, run_info_2):
    """Test get_matching_allele_files with completely different sequences."""
    run_info_2["results"]["refs"]["Reference"]["sequence"] = "TTTTTTTTTTTTTTTTTTTTTTTTTT"
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_empty_allele_files(run_info_1, run_info_2):
    """Test get_matching_allele_files with empty allele_frequency_files."""
    run_info_1["results"]["refs"]["Reference"]["allele_frequency_files"] = []
    run_info_2["results"]["refs"]["Reference"]["allele_frequency_files"] = []
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_empty_sgrna_sequences(run_info_1, run_info_2):
    """Test get_matching_allele_files with empty sgRNA sequences."""
    run_info_1["results"]["refs"]["Reference"]["sgRNA_orig_sequences"] = []
    run_info_2["results"]["refs"]["Reference"]["sgRNA_orig_sequences"] = []
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_multiple_amplicons_match(run_info_1, run_info_2):
    """Test get_matching_allele_files with multiple matching amplicons."""
    # Add a second amplicon to both
    run_info_1["results"]["refs"]["Amplicon2"] = deepcopy(run_info_1["results"]["refs"]["Reference"])
    run_info_1["results"]["refs"]["Amplicon2"]["sequence"] = "GGGGGGGGGGGGGGGGGGGGGGGGGG"
    run_info_1["results"]["refs"]["Amplicon2"]["sgRNA_orig_sequences"] = ["GGGGGGGGGG"]
    run_info_1["results"]["refs"]["Amplicon2"]["allele_frequency_files"] = ["Alleles_frequency_table_around_sgRNA_GGGGGGGGGG.txt"]

    run_info_2["results"]["refs"]["Amplicon2"] = deepcopy(run_info_1["results"]["refs"]["Amplicon2"])

    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert len(matching_allele_files) == 2


def test_get_matching_allele_files_order_independence(run_info_1, run_info_2):
    """Test that order of run_info arguments doesn't matter for matching."""
    matching_1_2 = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    matching_2_1 = CRISPRessoCompareCORE.get_matching_allele_files(run_info_2, run_info_1)
    # Both should find matching files (though order within tuples may differ)
    assert len(matching_1_2) == len(matching_2_1)
