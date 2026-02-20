"""Unit tests for CRISPResso2CORE."""
import os
import pytest
import pandas as pd
from pytest_check import check

from CRISPResso2 import CRISPRessoCORE, CRISPRessoShared, CRISPRessoCOREResources

def calc_score(seq, ref):
    score = 0
    for seq_i, ref_i in zip(seq, ref):
        if seq_i == ref_i:
            score += 1
    return (score / float(len(seq))) * 100.0


# =============================================================================
# Tests for get_consensus_alignment_from_pairs
# =============================================================================


def test_get_consensus_alignment_from_pairs():
    """Tests for generating consensus alignments from paired reads."""
    try:
        CRISPRessoCORE.get_consensus_alignment_from_pairs
    except AttributeError:
        pytest.xfail('get_consensus_alignment_from_pairs is not implemented yet!')

    print("testing Easy")

    # basic test
    qual1 = "AAAA"
    aln1_seq = "--CGAT----"
    aln1_ref = "ATCGATCGAT"
    aln2_seq = "-----TCGAT"
    aln2_ref = "ATCGATCGAT"
    qual2 = "AAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGATCGAT"
        assert ref_seq == "ATCGATCGAT"
        assert score == 80
        assert caching_ok

    # test quality difference
    qual1 = "AAAB"
    aln1_seq = "--CGAT----"
    aln1_ref = "ATCGATCGAT"
    aln2_seq = "-----GCGAT"
    aln2_ref = "ATCGATCGAT"
    qual2 = "AAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGATCGAT"
        assert ref_seq == "ATCGATCGAT"
        assert score == 80
        assert not caching_ok

    # test quality difference
    qual1 = "AAAA"
    aln1_seq = "--CGAT----"
    aln1_ref = "ATCGATCGAT"
    aln2_seq = "-----GCGAT"
    aln2_ref = "ATCGATCGAT"
    qual2 = "BAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGAGCGAT"
        assert ref_seq == "ATCGATCGAT"
        assert score == 70
        assert not caching_ok

    # gaps between r1 and r2
    qual1 = "AAAAAAAAAA"
    aln1_seq = "--CGA-----"
    aln1_ref = "ATCGATCGAT"
    aln2_seq = "-------GA-"
    aln2_ref = "ATCGATCGAT"
    qual2 = "AAAAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGANNGAN"
        assert ref_seq == "ATCGATCGAT"
        assert score == 50
        assert caching_ok

    print('Finished easy tests... now for the hard stuff')

    # insertion in r1
    qual1 = "AAAA"
    aln1_seq = "--CCGA-----".replace(" ", "")  # added replace for vertical alignment
    aln1_ref = "ATC-GATCGAT".replace(" ", "")
    aln2_seq = "--- ----GA-".replace(" ", "")
    aln2_ref = "ATC GATCGAT".replace(" ", "")
    qual2 = "AA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCCGANNGAN"
        assert ref_seq == "ATC-GATCGAT"
        assert score == 45.455
        assert caching_ok

    # deletion in r1
    qual1 = "AA"
    aln1_seq = "--C-A-----".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "-------GA-".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "AA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNC-ANNGAN"
        assert ref_seq == "ATCGATCGAT"
        assert score == 40
        assert caching_ok

    # insertion in r2
    qual1 = "AAAA"
    aln1_seq = "--CCGA-----".replace(" ", "")  # added replace for vertical alignment
    aln1_ref = "ATCCGATC AT".replace(" ", "")
    aln2_seq = "--------GA-".replace(" ", "")
    aln2_ref = "ATCCGATC-AT".replace(" ", "")
    qual2 = "AA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCCGANNGAN"
        assert ref_seq == "ATCCGATC-AT"
        assert score == 45.455
        assert caching_ok

    # deletion in r2
    qual1 = "AAA"
    aln1_seq = "--CGA-----".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "-----T-GA-".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "AAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGAT-GAN"
        assert ref_seq == "ATCGATCGAT"
        assert score == 60
        assert caching_ok

    # insertion in r1 and r2
    qual1 = "AAAAAA"
    aln1_seq = "--CGATCC---".replace(" ", "")
    aln1_ref = "ATCGAT-CGAT".replace(" ", "")
    aln2_seq = "----ATACGA-".replace(" ", "")
    aln2_ref = "ATCGAT-CGAT".replace(" ", "")
    qual2 = "AAAAAA ".replace(" ", "")

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGATCCGAN"
        assert ref_seq == "ATCGAT-CGAT"
        assert score == 63.636
        assert not caching_ok

    # insertion in r1 and r2, different positions
    qual1 = "AAAAAA"
    aln1_seq = "--CGATCC---".replace(" ", "")
    aln1_ref = "ATCGAT-CGAT".replace(" ", "")
    aln2_seq = "----ATATGA-".replace(" ", "")
    aln2_ref = "ATCGATC-GAT".replace(" ", "")
    qual2 = "AAAAAA ".replace(" ", "")

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "NNCGATCCTGAN"
        assert ref_seq == "ATCGAT-C-GAT"
        assert score == 58.333
        assert not caching_ok

    # insertion at beginning of r1
    qual1 = "AAAAA"
    aln1_seq = "TA-CGA----- ".replace(" ", "")
    aln1_ref = "-ATCGATCGAT ".replace(" ", "")
    aln2_seq = " --------AT ".replace(" ", "")
    aln2_ref = " ATCGATCGAT".replace(" ", "")
    qual2 = "AA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "A-CGANNNAT"
        assert ref_seq == "ATCGATCGAT"
        assert score == 60.0
        assert caching_ok

    # insertion at end of r2 and beginning of r1
    qual1 = "AAAAA"
    aln1_seq = "TA-CGA-----   ".replace(" ", "")
    aln1_ref = "-ATCGATCGAT   ".replace(" ", "")
    aln2_seq = " -----TCGATCCA".replace(" ", "")
    aln2_ref = " ATCGATCGAT---".replace(" ", "")
    qual2 = "AAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == "A-CGATCGAT"
        assert ref_seq == "ATCGATCGAT"
        expected_score = 90.0
        assert score == expected_score
        assert caching_ok

    qual1 = '>1>1A@DFAADAGGGGGGGGGGHHHHHHHHHHHHHHHGGHHHHHHHHGGGHHHHHHHHHGHHHHHHHHHHHHGHGGGGGGGGGGHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHGHHGGGGGHHHHG'
    aln1_seq = 'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGCGGC-----CTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT'
    aln1_ref = 'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT'
    aln2_seq = (
        'AACCACAGCC-----GAGCCTCTTGAAGCCATTCTTACAGATGATGAAC-CAGG--CCGCGTTGGGAGCTCCAGAAGGGGATCATGACCT----CCTCACCTGTGGGCAGTGCCAGATGAACTTCCCATTGGGGGACATT'
    )
    aln2_ref = (
        'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGC--CGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT-----'
    )
    qual2 = 'BCCDCCDFDDDDGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHGHGGGGGGGHGHGGHHHHHHHHHGGGGGHHHHHHHHHHHHHHHHHHHGGHGHHHHHHHHHHHHHHHHHHHHHHHGGGGGHH'

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        expected_aln_seq = (
            'AACCACCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAACTCATGATCCCCTTCTGGAGCTCCAAAAGGGGATCATGACCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGTTCCCATTTGGTT'
        )
        expected_ref_seq = (
            'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGC--CGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT'
        )
        expected_score = 86.667
        assert aln_seq == expected_aln_seq
        assert ref_seq == expected_ref_seq
        assert score == expected_score
        assert not caching_ok

    # alternating qualities
    qual1 = "BABABABABA"
    aln1_seq = "ACCAACCAAT".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "TTGGTTGGTT".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "ABABABABAB"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1,
        aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2
    )
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)
    check.is_false(caching_ok)

    # large insertion in r1
    qual1 = "AAAAAA"
    aln1_seq = "ACGTGA---------".replace(" ", "")
    aln1_ref = "A-----TCGATCGAT".replace(" ", "")
    aln2_seq = "------CGAT".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "AAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1,
        aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2
    )
    check.equal(aln_seq, "ACGTGANNNNNCGAT")
    check.equal(ref_seq, "A-----TCGATCGAT")
    check.equal(score, 33.333)
    check.is_true(caching_ok)

    # large insertion in r2
    qual1 = "AAAAA"
    aln1_seq = "ATCGA-----".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "-----TTAGCT---".replace(" ", "")
    aln2_ref = "ATCGAT---C-GAT".replace(" ", "")
    qual2 = "AAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1,
        aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2
    )
    check.equal(aln_seq, "ATCGATTAGCTNNN")
    check.equal(ref_seq, "ATCGAT---C-GAT")
    check.equal(score, 50)
    check.is_true(caching_ok)

    # Conflicts with reference
    qual1 = "AAAAAAAAAA"
    aln1_seq = "TAGCTAGCTA".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "TAGCTAGCTA".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "AAAAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1,
        aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2
    )
    check.equal(aln_seq, "TAGCTAGCTA")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 0)
    check.is_true(caching_ok)

    # Conflicts between reads
    qual1 = "AAAAAAAAAA"
    aln1_seq = "TAGCTAGCTA".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "ATCGATCGAT".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "AAAAAAAAAA"

    aln_seq, _aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1,
        aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2
    )
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)
    check.is_false(caching_ok)

    # Alternating reads
    qual1 = "AAAAAAAAAA"
    aln1_seq = "AT--AT--AT".replace(" ", "")
    aln1_ref = "ATCGATCGAT".replace(" ", "")
    aln2_seq = "--CG--CG--".replace(" ", "")
    aln2_ref = "ATCGATCGAT".replace(" ", "")
    qual2 = "AAAAAAAAAA"

    aln_seq, _aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1,
        aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2
    )
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)
    check.is_true(caching_ok)


# =============================================================================
# Tests for split_quant_window_coordinates
# =============================================================================


def test_split_quant_window_coordinates_single():
    """Test split_quant_window_coordinates with single range."""
    assert [(5, 10)] == CRISPRessoCORE.split_quant_window_coordinates('5-10')


def test_split_quant_window_coordinates_multiple():
    """Test split_quant_window_coordinates with multiple ranges."""
    assert CRISPRessoCORE.split_quant_window_coordinates('2-5_10-12') == [(2, 5), (10, 12)]


def test_split_quant_window_coordinates_error():
    """Test split_quant_window_coordinates raises exception with invalid input."""
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('a-5')


def test_split_quant_window_coordinates_empty():
    """Test split_quant_window_coordinates raises exception with empty range."""
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('_')


def test_split_quant_window_coordinates_partially_empty():
    """Test split_quant_window_coordinates raises exception with partially empty range."""
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('1-3_')


def test_split_quant_window_coordinates_blank():
    """Test split_quant_window_coordinates raises exception with blank input."""
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('')


def test_split_quant_window_coordinates_single_position():
    """Test split_quant_window_coordinates with single position range."""
    result = CRISPRessoCORE.split_quant_window_coordinates("5-5")
    assert result == [(5, 5)]


def test_split_quant_window_coordinates_large_range():
    """Test split_quant_window_coordinates with large range."""
    result = CRISPRessoCORE.split_quant_window_coordinates("0-1000")
    assert result == [(0, 1000)]


def test_split_quant_window_coordinates_many_ranges():
    """Test split_quant_window_coordinates with many ranges (underscore-separated)."""
    result = CRISPRessoCORE.split_quant_window_coordinates("0-10_20-30_40-50")
    assert len(result) == 3
    assert result[0] == (0, 10)
    assert result[1] == (20, 30)
    assert result[2] == (40, 50)


def test_split_quant_window_coordinates_three_ranges():
    """Test split_quant_window_coordinates with three ranges."""
    result = CRISPRessoCORE.split_quant_window_coordinates("1-10_20-30_40-50")
    assert len(result) == 3
    assert result == [(1, 10), (20, 30), (40, 50)]


def test_split_quant_window_coordinates_large_numbers():
    """Test split_quant_window_coordinates with large numbers."""
    result = CRISPRessoCORE.split_quant_window_coordinates("1000-2000")
    assert result == [(1000, 2000)]


def test_split_quant_window_coordinates_zero_start():
    """Test split_quant_window_coordinates starting at 0."""
    result = CRISPRessoCORE.split_quant_window_coordinates("0-5")
    assert result == [(0, 5)]


# =============================================================================
# Tests for get_include_idxs_from_quant_window_coordinates
# =============================================================================


def test_get_include_idxs_from_quant_window_coordinates():
    """Test get_include_idxs_from_quant_window_coordinates function."""
    quant_window_coordinates = '1-10_12-20'
    expected_idxs = [*list(range(1, 11)), *list(range(12, 21))]
    assert CRISPRessoCORE.get_include_idxs_from_quant_window_coordinates(quant_window_coordinates) == expected_idxs


def test_get_include_idxs_single_position():
    """Test get_include_idxs with single position range."""
    result = CRISPRessoCORE.get_include_idxs_from_quant_window_coordinates("5-5")
    assert result == [5]


def test_get_include_idxs_multiple_ranges():
    """Test get_include_idxs with multiple ranges."""
    result = CRISPRessoCORE.get_include_idxs_from_quant_window_coordinates("0-2_10-12")
    assert result == [0, 1, 2, 10, 11, 12]


def test_get_include_idxs_large_range():
    """Test get_include_idxs with large range."""
    result = CRISPRessoCORE.get_include_idxs_from_quant_window_coordinates("0-100")
    assert len(result) == 101
    assert result[0] == 0
    assert result[-1] == 100


# =============================================================================
# Tests for get_cloned_include_idxs_from_quant_window_coordinates
# =============================================================================


def test_get_cloned_include_idxs_from_quant_window_coordinates():
    """Test get_cloned_include_idxs_from_quant_window_coordinates function."""
    quant_window_coordinates = '1-10_12-20'
    ref = 'TTACCGAGTGCACAAGTGCACGT'
    aln = 'TTACCGAGTGCACAAGTGCACGT'
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    expected_idxs = [*list(range(1, 11)), *list(range(12, 21))]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_beginning():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion at beginning."""
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp insertion at the beginning (left)
    # Ind:                1111111111222
    #           01234567890123456789012
    # QWC:       |        | |       |
    ref = '-----TTACCGAGTGCACAAGTGCACGT'
    aln = 'AAGGTTTACCGAGTGCACAAGTGCACGT'
    # QWC:       |        | |       |
    # Ind:           111111111122222222
    #      0123456789012345678901234567
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    expected_idxs = [*list(range(6, 16)), *list(range(17, 26))]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_beginning():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion at beginning."""
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp deletion at the beginning (left)
    # Ind:           1111111111222
    #      01234567890123456789012
    # QWC:  |        | |       |
    ref = 'TTACCGAGTGCACAAGTGCACGT'
    aln = '------AGTGCACAAGTGCACGT'
    # QWC:       |   | |       |
    # Ind:                 1111111
    #            01234567890123456
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    expected_idxs = [*list(range(0, 5)), *list(range(6, 15))]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion."""
    quant_window_coordinates = '10-20_35-40'
    ref = 'A' * 23 + 'T' * 7 + 'G' * 30
    aln = 'A' * 23 + '-' * 7 + 'G' * 30
    # represents a 7bp deletion in the middle
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    deletion_size = 7
    expected_idxs = [*list(range(10, 21)), *list(range(35 - deletion_size, 41 - deletion_size))]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_modified():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with modified deletion."""
    quant_window_coordinates = '10-25_35-40'
    # represents a 7bp deletion in the middle, where part of the QW is deleted
    # [0, 1, 3, 4, ... , 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 24, ... , 33]
    ref = 'A' * 23 + 'T' * 7 + 'G' * 30
    aln = 'A' * 23 + '-' * 7 + 'G' * 30
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    deletion_size = 7
    expected_idxs = [*list(range(10, 23)), *list(range(35 - deletion_size, 41 - deletion_size))]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_end_modified():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion at end."""
    # 5 bp deletion at end of 20 bp sequence
    quant_window_coordinates = '1-5_10-20'
    # Ind:           11111111112
    #      012345678901234567890
    # QWC:  |   |    |         |
    ref = 'AAAAAAAAAAAAAAAATTTTT'
    aln = 'AAAAAAAAAAAAAAAA-----'
    # QWC:  |   |    |    |
    # Ind:           111111
    #      0123456789012345
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    expected_idxs = [*list(range(1, 6)), *list(range(10, 16))]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion and deletion."""
    # 5 bp deletion and 5 bp insertion
    quant_window_coordinates = '1-5_10-18'
    # Ind:           1111     11111
    #      01234567890123     45678
    # QWC:  |   |    |            |
    ref = 'AAAAATTTTTGGGG-----AAAAA'
    aln = 'AAAAA-----GGGGCCCCCAAAAA'
    # QWC:  |                     |
    # Ind:                111111111
    #      01234     56789012345678
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 19))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion_modified():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with modified insertion and deletion."""
    quant_window_coordinates = '1-5_10-20'
    # Ind:                 11111111112
    #      012 34567     8901234567890
    # QWC:  |    |         |         |
    ref = 'AAA-CCCCC-----AAATTTTCCCCCC'
    aln = 'AAAT-CCCCGGGGGAAA----CCCCCC'
    # QWC:  |    |         |         |
    # Ind:           1111111    111122
    #      0123 456789012345    678901
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    expected_idxs = [*[1, 2, 3, 4, 5], *[15, 16, 17, 18, 19, 20, 21]]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == expected_idxs


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion."""
    quant_window_coordinates = '2-7'

    # Ind: 0123  456789
    # QWC:   |      |
    ref = 'AAAA--TTTTTT'
    aln = 'AAAAGGTTTTTT'
    # QWC:   |      |
    # Ind:           11
    #      012345678901
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(2, 10))


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_start():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion at start."""
    quant_window_coordinates = '1-3'

    # Ind: 0    123456
    # QWC:      | |
    ref = 'T----ACTGT'
    aln = 'TTCCCACTGT'
    # QWC:      | |
    # Ind: 0123456789
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [5, 6, 7]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_end():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion at end."""
    quant_window_coordinates = '1-4'

    # Ind: 0123    45678
    # QWC:  |      |
    ref = 'GGGT----ACTGT'
    aln = 'GGGTTCCCACTGT'
    # QWC:  |      |
    # Ind:           111
    #      0123456789012
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(1, 9))


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_before_qw():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion before quantification window."""
    quant_window_coordinates = '6-9'

    # Ind:           11
    #      012345678901
    # QWC:       |  |
    ref = 'GGGTACTGTCCA'
    aln = 'GGGTA-TGTCCA'
    # QWC:       |  |
    # Ind:            1
    #      01234 567890
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(5, 9))


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_before_qw():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion before quantification window."""
    quant_window_coordinates = '6-9'
    # Ind:            1
    #      01234 567890
    # QWC:        |  |
    ref = 'GGGTA-TGTCCA'
    aln = 'GGGTACTGTCCA'
    # QWC:        |  |
    # Ind:           11
    #      012345678901
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(7, 11))


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_deletion_outside_qw():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion and deletion outside quantification window."""
    quant_window_coordinates = '4-6_11-14'

    # Ind:               11111
    #      0123    45678901234
    # QWC:         | |    |  |
    ref = 'GGGT----ACTTTTGTCCA'
    aln = 'GGGTTCCCACT---GTCCA'
    # QWC:         | |    |  |
    # Ind:           1   11111
    #      01234567890   12345
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [8, 9, 10, 12, 13, 14, 15]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_across_qw():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with insertion across quantification window."""
    # 6 bp insertion in middle of 4 bp sequence
    quant_window_coordinates = '1-3'
    # Ind: 01      23
    # QWC:  |       |
    ref = 'AA------TT'
    aln = 'AAGGGGGGTT'
    # QWC:  |       |
    # Ind: 0123456789
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(1, 10))


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_overlap_start():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion overlapping start of quantification window."""
    quant_window_coordinates = '2-5'
    # Ind: 012345
    # QWC:   |  |
    ref = 'AATTTT'
    aln = 'A--TTT'
    # QWC:    | |
    # Ind: 0  123
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2, 3]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_overlap_end():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion overlapping end."""
    quant_window_coordinates = '1-3'
    # Ind: 012345
    # QWC:  | |
    ref = 'AATTTT'
    aln = 'AAT---'
    # QWC:  ||
    # Ind: 012
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_overlap_end_single_bp():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with single bp deletion overlapping end."""
    quant_window_coordinates = '2-3'
    # Ind: 012345
    # QWC:   ||
    ref = 'AATTTT'
    aln = 'AAT---'
    # QWC:   |
    # Ind: 012
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [2]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_entire_qw():
    """Test get_cloned_include_idxs_from_quant_window_coordinates with deletion of entire quantification window."""
    # 5 bp deletion of entire qw
    quant_window_coordinates = '1-4_7-10'
    ref = 'AAAAAATTTT'
    aln = 'AAAAAA----'
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2, 3, 4]


def test_get_cloned_include_idxs_from_quant_window_coordinates_include_zero():
    """Test get_cloned_include_idxs_from_quant_window_coordinates including zero index."""
    quant_window_coordinates = '0-4'
    ref = 'AAAAA'
    aln = 'AAAAA'
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [0, 1, 2, 3, 4]


# =============================================================================
# Tests for get_variant_cache_equal_boundaries (parallelization functions)
# =============================================================================


def test_regular_input():
    """Test get_variant_cache_equal_boundaries with typical input."""
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(100, 4) == [0, 25, 50, 75, 100]


def test_remainder_input():
    """Test get_variant_cache_equal_boundaries with remainder input."""
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(101, 4) == [0, 25, 50, 75, 101]


def test_similar_num_reads_input():
    """Test get_variant_cache_equal_boundaries with similar number of reads and processes."""
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(11, 10) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]


def test_large_similar_num_reads_input():
    """Test get_variant_cache_equal_boundaries with large similar number of reads and processes."""
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(101, 100) == [*list(range(0, 100)), 101]


def test_more_processes_than_reads():
    """Test get_variant_cache_equal_boundaries raises exception when more processes than reads."""
    with pytest.raises(ValueError):
        CRISPRessoCORE.get_variant_cache_equal_boundaries(3, 5)


def test_single_process():
    """Test get_variant_cache_equal_boundaries with a single process."""
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(50, 1) == [0, 50]


def test_zero_sequences():
    """Test get_variant_cache_equal_boundaries raises exception with zero sequences."""
    with pytest.raises(ValueError):
        CRISPRessoCORE.get_variant_cache_equal_boundaries(0, 3)


def test_large_numbers():
    """Test get_variant_cache_equal_boundaries with large number of processes and sequences."""
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(10000, 10)
    expected_boundaries_count = 11  # n_processes + 1
    assert len(boundaries) == expected_boundaries_count


def test_sublist_generation():
    """Test sublist generation from variant cache boundaries."""
    n_processes = 4
    unique_reads = 100
    mock_variant_cache = [i for i in range(unique_reads)]
    assert len(mock_variant_cache) == unique_reads
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(unique_reads, n_processes)
    assert boundaries == [0, 25, 50, 75, 100]
    sublists = []
    for i in range(n_processes):
        left_sublist_index = boundaries[i]
        right_sublist_index = boundaries[i + 1]
        sublist = mock_variant_cache[left_sublist_index:right_sublist_index]
        sublists.append(sublist)
    assert [len(sublist) for sublist in sublists] == [25, 25, 25, 25]
    assert [s for sublist in sublists for s in sublist] == mock_variant_cache


def test_irregular_sublist_generation():
    """Test sublist generation with irregular number of reads."""
    n_processes = 4
    unique_reads = 113
    mock_variant_cache = [i for i in range(unique_reads)]
    assert len(mock_variant_cache) == unique_reads
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(unique_reads, n_processes)
    # assert boundaries == [0, 25, 50, 75, 100]
    sublists = []
    for i in range(n_processes):
        left_sublist_index = boundaries[i]
        right_sublist_index = boundaries[i + 1]
        sublist = mock_variant_cache[left_sublist_index:right_sublist_index]
        sublists.append(sublist)
    assert [len(sublist) for sublist in sublists] == [28, 28, 28, 29]
    assert [s for sublist in sublists for s in sublist] == mock_variant_cache


def test_get_variant_cache_equal_boundaries_exact_division():
    """Test boundaries with exact division."""
    result = CRISPRessoCORE.get_variant_cache_equal_boundaries(100, 5)
    assert result == [0, 20, 40, 60, 80, 100]


def test_get_variant_cache_equal_boundaries_two_processes():
    """Test boundaries with two processes."""
    result = CRISPRessoCORE.get_variant_cache_equal_boundaries(100, 2)
    assert result == [0, 50, 100]


def test_get_variant_cache_equal_boundaries_many_processes():
    """Test boundaries with many processes."""
    result = CRISPRessoCORE.get_variant_cache_equal_boundaries(1000, 8)
    assert len(result) == 9  # n_processes + 1


def test_get_variant_cache_equal_boundaries_uneven():
    """Test boundaries with uneven division."""
    result = CRISPRessoCORE.get_variant_cache_equal_boundaries(17, 4)
    # Should handle 17 reads across 4 processes
    assert len(result) == 5  # n+1 boundaries
    assert result[0] == 0
    assert result[-1] == 17


def test_get_variant_cache_equal_boundaries_prime():
    """Test boundaries with prime number of reads."""
    result = CRISPRessoCORE.get_variant_cache_equal_boundaries(37, 4)
    assert len(result) == 5
    assert result[-1] == 37


def test_get_variant_cache_equal_boundaries_small():
    """Test boundaries with small number of reads."""
    result = CRISPRessoCORE.get_variant_cache_equal_boundaries(5, 5)
    assert result == [0, 1, 2, 3, 4, 5]


# =============================================================================
# Tests for get_base_edit_target_sequence
# =============================================================================


def test_get_base_edit_target_sequence():
    """Test get_base_edit_target_sequence function."""
    df_alleles = pd.read_csv('tests/df_alleles.txt')
    ref_seq = (
        'CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG'
    )
    base_editor_target_ref_skip_allele_count = 0

    target_seq = CRISPRessoCORE.get_base_edit_target_sequence(
        ref_seq,
        df_alleles,
        base_editor_target_ref_skip_allele_count
    )

    expected_seq = (
        'AATACGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCCGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGGCATGGCCCCATTCGCACGGCTCTGGAGCGGC'
    )
    assert target_seq == expected_seq


def test_get_base_edit_target_sequence_alt():
    """Test get_base_edit_target_sequence function with alternative data."""
    df_alleles = pd.read_csv('tests/test_be_df.txt')
    ref_seq = 'AAAA'
    base_editor_target_ref_skip_allele_count = 0

    target_seq = CRISPRessoCORE.get_base_edit_target_sequence(
        ref_seq,
        df_alleles,
        base_editor_target_ref_skip_allele_count
    )

    assert target_seq == 'AAGA'


# =============================================================================
# Tests for get_bp_substitutions
# =============================================================================


def test_get_bp_subs_one_sub():
    """Test get_bp_substitutions with one substitution."""
    ref_seq = 'AAA'
    aln_ref_seq = 'AAA'
    aln_target_seq = 'CAA'
    ref_positions_to_include = [0, 1, 2]
    ref_changes_dict = CRISPRessoCORE.get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include)

    assert len(bp_substitutions_arr) == 1
    assert bp_substitutions_arr[0][0] == 0
    assert bp_substitutions_arr[0][1] == 'A'
    assert bp_substitutions_arr[0][2] == 'C'


def test_get_bp_subs_two_subs():
    """Test get_bp_substitutions with two substitutions."""
    ref_seq = 'AAA'
    aln_ref_seq = 'AAA'
    aln_target_seq = 'CCA'
    ref_positions_to_include = [0, 1, 2]
    ref_changes_dict = CRISPRessoCORE.get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include)

    expected_count = 2
    assert len(bp_substitutions_arr) == expected_count
    assert bp_substitutions_arr[0] == (0, 'A', 'C')
    assert bp_substitutions_arr[1] == (1, 'A', 'C')


def test_get_bp_subs_insertions():
    """Test get_bp_substitutions with insertions."""
    ref_seq = 'AAAAAA'
    aln_ref_seq = 'AAA-AAA-----'
    aln_target_seq = 'AAACAAACCCCC'
    ref_positions_to_include = [0, 1, 2, 3, 4, 5]
    ref_changes_dict = CRISPRessoCORE.get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include)

    expected_count = 2
    assert len(bp_substitutions_arr) == expected_count
    assert bp_substitutions_arr[0] == (2, 'A', 'AC')
    assert bp_substitutions_arr[1] == (5, 'A', 'ACCCCC')


def test_get_bp_substitutions_no_changes():
    """Test get_bp_substitutions with no substitutions."""
    ref_seq = 'ATCG'
    ref_changes_dict = CRISPRessoCORE.get_refpos_values('ATCG', 'ATCG')
    ref_positions = [0, 1, 2, 3]

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    assert result == []


def test_get_bp_substitutions_single_sub():
    """Test get_bp_substitutions with single A->G substitution at position 0."""
    ref_seq = 'ATCG'
    ref_changes_dict = CRISPRessoCORE.get_refpos_values('ATCG', 'GTCG')
    ref_positions = [0, 1, 2, 3]

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    assert result == [(0, 'A', 'G')]


def test_get_bp_substitutions_multiple_subs():
    """Test get_bp_substitutions with multiple substitutions."""
    ref_seq = 'ATCG'
    ref_changes_dict = CRISPRessoCORE.get_refpos_values('ATCG', 'GACT')
    ref_positions = [0, 1, 2, 3]

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    assert result == [(0, 'A', 'G'), (1, 'T', 'A'), (3, 'G', 'T')]


def test_get_bp_substitutions_partial_positions():
    """Test get_bp_substitutions only checks specified positions."""
    ref_seq = 'ATCG'
    ref_changes_dict = CRISPRessoCORE.get_refpos_values('ATCG', 'GTCA')
    ref_positions = [0, 3]  # Only check positions 0 and 3

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    assert result == [(0, 'A', 'G'), (3, 'G', 'A')]


def test_get_bp_substitutions_all_match():
    """Test get_bp_substitutions with identical read and reference."""
    ref_seq = 'ATCG'
    ref_changes_dict = CRISPRessoCORE.get_refpos_values('ATCG', 'ATCG')
    ref_positions = [0, 1, 2, 3]

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    assert result == []


def test_get_bp_substitutions_all_different():
    """Test get_bp_substitutions with all positions substituted."""
    ref_seq = 'ATCG'
    ref_changes_dict = CRISPRessoCORE.get_refpos_values('ATCG', 'TAGC')
    ref_positions = [0, 1, 2, 3]

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    assert result == [(0, 'A', 'T'), (1, 'T', 'A'), (2, 'C', 'G'), (3, 'G', 'C')]


def test_get_bp_substitutions_with_insertion():
    """Test get_bp_substitutions recognizes insertions."""
    ref_seq = 'ATCG'
    ref_changes_dict = {0: 'A', 1: 'TG', 2: 'C', 3: 'G'}  # Insertion at position 1
    ref_positions = [0, 1, 2, 3]

    result = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions)

    # Position 1 should show as substitution/insertion
    assert any(sub[0] == 1 for sub in result)


# =============================================================================
# Tests for get_refpos_values
# =============================================================================


def test_get_refpos_values_no_gaps():
    """Test get_refpos_values with no gaps."""
    ref_seq = "ATCG"
    read_seq = "ATCG"
    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)
    assert result[0] == "A"
    assert result[1] == "T"
    assert result[2] == "C"
    assert result[3] == "G"


def test_get_refpos_values_gap_in_ref():
    """Test get_refpos_values with gap in reference."""
    ref_seq = "--ATGC"
    read_seq = "GGATGC"
    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)
    assert result[0] == "GGA"  # Initial insertions go to position 0
    assert result[1] == "T"
    assert result[2] == "G"
    assert result[3] == "C"


def test_get_refpos_values_gap_in_read():
    """Test get_refpos_values with gap in read."""
    ref_seq = "ATGC"
    read_seq = "A-GC"
    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)
    assert result[0] == "A"
    assert result[1] == "-"
    assert result[2] == "G"
    assert result[3] == "C"


def test_get_refpos_values_example_from_docstring():
    """Test get_refpos_values with example from docstring."""
    ref_seq = "--A-TGC-"
    read_seq = "GGAGTCGA"
    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)
    assert result[0] == "GGAG"
    assert result[1] == "T"
    assert result[2] == "C"
    assert result[3] == "GA"


def test_get_refpos_values_insertion_middle():
    """Test get_refpos_values with insertion in middle."""
    ref_seq = "AT-CG"
    read_seq = "ATGCG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    assert result[0] == "A"
    assert result[1] == "TG"  # Insertion attached to previous position
    assert result[2] == "C"
    assert result[3] == "G"


def test_get_refpos_values_deletion():
    """Test get_refpos_values with deletion."""
    ref_seq = "ATCG"
    read_seq = "A--G"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    assert result[0] == "A"
    assert result[1] == "-"
    assert result[2] == "-"
    assert result[3] == "G"


def test_get_refpos_values_complex():
    """Test get_refpos_values with complex alignment."""
    ref_seq = "A-TC--G"
    read_seq = "AGTCAAG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    assert result[0] == "AG"  # A + inserted G
    assert result[1] == "T"
    assert result[2] == "CAA"  # C + inserted AA
    assert result[3] == "G"


def test_get_refpos_values_all_matches():
    """Test get_refpos_values with all matching positions."""
    ref_seq = "ATCGATCG"
    read_seq = "ATCGATCG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    for i, base in enumerate("ATCGATCG"):
        assert result[i] == base


def test_get_refpos_values_insertion_at_start():
    """Test get_refpos_values with insertion at start."""
    ref_seq = "--ATCG"
    read_seq = "GGATCG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    # First ref position should have leading insertions + aligned base
    assert "GGA" == result[0]


def test_get_refpos_values_insertion_at_end():
    """Test get_refpos_values with insertion at end."""
    ref_seq = "ATCG--"
    read_seq = "ATCGGG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    # Last position should have aligned base + trailing insertions
    assert "GGG" == result[3]


def test_get_refpos_values_deletions():
    """Test get_refpos_values with deletions."""
    ref_seq = "ATCGATCG"
    read_seq = "AT--ATCG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    assert result[2] == "-"
    assert result[3] == "-"


def test_get_refpos_values_mixed_indels():
    """Test get_refpos_values with mixed insertions and deletions."""
    ref_seq = "AT-CGATCG"
    read_seq = "ATG--ATCG"

    result = CRISPRessoCORE.get_refpos_values(ref_seq, read_seq)

    # Check that we got expected length
    assert len(result) > 0


# =============================================================================
# Tests for get_greater_qual_nuc
# =============================================================================


def test_get_greater_qual_nuc_same_nucleotides():
    """Test get_greater_qual_nuc when nucleotides are the same."""
    nuc, decision_made, qual = CRISPRessoCORE.get_greater_qual_nuc("A", "I", "A", "I", True)
    assert nuc == "A"
    assert decision_made is False
    assert qual == "I"


def test_get_greater_qual_nuc_different_quality():
    """Test get_greater_qual_nuc with different quality scores."""
    # Higher quality is 'J' (ASCII 74) vs 'I' (ASCII 73)
    nuc, decision_made, qual = CRISPRessoCORE.get_greater_qual_nuc("A", "J", "T", "I", True)
    assert nuc == "A"
    assert decision_made is True
    assert qual == "J"


def test_get_greater_qual_nuc_r2_better():
    """Test get_greater_qual_nuc when R2 has better quality."""
    nuc, decision_made, qual = CRISPRessoCORE.get_greater_qual_nuc("A", "I", "T", "J", True)
    assert nuc == "T"
    assert decision_made is True
    assert qual == "J"


def test_get_greater_qual_nuc_equal_quality_r1_best():
    """Test get_greater_qual_nuc with equal quality, R1 has better alignment."""
    nuc, decision_made, qual = CRISPRessoCORE.get_greater_qual_nuc("A", "I", "T", "I", True)
    assert nuc == "A"
    assert decision_made is True


def test_get_greater_qual_nuc_equal_quality_r2_best():
    """Test get_greater_qual_nuc with equal quality, R2 has better alignment."""
    nuc, decision_made, qual = CRISPRessoCORE.get_greater_qual_nuc("A", "I", "T", "I", False)
    assert nuc == "T"
    assert decision_made is True


def test_get_greater_qual_nuc_equal_nucs_different_qual():
    """Test when nucleotides are same but quality differs."""
    nuc, decision, qual = CRISPRessoCORE.get_greater_qual_nuc("A", "H", "A", "J", True)
    assert nuc == "A"
    assert decision is False  # No decision needed since nucs are same
    assert qual == "J"  # Higher quality


def test_get_greater_qual_nuc_different_nucs_equal_qual():
    """Test when nucleotides differ but quality is equal."""
    # When equal quality, use is_best_aln_r1 to decide
    nuc1, decision1, qual1 = CRISPRessoCORE.get_greater_qual_nuc("A", "I", "T", "I", True)
    nuc2, decision2, qual2 = CRISPRessoCORE.get_greater_qual_nuc("A", "I", "T", "I", False)

    assert nuc1 == "A"  # R1 is best
    assert nuc2 == "T"  # R2 is best
    assert decision1 is True
    assert decision2 is True


# =============================================================================
# Tests for normalize_name
# =============================================================================


def test_normalize_name_provided_name():
    """Test normalize_name with provided name."""
    result = CRISPRessoCORE.normalize_name("test_sample", None, None, None)
    assert result == "test_sample"


def test_normalize_name_from_r1():
    """Test normalize_name derives name from R1 fastq."""
    result = CRISPRessoCORE.normalize_name(None, "/path/to/sample.fastq.gz", None, None)
    assert result == "sample"


def test_normalize_name_from_r1_r2():
    """Test normalize_name derives name from R1 and R2 fastq."""
    result = CRISPRessoCORE.normalize_name(None, "/path/to/sample_R1.fastq", "/path/to/sample_R2.fastq", None)
    assert result == "sample_R1_sample_R2"


def test_normalize_name_from_bam():
    """Test normalize_name derives name from BAM file."""
    result = CRISPRessoCORE.normalize_name(None, None, None, "/path/to/sample.bam")
    assert result == "sample"


def test_normalize_name_cleans_special_chars():
    """Test normalize_name cleans special characters."""
    result = CRISPRessoCORE.normalize_name("test:sample/with*special", None, None, None)
    assert result == "test_sample_with_special"


def test_normalize_name_empty_inputs():
    """Test normalize_name with all empty inputs returns None."""
    result = CRISPRessoCORE.normalize_name(None, None, None, None)
    assert result is None


def test_normalize_name_complex_special_chars():
    """Test normalize_name handles multiple special characters."""
    result = CRISPRessoCORE.normalize_name("test:sample/with*[special]chars", None, None, None)
    assert result == "test_sample_with_special_chars"


def test_normalize_name_preserves_underscores():
    """Test normalize_name preserves underscores."""
    result = CRISPRessoCORE.normalize_name("test_sample_name", None, None, None)
    assert result == "test_sample_name"


def test_normalize_name_fastq_extensions():
    """Test normalize_name handles various fastq extensions."""
    result = CRISPRessoCORE.normalize_name(None, "/path/sample.fq.gz", None, None)
    assert result == "sample"


def test_normalize_name_r1_r2_suffix():
    """Test normalize_name with R1/R2 suffix files."""
    result = CRISPRessoCORE.normalize_name(
        None,
        "/path/sample_L001_R1_001.fastq.gz",
        "/path/sample_L001_R2_001.fastq.gz",
        None
    )
    assert result == "sample_L001_R1_001_sample_L001_R2_001"


def test_normalize_name_unicode_decomposable():
    """Test normalize_name decomposes accented characters to base letters."""
    # ë -> e, ä -> a via NFKD normalization in slugify
    result = CRISPRessoCORE.normalize_name("tëst_sämple", None, None, None)
    assert result == "test_sample"


def test_normalize_name_unicode_non_decomposable():
    """Test normalize_name drops non-decomposable unicode characters."""
    result = CRISPRessoCORE.normalize_name("sample_日本語", None, None, None)
    assert result == "sample_"


# =============================================================================
# Tests for to_numeric_ignore_columns
# =============================================================================


def test_to_numeric_ignore_columns_basic():
    """Test to_numeric_ignore_columns with basic dataframe."""
    df = pd.DataFrame({
        "name": ["a", "b", "c"],
        "value": ["1", "2", "3"],
        "count": ["10", "20", "30"]
    })
    result = CRISPRessoCORE.to_numeric_ignore_columns(df, {"name"})
    assert result["value"].dtype in [int, float, "int64", "float64"]
    assert result["count"].dtype in [int, float, "int64", "float64"]
    assert result["name"].dtype == object


def test_to_numeric_ignore_columns_multiple_ignore():
    """Test to_numeric_ignore_columns ignoring multiple columns."""
    df = pd.DataFrame({
        "name": ["a", "b"],
        "label": ["x", "y"],
        "value": ["1", "2"]
    })
    result = CRISPRessoCORE.to_numeric_ignore_columns(df, {"name", "label"})
    assert result["name"].dtype == object
    assert result["label"].dtype == object
    assert result["value"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_empty_ignore():
    """Test to_numeric_ignore_columns with empty ignore set."""
    df = pd.DataFrame({
        "a": ["1", "2"],
        "b": ["3", "4"]
    })
    result = CRISPRessoCORE.to_numeric_ignore_columns(df, set())
    assert result["a"].dtype in [int, float, "int64", "float64"]
    assert result["b"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_all_numeric():
    """Test to_numeric_ignore_columns with all numeric columns."""
    df = pd.DataFrame({
        "a": ["1", "2", "3"],
        "b": ["4.5", "5.5", "6.5"]
    })

    result = CRISPRessoCORE.to_numeric_ignore_columns(df, set())

    assert result["a"].dtype in [int, float, "int64", "float64"]
    assert result["b"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_preserve_strings():
    """Test to_numeric_ignore_columns preserves string columns in ignore list."""
    df = pd.DataFrame({
        "name": ["abc", "def", "ghi"],
        "value": ["1", "2", "3"]
    })

    result = CRISPRessoCORE.to_numeric_ignore_columns(df, {"name"})

    assert result["name"].dtype == object
    assert list(result["name"]) == ["abc", "def", "ghi"]


def test_to_numeric_ignore_columns_float_strings():
    """Test to_numeric_ignore_columns with float strings."""
    df = pd.DataFrame({
        "val": ["1.5", "2.5", "3.5"]
    })

    result = CRISPRessoCORE.to_numeric_ignore_columns(df, set())

    assert result["val"].dtype == float


def test_to_numeric_ignore_columns_mixed_types():
    """Test to_numeric_ignore_columns with mixed data types."""
    df = pd.DataFrame({
        "name": ["a", "b", "c"],
        "int_val": ["1", "2", "3"],
        "float_val": ["1.1", "2.2", "3.3"]
    })

    result = CRISPRessoCORE.to_numeric_ignore_columns(df, {"name"})

    assert result["name"].dtype == object
    assert result["int_val"].dtype in [int, float, "int64", "float64"]
    assert result["float_val"].dtype in [float, "float64"]


# =============================================================================
# Tests for check_library
# =============================================================================


def test_check_library_installed():
    """Test check_library returns module for installed library."""
    result = CRISPRessoCORE.check_library("os")
    assert result is not None
    import os as os_module
    assert result == os_module


def test_check_library_pandas():
    """Test check_library returns pandas module."""
    result = CRISPRessoCORE.check_library("pandas")
    assert result is not None


def test_check_library_numpy():
    """Test check_library with numpy."""
    result = CRISPRessoCORE.check_library("numpy")
    assert result is not None


def test_check_library_sys():
    """Test check_library with sys module."""
    result = CRISPRessoCORE.check_library("sys")
    assert result is not None
    import sys
    assert result == sys


# =============================================================================
# Tests for which
# =============================================================================


def test_which_existing_program():
    """Test which returns path for existing program."""
    result = CRISPRessoCORE.which("python")
    assert result is not None
    assert "python" in result


def test_which_nonexistent_program():
    """Test which returns None for nonexistent program."""
    result = CRISPRessoCORE.which("nonexistent_program_xyz_12345")
    assert result is None


def test_which_with_path():
    """Test which with full path."""
    import sys
    result = CRISPRessoCORE.which(sys.executable)
    assert result == sys.executable


def test_which_ls():
    """Test which finds ls command."""
    result = CRISPRessoCORE.which("ls")
    assert result is not None


def test_which_bash():
    """Test which finds bash."""
    result = CRISPRessoCORE.which("bash")
    assert result is not None


# =============================================================================
# Tests for get_scores_and_counts
# =============================================================================


def test_get_scores_and_counts_basic():
    """Test get_scores_and_counts with basic variant dictionary."""
    variant_dict = {
        'seq1': {'aln_scores': [95.0, 90.0], 'count': 100},
        'seq2': {'aln_scores': [80.0], 'count': 50},
    }

    homology_scores, counts, alleles_data = CRISPRessoCORE.get_scores_and_counts(variant_dict)

    assert len(homology_scores) == 2
    assert 95.0 in homology_scores
    assert 80.0 in homology_scores
    assert sum(counts) == 150


def test_get_scores_and_counts_empty():
    """Test get_scores_and_counts with empty dictionary."""
    variant_dict = {}

    homology_scores, counts, alleles_data = CRISPRessoCORE.get_scores_and_counts(variant_dict)

    assert len(homology_scores) == 0
    assert len(counts) == 0


def test_get_scores_and_counts_single_entry():
    """Test get_scores_and_counts with single variant."""
    variant_dict = {
        'ATCG': {'aln_scores': [100.0], 'count': 1000},
    }

    homology_scores, counts, alleles_data = CRISPRessoCORE.get_scores_and_counts(variant_dict)

    assert homology_scores == [100.0]
    assert counts == [1000]


def test_get_scores_and_counts_high_scores():
    """Test get_scores_and_counts with high alignment scores."""
    variant_dict = {
        'perfect': {'aln_scores': [100.0], 'count': 500},
        'good': {'aln_scores': [95.0, 90.0], 'count': 300},
    }

    scores, counts, data = CRISPRessoCORE.get_scores_and_counts(variant_dict)

    assert 100.0 in scores
    assert 95.0 in scores
    assert sum(counts) == 800


def test_get_scores_and_counts_low_scores():
    """Test get_scores_and_counts with low alignment scores."""
    variant_dict = {
        'poor': {'aln_scores': [50.0], 'count': 10},
    }

    scores, counts, data = CRISPRessoCORE.get_scores_and_counts(variant_dict)

    assert scores == [50.0]
    assert counts == [10]


def test_get_scores_and_counts_returns_alleles_data():
    """Test get_scores_and_counts returns proper alleles data structure."""
    variant_dict = {
        'seq1': {'aln_scores': [90.0], 'count': 100},
    }

    scores, counts, data = CRISPRessoCORE.get_scores_and_counts(variant_dict)

    assert len(data) == 1
    assert data[0]['sequence'] == 'seq1'
    assert data[0]['homology_score'] == 90.0
    assert data[0]['count'] == 100


# =============================================================================
# Tests for get_upset_plot_counts
# =============================================================================


def test_get_upset_plot_counts():
    """Test get_upset_plot_counts function."""
    df_alleles = pd.read_csv('tests/test_be_df.txt')
    bp_substitutions_arr = [(3, 'A', 'G')]

    wt_ref_name = 'TEST'

    counts_dict = CRISPRessoCORE.get_upset_plot_counts(
        df_alleles,
        bp_substitutions_arr,
        wt_ref_name
    )

    expected_dict_len = 19
    expected_total_alleles = 3
    expected_total_reads = 100
    assert len(counts_dict) == expected_dict_len
    assert counts_dict['total_alleles'] == expected_total_alleles
    assert counts_dict['total_alleles_reads'] == expected_total_reads
    assert sum(counts_dict['binary_allele_counts'].values()) == expected_total_reads


# =============================================================================
# Tests for write_base_edit_counts
# =============================================================================


def test_write_base_edit_counts():
    """Test write_base_edit_counts function."""
    output_directory = '.'
    clean_file_prefix = ""

    def _jp(filename):
        return os.path.join(output_directory, clean_file_prefix + filename)
    ref_name = 'TEST'
    bp_substitutions_arr = [(3, 'A', 'G')]
    counts_dict = CRISPRessoCORE.get_upset_plot_counts(
        pd.read_csv('tests/test_be_df.txt'),
        bp_substitutions_arr,
        ref_name,
    )

    expected_files = [
        '10i.TEST.arrays.txt',
        '10i.TEST.binary_allele_counts.txt',
        '10i.TEST.category_allele_counts.txt',
        '10i.TEST.counts.txt',
        '10i.TEST.precise_allele_counts.txt',
    ]

    CRISPRessoCORE.write_base_edit_counts(
        ref_name,
        counts_dict,
        bp_substitutions_arr,
        _jp
    )

    for filename in expected_files:
        if os.path.exists(_jp(filename)):
            os.remove(_jp(filename))
        else:
            raise AssertionError()

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "NNC-ANNGAN"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 40
        assert caching_ok

    #insertion in r2
    qual1                =   "AAAA"
    aln1_seq             = "--CCGA-----".replace(" ","") #added replace for vertical alignment
    aln1_ref             = "ATCCGATC AT".replace(" ","")
    aln2_seq             = "--------GA-".replace(" ","")
    aln2_ref             = "ATCCGATC-AT".replace(" ","")
    qual2                =         "AA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "NNCCGANNGAN"
        assert ref_seq ==      "ATCCGATC-AT"
        assert score == 45.455
        assert caching_ok

    # deletion in r2
    qual1                =   "AAA"
    aln1_seq             = "--CGA-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-----T-GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =      "AAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "NNCGAT-GAN"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 60
        assert caching_ok

    # insertion in r1 and r2
    qual1                =   "AAAAAA"
    aln1_seq             = "--CGATCC---".replace(" ","")
    aln1_ref             = "ATCGAT-CGAT".replace(" ","")
    aln2_seq             = "----ATACGA-".replace(" ","")
    aln2_ref             = "ATCGAT-CGAT".replace(" ","")
    qual2                =     "AAAAAA ".replace(" ", "")

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "NNCGATCCGAN"
        assert ref_seq ==      "ATCGAT-CGAT"
        assert score == 63.636
        assert not caching_ok

    # insertion in r1 and r2, different positions
    qual1                =   "AAAAAA"
    aln1_seq             = "--CGATCC---".replace(" ","")
    aln1_ref             = "ATCGAT-CGAT".replace(" ","")
    aln2_seq             = "----ATATGA-".replace(" ","")
    aln2_ref             = "ATCGATC-GAT".replace(" ","")
    qual2                =     "AAAAAA ".replace(" ", "")

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "NNCGATCCTGAN"
        assert ref_seq ==      "ATCGAT-C-GAT"
        assert score == 58.333
        assert not caching_ok

    # insertion at beginning of r1
    qual1                = "AAAAA"
    aln1_seq             = "TA-CGA----- ".replace(" ","")
    aln1_ref             = "-ATCGATCGAT ".replace(" ","")
    aln2_seq             = " --------AT ".replace(" ","")
    aln2_ref             = " ATCGATCGAT".replace(" ","")
    qual2                =          "AA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "A-CGANNNAT"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 60.0
        assert caching_ok

    # insertion at end of r2 and beginning of r1
    qual1                = "AAAAA"
    aln1_seq             = "TA-CGA-----   ".replace(" ","")
    aln1_ref             = "-ATCGATCGAT   ".replace(" ","")
    aln2_seq             = " -----TCGATCCA".replace(" ","")
    aln2_ref             = " ATCGATCGAT---".replace(" ","")
    qual2                =       "AAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq ==      "A-CGATCGAT"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 90.0
        assert caching_ok

    qual1    = '>1>1A@DFAADAGGGGGGGGGGHHHHHHHHHHHHHHHGGHHHHHHHHGGGHHHHHHHHHGHHHHHHHHHHHHGHGGGGGGGGGGHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHGHHGGGGGHHHHG'
    aln1_seq = 'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGCGGC-----CTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT'
    aln1_ref = 'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT'
    aln2_seq = 'AACCACAGCC-----GAGCCTCTTGAAGCCATTCTTACAGATGATGAAC-CAGG--CCGCGTTGGGAGCTCCAGAAGGGGATCATGACCT----CCTCACCTGTGGGCAGTGCCAGATGAACTTCCCATTGGGGGACATT'
    aln2_ref = 'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGC--CGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT-----'
    qual2    = 'BCCDCCDFDDDDGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHGHGGGGGGGHGHGGHHHHHHHHHGGGGGHHHHHHHHHHHHHHHHHHHGGHGHHHHHHHHHHHHHHHHHHHHHHHGGGGGHH'

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(
        aln1_seq,
        aln1_ref,
        calc_score(aln1_seq, aln1_ref),
        qual1,
        aln2_seq,
        aln2_ref,
        calc_score(aln2_seq, aln2_ref),
        qual2,
    )
    with check:
        assert aln_seq == 'AACCACCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAACTCATGATCCCCTTCTGGAGCTCCAAAAGGGGATCATGACCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGTTCCCATTTGGTT'
        assert ref_seq == 'AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGC--CGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT'
        assert score == 86.667
        assert not caching_ok


    # alternating qualities
    qual1                = "BABABABABA"
    aln1_seq             = "ACCAACCAAT".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "TTGGTTGGTT".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "ABABABABAB"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1, aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2)
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)
    check.is_false(caching_ok) #TODO: Should be false?

    # large insertion in r1
    qual1                = "AAAAAA"
    aln1_seq             = "ACGTGA---------".replace(" ","")
    aln1_ref             = "A-----TCGATCGAT".replace(" ","")
    aln2_seq             = "------CGAT".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =       "AAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1, aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2)
    check.equal(aln_seq, "ACGTGANNNNNCGAT")
    check.equal(ref_seq, "A-----TCGATCGAT")
    check.equal(score, 33.333)
    check.is_true(caching_ok)

    # large insertion in r2
    qual1                = "AAAAA"
    aln1_seq             = "ATCGA-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-----TTAGCT---".replace(" ","")
    aln2_ref             = "ATCGAT---C-GAT".replace(" ","")
    qual2                =      "AAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1, aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2)
    check.equal(aln_seq, "ATCGATTAGCTNNN")
    check.equal(ref_seq, "ATCGAT---C-GAT")
    check.equal(score, 50)
    check.is_true(caching_ok)

    # Conflicts with reference
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "TAGCTAGCTA".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "TAGCTAGCTA".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "AAAAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1, aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2)
    check.equal(aln_seq, "TAGCTAGCTA")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 0)
    check.is_true(caching_ok)

    # Conflicts between reads
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "TAGCTAGCTA".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "ATCGATCGAT".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "AAAAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1, aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2)
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)
    check.is_false(caching_ok) #TODO: Should this be false?

    # Alternating reads
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "AT--AT--AT".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "--CG--CG--".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "AAAAAAAAAA"

    aln_seq, aln_qual, ref_seq, score, caching_ok = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, calc_score(aln1_seq, aln1_ref), qual1, aln2_seq, aln2_ref, calc_score(aln2_seq, aln2_ref), qual2)
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)
    check.is_true(caching_ok)


def test_split_quant_window_coordinates_single():
    assert [(5, 10)] == CRISPRessoCORE.split_quant_window_coordinates('5-10')


def test_split_quant_window_coordinates_multiple():
    assert CRISPRessoCORE.split_quant_window_coordinates('2-5_10-12') == [(2, 5), (10, 12)]


def test_split_quant_window_coordinates_error():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('a-5')


def test_split_quant_window_coordinates_empty():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('_')


def test_split_quant_window_coordinates_partially_empty():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('1-3_')


def test_split_quant_window_coordinates_blank():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('')


def test_get_include_idxs_from_quant_window_coordinates():
    quant_window_coordinates = '1-10_12-20'
    assert CRISPRessoCORE.get_include_idxs_from_quant_window_coordinates(quant_window_coordinates) == [*list(range(1, 11)), *list(range(12, 21))]


def test_get_cloned_include_idxs_from_quant_window_coordinates():
    quant_window_coordinates = '1-10_12-20'
    ref = 'TTACCGAGTGCACAAGTGCACGT'
    aln = 'TTACCGAGTGCACAAGTGCACGT'
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 11)), *list(range(12, 21))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_beginning():
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp insertion at the beginning (left)
    # Ind:                1111111111222
    #           01234567890123456789012
    # QWC:       |        | |       |
    ref = '-----TTACCGAGTGCACAAGTGCACGT'
    aln = 'AAGGTTTACCGAGTGCACAAGTGCACGT'
    # QWC:       |        | |       |
    # Ind:           111111111122222222
    #      0123456789012345678901234567
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(6, 16)), *list(range(17, 26))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_beginning():
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp deletion at the beginning (left)
    # Ind:           1111111111222
    #      01234567890123456789012
    # QWC:  |        | |       |
    ref = 'TTACCGAGTGCACAAGTGCACGT'
    aln = '------AGTGCACAAGTGCACGT'
    # QWC:       |   | |       |
    # Ind:                 1111111
    #            01234567890123456
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(0, 5)), *list(range(6, 15))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion():
    quant_window_coordinates = '10-20_35-40'
    ref = 'A' * 23 + 'T' * 7 + 'G' * 30
    aln = 'A' * 23 + '-' * 7 + 'G' * 30
    # represents a 7bp deletion in the middle
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(10, 21)), *list(range(35-7, 41-7))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_modified():
    quant_window_coordinates = '10-25_35-40'
    # represents a 7bp deletion in the middle, where part of the QW is deleted
    # [0, 1, 3, 4, ... , 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 24, ... , 33]
    ref = 'A' * 23 + 'T' * 7 + 'G' * 30
    aln = 'A' * 23 + '-' * 7 + 'G' * 30
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(10, 23)), *list(range(35-7, 41-7))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_end_modified():
    # 5 bp deletion at end of 20 bp sequence
    quant_window_coordinates = '1-5_10-20'
    # Ind:           11111111112
    #      012345678901234567890
    # QWC:  |   |    |         |
    ref = 'AAAAAAAAAAAAAAAATTTTT'
    aln = 'AAAAAAAAAAAAAAAA-----'
    # QWC:  |   |    |    |
    # Ind:           111111
    #      0123456789012345
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 6)), *list(range(10, 16))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion():
    # 5 bp deletion and 5 bp insertion
    quant_window_coordinates = '1-5_10-18'
    # Ind:           1111     11111
    #      01234567890123     45678
    # QWC:  |   |    |            |
    ref = 'AAAAATTTTTGGGG-----AAAAA'
    aln = 'AAAAA-----GGGGCCCCCAAAAA'
    # QWC:  |                     |
    # Ind:                111111111
    #      01234     56789012345678
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 19))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion_modified():
    quant_window_coordinates = '1-5_10-20'
    # Ind:                 11111111112
    #      012 34567     8901234567890
    # QWC:  |    |         |         |
    ref = 'AAA-CCCCC-----AAATTTTCCCCCC'
    aln = 'AAAT-CCCCGGGGGAAA----CCCCCC'
    # QWC:  |    |         |         |
    # Ind:           1111111    111122
    #      0123 456789012345    678901
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*[1, 2, 3, 4, 5], *[15, 16, 17, 18, 19, 20, 21]]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion():
    quant_window_coordinates = '2-7'

    # Ind: 0123  456789
    # QWC:   |      |
    ref = 'AAAA--TTTTTT'
    aln = 'AAAAGGTTTTTT'
    # QWC:   |      |
    # Ind:           11
    #      012345678901
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(2, 10))


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_start():
    quant_window_coordinates = '1-3'

    # Ind: 0    123456
    # QWC:      | |
    ref = 'T----ACTGT'
    aln = 'TTCCCACTGT'
    # QWC:      | |
    # Ind: 0123456789
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [5, 6, 7]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_end():
    quant_window_coordinates = '1-4'

    # Ind: 0123    45678
    # QWC:  |      |
    ref = 'GGGT----ACTGT'
    aln = 'GGGTTCCCACTGT'
    # QWC:  |      |
    # Ind:           111
    #      0123456789012
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(1, 9))


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_before_qw():
    quant_window_coordinates = '6-9'

    # Ind:           11
    #      012345678901
    # QWC:       |  |
    ref = 'GGGTACTGTCCA'
    aln = 'GGGTA-TGTCCA'
    # QWC:       |  |
    # Ind:            1
    #      01234 567890
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(5, 9))


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_before_qw():
    quant_window_coordinates = '6-9'
    # Ind:            1
    #      01234 567890
    # QWC:        |  |
    ref = 'GGGTA-TGTCCA'
    aln = 'GGGTACTGTCCA'
    # QWC:        |  |
    # Ind:           11
    #      012345678901
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(7, 11))


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_deletion_outside_qw():
    quant_window_coordinates = '4-6_11-14'

    # Ind:               11111
    #      0123    45678901234
    # QWC:         | |    |  |
    ref = 'GGGT----ACTTTTGTCCA'
    aln = 'GGGTTCCCACT---GTCCA'
    # QWC:         | |    |  |
    # Ind:           1   11111
    #      01234567890   12345
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [8, 9, 10, 12, 13, 14, 15]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_across_qw():
    # 6 bp insertion in middle of 4 bp sequence
    quant_window_coordinates = '1-3'
    # Ind: 01      23
    # QWC:  |       |
    ref = 'AA------TT'
    aln = 'AAGGGGGGTT'
    # QWC:  |       |
    # Ind: 0123456789
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == list(range(1, 10))


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_overlap_start():
    quant_window_coordinates = '2-5'
    # Ind: 012345
    # QWC:   |  |
    ref = 'AATTTT'
    aln = 'A--TTT'
    # QWC:    | |
    # Ind: 0  123
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2, 3]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_overlap_end():
    quant_window_coordinates = '1-3'
    # Ind: 012345
    # QWC:  | |
    ref = 'AATTTT'
    aln = 'AAT---'
    # QWC:  ||
    # Ind: 012
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_overlap_end_single_bp():
    quant_window_coordinates = '2-3'
    # Ind: 012345
    # QWC:   ||
    ref = 'AATTTT'
    aln = 'AAT---'
    # QWC:   |
    # Ind: 012
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [2]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_entire_qw():
    # 5 bp deletion of entire qw
    quant_window_coordinates = '1-4_7-10'
    ref = 'AAAAAATTTT'
    aln = 'AAAAAA----'
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2, 3, 4]


def test_get_cloned_include_idxs_from_quant_window_coordinates_include_zero():
    quant_window_coordinates = '0-4'
    ref = 'AAAAA'
    aln = 'AAAAA'
    s1inds, _ = CRISPRessoShared.get_relative_coordinates(ref, aln)
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [0, 1, 2, 3, 4]


# Testing parallelization functions
def test_regular_input():
    # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(100, 4) == [0, 25, 50, 75, 100]

def test_remainder_input():
#     # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(101, 4) == [0, 25, 50, 75, 101]

def test_similar_num_reads_input():
#     # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(11, 10) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]

def test_large_similar_num_reads_input():
#     # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(101, 100) == list(range(0, 100)) + [101]

def test_more_processes_than_reads():
#     # Test with typical input
    # assert CRISPRessoCORE.get_variant_cache_equal_boundaries(3, 5)
    # assert that an exception is raised
    with pytest.raises(Exception):
        CRISPRessoCORE.get_variant_cache_equal_boundaries(3, 5)

def test_single_process():
    # Test with a single process
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(50, 1) == [0, 50]

def test_zero_sequences():
    # Test with zero unique sequences
    with pytest.raises(Exception):
        CRISPRessoCORE.get_variant_cache_equal_boundaries(0, 3)

def test_large_numbers():
    # Test with large number of processes and sequences
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(10000, 10)
    assert len(boundaries) == 11  # Check that there are 11 boundaries

def test_sublist_generation():
    n_processes = 4
    unique_reads = 100
    mock_variant_cache = [i for i in range(unique_reads)]
    assert len(mock_variant_cache) == 100
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(unique_reads, n_processes)
    assert boundaries == [0, 25, 50, 75, 100]
    sublists = []
    for i in range(n_processes):
        left_sublist_index = boundaries[i]
        right_sublist_index = boundaries[i+1]
        sublist = mock_variant_cache[left_sublist_index:right_sublist_index]
        sublists.append(sublist)
    assert [len(sublist) for sublist in sublists] == [25, 25, 25, 25]
    assert [s for sublist in sublists for s in sublist] == mock_variant_cache

def test_irregular_sublist_generation():
    n_processes = 4
    unique_reads = 113
    mock_variant_cache = [i for i in range(unique_reads)]
    assert len(mock_variant_cache) == 113
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(unique_reads, n_processes)
    # assert boundaries == [0, 25, 50, 75, 100]
    sublists = []
    for i in range(n_processes):
        left_sublist_index = boundaries[i]
        right_sublist_index = boundaries[i+1]
        sublist = mock_variant_cache[left_sublist_index:right_sublist_index]
        sublists.append(sublist)
    assert [len(sublist) for sublist in sublists] == [28,28,28,29]
    assert [s for sublist in sublists for s in sublist] == mock_variant_cache


# Test upset plot data functions
def test_get_base_edit_target_sequence():
    df_alleles = pd.read_csv('tests/df_alleles.txt')
    ref_seq = 'CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG'
    base_editor_target_ref_skip_allele_count = 0

    target_seq = CRISPRessoCORE.get_base_edit_target_sequence(
        ref_seq,
        df_alleles,
        base_editor_target_ref_skip_allele_count
    )

    assert target_seq == 'AATACGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCCGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGGCATGGCCCCATTCGCACGGCTCTGGAGCGGC'


def test_get_bp_subs_one_sub():
    ref_seq = 'AAA'
    aln_ref_seq = 'AAA'
    aln_target_seq = 'CAA'
    ref_positions_to_include = [0, 1, 2]
    ref_changes_dict = CRISPRessoCORE.get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include)

    assert len(bp_substitutions_arr) == 1
    assert bp_substitutions_arr[0][0] == 0
    assert bp_substitutions_arr[0][1] == 'A'
    assert bp_substitutions_arr[0][2] == 'C'

def test_get_bp_subs_two_subs():
    ref_seq = 'AAA'
    aln_ref_seq = 'AAA'
    aln_target_seq = 'CCA'
    ref_positions_to_include = [0, 1, 2]
    ref_changes_dict = CRISPRessoCORE.get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include)

    assert len(bp_substitutions_arr) == 2
    assert bp_substitutions_arr[0] == (0, 'A', 'C')
    assert bp_substitutions_arr[1] == (1, 'A', 'C')

def test_get_bp_subs_insertions():

    ref_seq = 'AAAAAA'
    aln_ref_seq = 'AAA-AAA-----'
    aln_target_seq = 'AAACAAACCCCC'
    ref_positions_to_include = [0, 1, 2, 3, 4, 5]
    ref_changes_dict = CRISPRessoCORE.get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = CRISPRessoCORE.get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include)

    assert len(bp_substitutions_arr) == 2
    assert bp_substitutions_arr[0] == (2, 'A', 'AC')
    assert bp_substitutions_arr[1] == (5, 'A', 'ACCCCC')


def test_get_base_edit_target_sequence():

    df_alleles = pd.read_csv('tests/test_be_df.txt')
    ref_seq = 'AAAA'
    base_editor_target_ref_skip_allele_count = 0

    target_seq = CRISPRessoCORE.get_base_edit_target_sequence(
        ref_seq,
        df_alleles,
        base_editor_target_ref_skip_allele_count
    )

    assert target_seq == 'AAGA'



def test_get_upset_plot_counts():
    df_alleles = pd.read_csv('tests/test_be_df.txt')
    target_seq = 'AAGA'
    bp_substitutions_arr = [(3, 'A', 'G')]

    wt_ref_name = 'TEST'

    counts_dict = CRISPRessoCORE.get_upset_plot_counts(
        df_alleles,
        bp_substitutions_arr,
        wt_ref_name
    )

    assert len(counts_dict) == 19
    assert counts_dict['total_alleles'] == 3
    assert counts_dict['total_alleles_reads'] == 100
    assert sum(counts_dict['binary_allele_counts'].values()) == 100



def test_write_base_edit_counts():


    OUTPUT_DIRECTORY = '.'
    clean_file_prefix = ""
    _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, clean_file_prefix + filename)
    ref_name = 'TEST'
    bp_substitutions_arr = [(3, 'A', 'G')]
    counts_dict = CRISPRessoCORE.get_upset_plot_counts(
        pd.read_csv('tests/test_be_df.txt'),
        bp_substitutions_arr,
        ref_name,
    )

    expected_files = [
        '10i.TEST.arrays.txt',
        '10i.TEST.binary_allele_counts.txt',
        '10i.TEST.category_allele_counts.txt',
        '10i.TEST.counts.txt',
        '10i.TEST.precise_allele_counts.txt',
    ]


    CRISPRessoCORE.write_base_edit_counts(
        ref_name,
        counts_dict,
        bp_substitutions_arr,
        _jp
    )

    for filename in expected_files:
        if os.path.exists(filename):
            os.remove(filename)
        else:
            assert False


if __name__ == "__main__":
    # execute only if run as a script
    test_get_consensus_alignment_from_pairs()
