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


def test_get_consensus_alignment_from_pairs():
    """Tests for generating consensus alignments from paired reads."""
    try:
        CRISPRessoCORE.get_consensus_alignment_from_pairs
    except AttributeError:
        pytest.xfail('get_consensus_alignment_from_pairs is not implemented yet!')

    print("testing Easy")

    #basic test
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----TCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

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
        assert aln_seq ==      "NNCGATCGAT"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 80
        assert caching_ok

    #test quality difference
    qual1                =   "AAAB"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

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
        assert aln_seq ==      "NNCGATCGAT"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 80
        assert not caching_ok

    #test quality difference
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "BAAAA"

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
        assert aln_seq ==      "NNCGAGCGAT"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 70
        assert not caching_ok

    #gaps between r1 and r2
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "--CGA-----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-------GA-"
    aln2_ref             = "ATCGATCGAT"
    qual2                = "AAAAAAAAAA"

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
        assert aln_seq ==      "NNCGANNGAN"
        assert ref_seq ==      "ATCGATCGAT"
        assert score == 50
        assert caching_ok

    print('Finished easy tests... now for the hard stuff')

    #insertion in r1
    qual1                =   "AAAA"
    aln1_seq             = "--CCGA-----".replace(" ","") #added replace for vertical alignment
    aln1_ref             = "ATC-GATCGAT".replace(" ","")
    aln2_seq             = "--- ----GA-".replace(" ","")
    aln2_ref             = "ATC GATCGAT".replace(" ","")
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
        assert ref_seq ==      "ATC-GATCGAT"
        assert score == 45.455
        assert caching_ok

    #deletion in r1
    qual1                =   "AA"
    aln1_seq             = "--C-A-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-------GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =        "AA"

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


def test_long_coding_seq_filename_not_too_long():
    """Verify that a long coding_seq doesn't produce filenames exceeding OS limits.

    This is a regression test for the issue where coding sequences >200bp
    were embedded verbatim in filenames, exceeding the 255-byte limit.
    """
    long_coding_seq = 'ATCGATCG' * 50  # 400 characters

    coding_seq_short = CRISPRessoShared._get_short_seq_id(long_coding_seq)

    # Simulate the two filename patterns from CRISPRessoCORE.py
    ref_plot_name = 'Reference.'  # typical ref_plot_name
    fig_filename = '9a.' + ref_plot_name + 'amino_acid_table_around_' + coding_seq_short + '.pdf'
    table_filename = ref_plot_name + 'amino_acid_table_for_' + coding_seq_short + '.txt'

    # 255 bytes is the typical max filename length on macOS/Linux
    assert len(fig_filename.encode('utf-8')) <= 255, \
        f"Figure filename too long: {len(fig_filename.encode('utf-8'))} bytes"
    assert len(table_filename.encode('utf-8')) <= 255, \
        f"Table filename too long: {len(table_filename.encode('utf-8'))} bytes"

    # The short id should be much shorter than the original
    assert len(coding_seq_short) < len(long_coding_seq)
    assert len(coding_seq_short) == 17  # 8 prefix + '_' + 8 hash


if __name__ == "__main__":
# execute only if run as a script
    test_get_consensus_alignment_from_pairs()
