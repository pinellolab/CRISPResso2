"""Unit tests for CRISPResso2CORE."""
import pytest
from pytest_check import check

from CRISPResso2 import CRISPRessoCORE, CRISPRessoShared

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
    s1inds = list(range(22))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 11)), *list(range(12, 21))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_beginning():
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp insertion at the beginning (left)
    s1inds = list(range(5, 27))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(6, 16)), *list(range(17, 26))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_beginning():
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp deletion at the beginning (left)
    s1inds = [-1, -1, -1, -1, -1 ] + list(range(26))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(0, 6)), *list(range(7, 16))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion():
    quant_window_coordinates = '10-20_35-40'
    # represents a 7bp deletion in the middle
    s1inds = list(range(23)) + [22, 22, 22, 22, 22, 22, 22] + list(range(23, 34))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(10, 21)), *list(range(35-7, 41-7))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_modified():
    quant_window_coordinates = '10-25_35-40'
    # represents a 7bp deletion in the middle, where part of the QW is deleted
    # [0, 1, 3, 4, ... , 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 24, ... , 33]
    s1inds = list(range(23)) + [22, 22, 22, 22, 22, 22, 22] + list(range(23, 34))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(10, 23)), *list(range(35-7, 41-7))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_end_modified(): 
    # 5 bp deletion at end of 20 bp sequence
    quant_window_coordinates = '1-5_10-20'
    s1inds = [*list(range(16)), *[15, 15, 15, 15, 15]]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 6)), *list(range(10, 16))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion():
    # 5 bp deletion and 5 bp insertion
    quant_window_coordinates = '1-5_10-20'
    s1inds = [0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 6, 7, 8, 9, 15, 16, 17, 18, 19, 20]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 6)), *[6, 7, 8, 9, 15, 16, 17, 18, 19, 20]]

def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion_modified():
    quant_window_coordinates = '1-5_10-20'
    s1inds = [0, 1, 2, 2, 4, 5, 6, 7, 7, 7, 7, 7, 7, 8, 9, 10, 15, 16, 17, 18, 19]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*[1,2,4,5], *[8, 9, 10, 15, 16, 17, 18, 19]]

def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_across_qw():
    # 6 bp insertion in middle of 4 bp sequence
    quant_window_coordinates = '1-4'
    s1inds = [0,1,2,9,10]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1,2,9,10]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_entire_qw():
    # 5 bp deletion of entire qw
    quant_window_coordinates = '1-4_7-10'
    s1inds = [0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2, 3, 4]

def test_get_cloned_include_idxs_from_quant_window_coordinates_include_zero():
    quant_window_coordinates = '0-5'
    s1inds = [0, 1, 2, 3, 4, 5]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [0, 1, 2, 3, 4, 5]


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



    
if __name__ == "__main__":
# execute only if run as a script
    test_get_consensus_alignment_from_pairs()
