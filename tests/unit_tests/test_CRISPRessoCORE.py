"""Unit tests for CRISPRessoCORE."""
import pytest

from CRISPResso2 import CRISPRessoCORE, CRISPRessoShared

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

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGATCGAT"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 80

    #test quality difference
    qual1                =   "AAAB"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGATCGAT"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 80

    #test quality difference
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "BAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    print('got aln_seq ' + str(aln_seq))
    assert aln_seq ==      "NNCGAGCGAT"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 70

    #gaps between r1 and r2
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "--CGA-----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-------GA-"
    aln2_ref             = "ATCGATCGAT"
    qual2                = "AAAAAAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGANNGAN"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 50

    print('Finished easy tests... now for the hard stuff')

    #insertion in r1
    qual1                =   "AAAA"
    aln1_seq             = "--CCGA-----".replace(" ","") #added replace for vertical alignment
    aln1_ref             = "ATC-GATCGAT".replace(" ","")
    aln2_seq             = "--- ----GA-".replace(" ","")
    aln2_ref             = "ATC GATCGAT".replace(" ","")
    qual2                =         "AA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCCGANNGAN"
    assert ref_seq ==      "ATC-GATCGAT"
    assert score == 45 #double check this score... should be 5/11

    #deletion in r1
    qual1                =   "AA"
    aln1_seq             = "--C-A-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-------GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =        "AA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNC-ANNGAN"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 50 #double check this score... should be 5/10


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
    
if __name__ == "__main__":
# execute only if run as a script
    test_get_consensus_alignment_from_pairs()
