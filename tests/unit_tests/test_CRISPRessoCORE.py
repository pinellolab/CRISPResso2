"""Unit tests for CRISPResso2CORE."""
import pytest

from CRISPResso2 import CRISPRessoCORE

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
    assert score == 40 #double check this score... should be 4/10

    # deletion in r2
    qual1                =   "AAA"
    aln1_seq             = "--CGA-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-----T-GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =      "AAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGAT-GAN"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 60 #double check this score... should be 6/10

    # insertion at beginning of r1
    qual1                = "AAAA"
    aln1_seq             = "TA-CGA----- ".replace(" ","")
    aln1_ref             = "-ATCGATCGAT ".replace(" ","")
    aln2_seq             = " --------AT ".replace(" ","")
    aln2_ref             = " ATCGATCGAT".replace(" ","")
    qual2                =      "AAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "TA-CGANNNAT"
    assert ref_seq ==      "-ATCGATCGAT"
    assert score == 54 #double check this score... should be 6/11

if __name__ == "__main__":
# execute only if run as a script
    test_get_consensus_alignment_from_pairs()
