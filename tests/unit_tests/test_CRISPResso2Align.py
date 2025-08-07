"""Unit tests for CRISPResso2Align."""

from CRISPResso2 import CRISPResso2Align
import numpy as np

ALN_MATRIX = CRISPResso2Align.read_matrix('./CRISPResso2/EDNAFULL')
AA_MATRIX = CRISPResso2Align.read_matrix('./CRISPResso2/BLOSUM62')


def test_global_align():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,0],dtype=int))
    assert seq1 == 'ATTA'
    assert seq2 == 'ATTA'
    assert score == 100

def test_global_align_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,0,0,0,0],dtype=int))
    assert seq1 == 'ATTA'
    assert seq2 == 'ATTA'
    assert score == 100

def test_global_align_gap_incentive_pos1_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,1,0,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'A-TTA'
    assert round(score,3) == round(100*4/5.0,3)

def test_global_align_gap_incentive_pos2_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,0,1,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'AT-TA'
    assert round(score,3) == round(100*4/5.0,3)

def test_global_align_gap_incentive_pos3_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,0,0,1,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATT-A'
    assert round(score,3) == round(100*4/5.0,3)


def test_global_align_gap_incentive_s1():
    """Test the global_align gap incentives for gaps in sequence 1 (the first sequence)."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,0,0],dtype=int))
#    print('seq1: ' + seq1 + ' seq2: ' + seq2 + ' score ' + str(score))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATTTA'
    assert score == 100

    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([1,0,0,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATT-A'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,1,0,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'A-TTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,1,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'AT-TA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,1,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATT-A'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,1,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATT-A'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,1],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATT-A'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTTT', 'TTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,0],dtype=int))
    assert seq1 == 'TTTTT'
    assert seq2 == 'TTTT-'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTTT', 'TTTT', matrix=ALN_MATRIX, gap_incentive=np.array([1,0,0,0,0],dtype=int))
    assert seq1 == 'TTTTT'
    assert seq2 == '-TTTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTTT', 'TTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,1,0,0,0],dtype=int))
    assert seq1 == 'TTTTT'
    assert seq2 == 'T-TTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTTT', 'TTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,1,0,0],dtype=int))
    assert seq1 == 'TTTTT'
    assert seq2 == 'TT-TT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTTT', 'TTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,1,0],dtype=int))
    assert seq1 == 'TTTTT'
    assert seq2 == 'TTT-T'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTTT', 'TTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,1],dtype=int))
    assert seq1 == 'TTTTT'
    assert seq2 == 'TTTT-'
    assert round(score,3) == round(100*4/5.0,3)

def test_global_align_gap_incentive_s2():
    """Test the global_align gap incentives for gaps in sequence 2 (the second sequence)."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([1,0,0,0,0,0],dtype=int))
    assert seq1 == 'ATT-A'
    assert seq2 == 'ATTTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,1,0,0,0,0],dtype=int))
    assert seq1 == 'A-TTA'
    assert seq2 == 'ATTTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,1,0,0,0],dtype=int))
    assert seq1 == 'AT-TA'
    assert seq2 == 'ATTTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,1,0,0],dtype=int))
    assert seq1 == 'ATT-A'
    assert seq2 == 'ATTTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,1,0],dtype=int))
    assert seq1 == 'ATT-A'
    assert seq2 == 'ATTTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTTA', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,0,1],dtype=int))
    assert seq1 == 'ATT-A'
    assert seq2 == 'ATTTA'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTT', 'TTTTT', matrix=ALN_MATRIX, gap_incentive=np.array([1,0,0,0,0,0],dtype=int))
    assert seq1 == '-TTTT'
    assert seq2 == 'TTTTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTT', 'TTTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,1,0,0,0,0],dtype=int))
    assert seq1 == 'T-TTT'
    assert seq2 == 'TTTTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTT', 'TTTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,1,0,0,0],dtype=int))
    assert seq1 == 'TT-TT'
    assert seq2 == 'TTTTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTT', 'TTTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,1,0,0],dtype=int))
    assert seq1 == 'TTT-T'
    assert seq2 == 'TTTTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTT', 'TTTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,1,0],dtype=int))
    assert seq1 == 'TTTT-'
    assert seq2 == 'TTTTT'
    assert round(score,3) == round(100*4/5.0,3)

    seq1, seq2, score = CRISPResso2Align.global_align('TTTT', 'TTTTT', matrix=ALN_MATRIX, gap_incentive=np.array([0,0,0,0,0,1],dtype=int))
    assert seq1 == 'TTTT-'
    assert seq2 == 'TTTTT'
    assert round(score,3) == round(100*4/5.0,3)


if __name__ == "__main__":
    # execute only if run as a script
    test_global_align()
    test_global_align_gap_incentive_s1()
    test_global_align_gap_incentive_s2()
