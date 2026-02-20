"""Unit tests for CRISPResso2Align."""

from CRISPResso2 import CRISPResso2Align
import numpy as np

ALN_MATRIX = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")
AA_MATRIX = CRISPResso2Align.read_matrix("./CRISPResso2/BLOSUM62")


# =============================================================================
# Tests for read_matrix
# =============================================================================


def test_read_matrix_ednafull():
    """Test reading EDNAFULL matrix returns correct scores."""
    matrix = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")
    # Match score for A-A should be 5, mismatch A-T should be -4
    assert matrix[ord('A'), ord('A')] == 5
    assert matrix[ord('A'), ord('T')] == -4


def test_read_matrix_blosum62():
    """Test reading BLOSUM62 matrix returns correct scores."""
    matrix = CRISPResso2Align.read_matrix("./CRISPResso2/BLOSUM62")
    # Self-score for M should be 5 in BLOSUM62
    assert matrix[ord('M'), ord('M')] == 5


# =============================================================================
# Tests for global_align - basic alignment
# =============================================================================


def test_global_align():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "ATTA"
    assert seq2 == "ATTA"
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


def test_global_align_BLOSUM62():
    """Test amino acid alignment with BLOSUM62 matrix."""
    seq1, seq2, score = CRISPResso2Align.global_align("MRWY", "MRWY", matrix=AA_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "MRWY"
    assert seq2 == "MRWY"
    assert score == 100


def test_global_align_score_perfect_match():
    """Test that perfect match gives score of 100."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "ATCGATCG", "ATCGATCG", matrix=ALN_MATRIX,
        gap_incentive=np.array([0] * 9, dtype=int)
    )
    assert score == 100


def test_global_align_symmetry():
    """Test that alignment is symmetric (swapping inputs gives same alignment)."""
    seq1_a, seq2_a, score_a = CRISPResso2Align.global_align(
        "ATCG", "ATCC", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    seq1_b, seq2_b, score_b = CRISPResso2Align.global_align(
        "ATCC", "ATCG", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    # Scores should be equal
    assert score_a == score_b


# =============================================================================
# Tests for global_align - gap incentives with BLOSUM62
# =============================================================================


def test_global_align_gap_incentive_pos1_BLOSUM62():
    """Test amino acid alignment with gap incentive at position 1."""
    # Use MRRRY vs MRRY - repeated R allows gap at different positions
    seq1, seq2, score = CRISPResso2Align.global_align("MRRRY", "MRRY", matrix=AA_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "MRRRY"
    assert seq2 == "M-RRY"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


def test_global_align_gap_incentive_pos2_BLOSUM62():
    """Test amino acid alignment with gap incentive at position 2."""
    # Use MRRRY vs MRRY - repeated R allows gap at different positions
    seq1, seq2, score = CRISPResso2Align.global_align("MRRRY", "MRRY", matrix=AA_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "MRRRY"
    assert seq2 == "MR-RY"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


def test_global_align_gap_incentive_pos3_BLOSUM62():
    """Test amino acid alignment with gap incentive at position 3."""
    # Use MRRRY vs MRRY - repeated R allows gap at different positions
    seq1, seq2, score = CRISPResso2Align.global_align("MRRRY", "MRRY", matrix=AA_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "MRRRY"
    assert seq2 == "MRR-Y"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


# =============================================================================
# Tests for global_align - gap incentives in sequence 1
# =============================================================================


def test_global_align_gap_incentive_s1():
    """Test the global_align gap incentives for gaps in sequence 1 (the first sequence)."""
    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0, 0], dtype=int))
    #    print('seq1: ' + seq1 + ' seq2: ' + seq2 + ' score ' + str(score))
    assert seq1 == "ATTTA"
    assert seq2 == "ATTTA"
    assert score == 100

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "A-TTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "AT-TA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TTTT-"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "-TTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "T-TTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TT-TT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TTT-T"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TTTT-"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


# =============================================================================
# Tests for global_align - gap incentives in sequence 2
# =============================================================================


def test_global_align_gap_incentive_s2():
    """Test the global_align gap incentives for gaps in sequence 2 (the second sequence)."""
    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "A-TTA"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "AT-TA"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "-TTTT"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "T-TTT"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "TT-TT"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "TTT-T"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "TTTT-"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "TTTT-"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


# =============================================================================
# Tests for global_align - edge cases
# =============================================================================


def test_global_align_single_base():
    """Test alignment with single base sequences."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "A", "A", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0], dtype=int)
    )
    assert seq1 == "A"
    assert seq2 == "A"
    assert score == 100


def test_global_align_with_n():
    """Test alignment with N bases."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "ANNG", "ATCG", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    assert seq1 == "A-NNG"
    assert seq2 == "ATC-G"
    assert score == 40.0


def test_global_align_all_n():
    """Test alignment with all N bases."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "NNNN", "NNNN", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    assert seq1 == "NNNN"
    assert seq2 == "NNNN"


def test_global_align_repeated_base():
    """Test alignment with repeated bases (homopolymer)."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "AAAAAAAAAA", "AAAAAAAAAA", matrix=ALN_MATRIX,
        gap_incentive=np.array([0] * 11, dtype=int)
    )
    assert seq1 == "AAAAAAAAAA"
    assert seq2 == "AAAAAAAAAA"
    assert score == 100


def test_global_align_completely_different():
    """Test alignment with completely different sequences."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "AAAA", "TTTT", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    assert seq1 == "---AAAA"
    assert seq2 == "TTTT---"
    assert score == 0.0



if __name__ == "__main__":
    # execute only if run as a script
    test_global_align()
    test_global_align_gap_incentive_s1()
    test_global_align_gap_incentive_s2()
