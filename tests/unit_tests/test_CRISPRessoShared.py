from CRISPResso2 import CRISPResso2Align, CRISPRessoShared

ALN_MATRIX = CRISPResso2Align.read_matrix('./CRISPResso2/EDNAFULL')


def test_get_mismatches():
    mismatch_cords = CRISPRessoShared.get_mismatches(
        'ATTA',
        'ATTA',
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 0

    mismatch_cords = CRISPRessoShared.get_mismatches(
        'GCAGTGGGCGCGCTA',
        'CCCACTGAAGGCCC',
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 6 

def test_get_relative_coordinates():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATCGT', 'TTCGT')
    assert s1inds_gap_left == [0, 1, 2, 3, 4]
    assert s1inds_gap_right == [0, 1, 2, 3, 4]


def test_get_relative_coordinates_to_gap():
    # unaligned sequences
    seq_1 = 'TTCGT'
    seq_2 = 'TTCT'

    # aligned_sequences
    to_sequence = 'TTC-T'
    from_sequence = 'TTCGT'

    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates(to_sequence, from_sequence)
    assert s1inds_gap_left == [0, 1, 2, 2, 3]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]


    assert seq_1[0] == seq_2[s1inds_gap_left[0]]
    assert seq_1[1] == seq_2[s1inds_gap_left[1]]
    assert seq_1[2] == seq_2[s1inds_gap_left[2]]
    assert seq_1[4] == seq_2[s1inds_gap_left[4]]


def test_get_relative_coordinates_start_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('--CGT', 'TTCGT')
    assert s1inds_gap_left == [-1, -1, 0, 1, 2]
    assert s1inds_gap_right == [0, 0, 0, 1, 2]


def test_get_relative_coordinates_from_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATCGT', 'ATC-T')
    assert s1inds_gap_left == [0, 1, 2, 4]
    assert s1inds_gap_right == [0, 1, 2, 4]

def test_get_relative_coordinates_end_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATC--', 'ATCGT')
    assert s1inds_gap_left == [0, 1, 2, 2, 2]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]

def test_get_relative_coordinates_multiple_gaps():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('AT--TC--G--CC', 'ATCGTCGCGTTCC')
    assert s1inds_gap_left == [0, 1, 1, 1, 2, 3, 3, 3, 4, 4, 4, 5, 6]
    assert s1inds_gap_right == [0, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 6]

def test_get_relative_coordinates_ind_and_dels():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATG--C', 'A-GCTC')
    assert s1inds_gap_left == [0, 2, 2, 2, 3]
    assert s1inds_gap_right == [0, 2, 3, 3, 3]

