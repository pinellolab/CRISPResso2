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
