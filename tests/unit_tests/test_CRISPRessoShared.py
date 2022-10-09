from CRISPResso2 import CRISPResso2Align, CRISPRessoShared, CRISPRessoObjects

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


    def test_get_amplicon_info_for_guides():
        #                                                                                     AGCCTCACTTCCTTCCCTAACAA
        amplicon_seq = "TTCCTTCCCTAACAACCATTAGATATCACATACATAGACTGAAGAAGAACCTGTCTTCTCTAAAAGACAAAGCCTCACTTCCTTCCCTAACAACCATTAGATATCACATACATAGAGATTAGCAGTATATTTCTCTAAGATTATATATCATATATAGAGAGATTACC"
        guide = CRISPRessoObjects.Guide(
            seq='AGCCTCACTTCCTTCCCTAACAA',name='testname',orig_seq='AGCCTCACTTCCTTCCCTAACAA',mismatches=[],
            qw_center=-3,qw_size=2,qw_size_5prime=5,qw_size_3prime=0,
            show_cut_point=True,pw_size=10,pw_size_5prime=5,pw_size_3prime=15
        )

        (guides_in_amplicon, this_include_idxs, this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(amplicon_seq, [guide],
                        quantification_window_coordinates=None, exclude_bp_from_left=5, exclude_bp_from_right=5,
                        discard_guide_positions_overhanging_amplicon_edge=False, shrink_plot_window_to_amplicon_size=False, shrink_quantification_window_to_included_bases=False)

        assert(len(guides_in_amplicon) == 1)
        assert(guides_in_amplicon[0].get_cut_point() == 89)
        assert(guides_in_amplicon[0].get_interval() == (70,92))
        assert(guides_in_amplicon[0].get_orientation() == 'F')
        assert(guides_in_amplicon[0].get_quantification_idx() == [85, 86, 87, 88, 89])

        #                                                                                     AGCCTCACTTCCTTCCCTAACAA
        #               TTCCTTCCCTAACAACCATTAGATAT
        #                                                                                                                                                                 ATCATATATAGAGAGATTACC
        amplicon_seq = "TTCCTTCCCTAACAACCATTAGATATCACATACATAGACTGAAGAAGAACCTGTCTTCTCTAAAAGACAAAGCCTCACTTCCTTCCCTAACAACCATTAGATATCACATACATAGAGATTAGCAGTATATTTCTCTAAGATTATATATCATATATAGAGAGATTACC"
        guide2 = CRISPRessoObjects.Guide(
            seq=CRISPRessoShared.reverse_complement('TTCCTTCCCTAACAACCATTAGATAT'),name='testname_2',orig_seq='TTCCTTCCCTAACAACCATTAGATAT',mismatches=[],
            qw_center=-3,qw_size=2,qw_size_5prime=5,qw_size_3prime=0,
            show_cut_point=True,pw_size=10,pw_size_5prime=5,pw_size_3prime=15
        )
        guide3 = CRISPRessoObjects.Guide(
            seq='ATCATATATAGAGAGATTACC',name='testname_2',orig_seq='ATCATATATAGAGAGATTACC',mismatches=[],
            qw_center=-1,qw_size=2,qw_size_5prime=5,qw_size_3prime=4,
            show_cut_point=True,pw_size=10,pw_size_5prime=5,pw_size_3prime=15
        )


        (guides_in_amplicon, this_include_idxs, this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(amplicon_seq, [guide,guide2,guide3],
                        quantification_window_coordinates=None, exclude_bp_from_left=5, exclude_bp_from_right=5,
                        discard_guide_positions_overhanging_amplicon_edge=False, shrink_plot_window_to_amplicon_size=True, shrink_quantification_window_to_included_bases=True)

        assert(len(guides_in_amplicon) == 4)
        assert(guides_in_amplicon[3].get_plot_idx() == (161, 166))
        assert(list(this_include_idxs) == [5,   6,   7,  81,  82,  83,  84,  85,  86,  87,  88,  89, 161])
        assert(list(this_exclude_idxs) == [  0,   1,   2,   3,   4, 162, 163, 164, 165, 166])

