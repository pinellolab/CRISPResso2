from CRISPResso2 import CRISPResso2Align, CRISPRessoShared
import re
import tempfile
import os
import gzip

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


def test_get_n_reads_fastq():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)

def test_get_n_reads_fastq_gzip():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # gzip file
    with open(f.name, 'rb') as f_in, gzip.open(f.name + '.gz', 'wb') as f_out:
        f_out.writelines(f_in)

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name + '.gz') == 1

    # Clean up: delete the file after the test
    os.remove(f.name)
    os.remove(f.name + '.gz')


def test_get_n_reads_fastq_three_extra_newlines():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.write("\n\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_four_extra_newlines():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.write("\n\n\n\n\n\n\n\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_100_reads():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        for i in range(100):
            f.write("@SEQ_ID\n")
            f.write("GATTACA\n")
            f.write("+\n")
            f.write("AAAAAAA\n")
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 100

    # Clean up: delete the file after the test
    os.remove(f.name)

def test_get_n_reads_fastq_no_newline():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_empty_file():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fastq') as f:
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 0

    # Clean up: delete the file after the test
    os.remove(f.name)


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


def test_get_quant_window_ranges_from_include_idxs():
    include_idxs = [0, 1, 2, 10, 11, 12]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(0, 2), (10, 12)]


def test_get_quant_window_ranges_from_include_idxs_empty():
    include_idxs = []
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == []


def test_get_quant_window_ranges_from_include_idxs_single():
    include_idxs = [50, 51, 52, 53]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53)]


def test_get_quant_window_ranges_from_include_idxs_single_gap():
    include_idxs = [50, 51, 52, 53, 55]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53), (55, 55)]


def test_get_quant_window_ranges_from_include_idxs_multiple_gaps():
    include_idxs = [50, 51, 52, 53, 55, 56, 57, 58, 60]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53), (55, 58), (60, 60)]


def test_get_silent_edits():
    ref_seq = 'AGS'
    seq = 'AGS'
    ref_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT')]
    seq_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGC')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'AGs'


def test_get_silent_edits_multiple_silent():
    ref_seq = 'PHATKIDS'
    seq = 'PHATKIDS'
    ref_codons = [('P', 'CCT'), ('H', 'CAT'), ('A', 'GCT'), ('T', 'ACT'), ('K', 'AAA'), ('I', 'ATT'), ('D', 'GAT'), ('S', 'AGC')]
    seq_codons = [('P', 'CCC'), ('H', 'CAC'), ('A', 'GCT'), ('T', 'ACT'), ('K', 'AAA'), ('I', 'ATT'), ('D', 'GAC'), ('S', 'AGT')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'phATKIds'


def test_get_silent_edits_no_silent():
    ref_seq = 'AGS'
    seq = 'AGT'
    ref_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT')]
    seq_codons = [('A', 'GCT'), ('G', 'GGT'), ('T', 'ACT')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'AGT'


def test_get_silent_edits_deletion():
    ref_seq = 'AGS'
    seq = 'AG-'
    ref_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT')]
    seq_codons = [('A', 'GCT'), ('G', 'GGT')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'AG-'


def test_get_silent_edits_multiple_deletions():
    ref_seq = 'AGS'
    seq = 'A--'
    ref_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT')]
    seq_codons = [('A', 'GCT')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'A--'


def test_get_silent_edits_insertion():
    ref_seq = 'AGS-'
    seq = 'AGST'
    ref_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT')]
    seq_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT'), ('T', 'ACT')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'AGST'

def test_get_silent_edits_middle_insertion():
    ref_seq = 'AG-S'
    seq = 'AGTS'
    ref_codons = [('A', 'GCT'), ('G', 'GGT'), ('S', 'AGT')]
    seq_codons = [('A', 'GCT'), ('G', 'GGT'),  ('T', 'ACT'), ('S', 'AGC')]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == 'AGTs'


# ── get_sgRNA_mismatch_vals tests ─────────────────────────────────────────────
# See PR #615: boundary condition when end_loc >= len(coords_r), which happens
# when a regex match extends to the very end of the from_sequence.

GAP_OPEN = -20
GAP_EXTEND = -2


def _get_mismatch_coords(seq1, seq2):
    """Get forward and reverse alignment coordinates for get_sgRNA_mismatch_vals."""
    coords_l, coords_r = CRISPRessoShared.get_alignment_coordinates(
        to_sequence=seq1, from_sequence=seq2,
        aln_matrix=ALN_MATRIX,
        needleman_wunsch_gap_open=GAP_OPEN,
        needleman_wunsch_gap_extend=GAP_EXTEND,
    )
    rev_coords_l, rev_coords_r = CRISPRessoShared.get_alignment_coordinates(
        to_sequence=seq2, from_sequence=seq1,
        aln_matrix=ALN_MATRIX,
        needleman_wunsch_gap_open=GAP_OPEN,
        needleman_wunsch_gap_extend=GAP_EXTEND,
    )
    return coords_l, coords_r, rev_coords_l, rev_coords_r


def test_sgRNA_mismatch_identical_no_mismatches():
    """Identical sequences should have no mismatches."""
    seq1 = 'ATCGTACG'
    seq2 = 'ATCGTACG'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 2, 5, coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_single_substitution():
    """A single substitution should be detected as a mismatch."""
    seq1 = 'ATCGTACG'
    seq2 = 'ATCATACG'  # position 3: G->A
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 2, 6, coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert sorted(mismatches) == [1]  # relative position 1 from start_loc=2


def test_sgRNA_mismatch_full_range_exclusive_end():
    """Full range with exclusive end (like match.end()) triggers the boundary."""
    seq1 = 'ATCGT'
    seq2 = 'ATCGT'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    # match.end() == len(seq2) == 5 — original bug: coords_r[5] was IndexError
    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 0, len(seq2), coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_interior_range():
    """Interior range of identical sequences should have no mismatches."""
    seq1 = 'ATCGTACGATCG'
    seq2 = 'ATCGTACGATCG'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 3, 8, coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_end_boundary_works():
    """end_loc == len(coords_r) should work, not IndexError."""
    seq1 = 'ATCGT'
    seq2 = 'ATCGT'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 0, len(coords_r), coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_last_position_detected_with_exclusive_end():
    """Mismatch at the last position should be detected when end_loc == len(coords_r)."""
    seq1 = 'ATCGA'  # last base A
    seq2 = 'ATCGT'  # last base T
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 0, len(coords_r), coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert sorted(mismatches) == [4]


def test_sgRNA_mismatch_last_position_missed_with_adjusted_end():
    """When end_loc = len(coords_r) - 1, range stops before the last position.

    This demonstrates that adjusting end_loc in the caller (as the PR does in
    get_prime_editing_guides) causes the last position to be skipped.
    """
    seq1 = 'ATCGA'
    seq2 = 'ATCGT'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    # Caller-adjusted end_loc: range(0, coords_r[4]) = range(0, 4)
    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 0, len(coords_r) - 1, coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []  # position 4 is not checked, so no mismatches found


def test_sgRNA_mismatch_longer_seq2_at_end():
    """When seq2 is longer and the match extends to its end.

    seq2='ATCGTAA' has extra bases that map back to seq1's last position.
    All bases in seq1[3:] match their mapped seq2 counterparts, so no mismatches.
    """
    seq1 = 'ATCGT'
    seq2 = 'ATCGTAA'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 3, len(coords_r), coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_deletion_in_seq1():
    """seq1 has a deletion relative to seq2 (C at position 2 is missing).

    seq1='ATGT' aligns to seq2='ATCGT' as AT-GT vs ATCGT.
    Position 2 in seq1 (G) maps to position 3 in seq2 (G) — a match, but
    the skipped seq2 position 2 (C) is detected as a deletion at relative
    position 2.
    """
    seq1 = 'ATGT'
    seq2 = 'ATCGT'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 0, 3, coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert sorted(mismatches) == [2]


def test_sgRNA_mismatch_insertion_at_end_boundary():
    """seq1 has an extra base (X) at the end beyond what seq2 covers.

    end_loc == len(coords_r) triggers the boundary, so end_loop = len(seq1) = 6.
    Position 5 in seq1 ('X') maps back to seq2 position 4 ('T') — a mismatch,
    and it's also a repeated rev_coords_l index (insertion), so relative
    position 5 is reported.
    """
    seq1 = 'ATCGTX'
    seq2 = 'ATCGT'
    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, 0, len(coords_r), coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert sorted(mismatches) == [5]


def test_sgRNA_mismatch_pe_extension_at_end():
    """Regression: PE extension at the very end of the amplicon.

    match.end() == len(amplicon), which triggered IndexError before the fix.
    """
    amplicon        = 'GCATGCATGCATCGTACG'
    edited_amplicon = 'GCATGCATGCATCGTACG'
    extension = 'CGTACG'

    match = re.search(extension, edited_amplicon)
    assert match.end() == len(edited_amplicon)

    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(
        amplicon, edited_amplicon,
    )
    assert match.end() >= len(coords_r)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        amplicon, edited_amplicon, match.start(), match.end(),
        coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_pe_extension_at_end_with_edit():
    """PE extension at end of amplicon where amplicons differ at last base."""
    amplicon   = 'GCATGCATGCATCGTACT'  # last base T
    edited_amp = 'GCATGCATGCATCGTACG'  # last base G
    extension = 'CGTACG'

    match = re.search(extension, edited_amp)
    assert match.end() == len(edited_amp)

    coords_l, coords_r = CRISPRessoShared.get_alignment_coordinates(
        to_sequence=amplicon, from_sequence=edited_amp,
        aln_matrix=ALN_MATRIX,
        needleman_wunsch_gap_open=GAP_OPEN,
        needleman_wunsch_gap_extend=GAP_EXTEND,
    )
    rev_coords_l, rev_coords_r = CRISPRessoShared.get_alignment_coordinates(
        to_sequence=edited_amp, from_sequence=amplicon,
        aln_matrix=ALN_MATRIX,
        needleman_wunsch_gap_open=GAP_OPEN,
        needleman_wunsch_gap_extend=GAP_EXTEND,
    )

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        amplicon, edited_amp, match.start(), match.end(),
        coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert sorted(mismatches) == [5]  # 6-base extension, mismatch at last base


def test_sgRNA_mismatch_pe_extension_not_at_end():
    """PE extension in the middle — normal case, no boundary issue."""
    amplicon = 'GCATGCATCGTACGAAAA'
    edited_amplicon = 'GCATGCATCGTACGAAAA'
    extension = 'CGTACG'

    match = re.search(extension, edited_amplicon)
    assert match.end() < len(edited_amplicon)

    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(
        amplicon, edited_amplicon,
    )
    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        amplicon, edited_amplicon, match.start(), match.end(),
        coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []


def test_sgRNA_mismatch_pe_regression_20bp():
    """Regression: 20bp amplicon with extension at the very end."""
    seq1 = 'AAAAATTTTTCCCCCGGGGG'
    seq2 = 'AAAAATTTTTCCCCCGGGGG'

    match = re.search('GGGGG', seq2)
    assert match.end() == 20

    coords_l, coords_r, rev_coords_l, rev_coords_r = _get_mismatch_coords(seq1, seq2)
    assert match.end() >= len(coords_r)

    mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(
        seq1, seq2, match.start(), match.end(),
        coords_l, coords_r, rev_coords_l, rev_coords_r,
    )
    assert mismatches == []
