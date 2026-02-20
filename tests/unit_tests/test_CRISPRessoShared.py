from CRISPResso2 import CRISPResso2Align, CRISPRessoShared
import tempfile
import os
import gzip

import numpy as np
import pandas as pd
import pytest

ALN_MATRIX = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")


# =============================================================================
# Tests for get_mismatches function
# =============================================================================


def test_get_mismatches():
    """Test get_mismatches function."""
    mismatch_cords = CRISPRessoShared.get_mismatches(
        "ATTA",
        "ATTA",
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 0

    mismatch_cords = CRISPRessoShared.get_mismatches(
        "GCAGTGGGCGCGCTA",
        "CCCACTGAAGGCCC",
        ALN_MATRIX,
        -5,
        -3,
    )
    expected_mismatch_count = 6
    assert len(mismatch_cords) == expected_mismatch_count


# =============================================================================
# Tests for get_relative_coordinates function
# =============================================================================


def test_get_relative_coordinates():
    """Test get_relative_coordinates function."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATCGT", "TTCGT")
    assert s1inds_gap_left == [0, 1, 2, 3, 4]
    assert s1inds_gap_right == [0, 1, 2, 3, 4]


# =============================================================================
# Tests for get_n_reads_fastq function
# =============================================================================


def test_get_n_reads_fastq():
    """Test get_n_reads_fastq function."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
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
    """Test get_n_reads_fastq function with gzipped file."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # gzip file
    with open(f.name, "rb") as f_in, gzip.open(f.name + ".gz", "wb") as f_out:
        f_out.writelines(f_in)

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name + ".gz") == 1

    # Clean up: delete the file after the test
    os.remove(f.name)
    os.remove(f.name + ".gz")


def test_get_n_reads_fastq_three_extra_newlines():
    """Test get_n_reads_fastq with three extra newlines."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
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
    """Test get_n_reads_fastq with four extra newlines."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
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
    """Test get_n_reads_fastq with 100 reads."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        for _i in range(100):
            f.write("@SEQ_ID\n")
            f.write("GATTACA\n")
            f.write("+\n")
            f.write("AAAAAAA\n")
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    expected_read_count = 100
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == expected_read_count

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_no_newline():
    """Test get_n_reads_fastq with no trailing newline."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
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
    """Test get_n_reads_fastq with empty file."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 0

    # Clean up: delete the file after the test
    os.remove(f.name)


# =============================================================================
# Tests for get_relative_coordinates function - gap handling
# =============================================================================


def test_get_relative_coordinates_to_gap():
    """Test get_relative_coordinates with gap in to_sequence."""
    # unaligned sequences
    seq_1 = "TTCGT"
    seq_2 = "TTCT"

    # aligned_sequences
    to_sequence = "TTC-T"
    from_sequence = "TTCGT"

    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates(to_sequence, from_sequence)
    assert s1inds_gap_left == [0, 1, 2, 2, 3]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]

    assert seq_1[0] == seq_2[s1inds_gap_left[0]]
    assert seq_1[1] == seq_2[s1inds_gap_left[1]]
    assert seq_1[2] == seq_2[s1inds_gap_left[2]]
    assert seq_1[4] == seq_2[s1inds_gap_left[4]]


def test_get_relative_coordinates_start_gap():
    """Test get_relative_coordinates with gap at start."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("--CGT", "TTCGT")
    assert s1inds_gap_left == [-1, -1, 0, 1, 2]
    assert s1inds_gap_right == [0, 0, 0, 1, 2]


def test_get_relative_coordinates_from_gap():
    """Test get_relative_coordinates with gap in from_sequence."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATCGT", "ATC-T")
    assert s1inds_gap_left == [0, 1, 2, 4]
    assert s1inds_gap_right == [0, 1, 2, 4]


def test_get_relative_coordinates_end_gap():
    """Test get_relative_coordinates with gap at end."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATC--", "ATCGT")
    assert s1inds_gap_left == [0, 1, 2, 2, 2]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]


def test_get_relative_coordinates_multiple_gaps():
    """Test get_relative_coordinates with multiple gaps."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("AT--TC--G--CC", "ATCGTCGCGTTCC")
    assert s1inds_gap_left == [0, 1, 1, 1, 2, 3, 3, 3, 4, 4, 4, 5, 6]
    assert s1inds_gap_right == [0, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 6]


def test_get_relative_coordinates_ind_and_dels():
    """Test get_relative_coordinates with insertions and deletions."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATG--C", "A-GCTC")
    assert s1inds_gap_left == [0, 2, 2, 2, 3]
    assert s1inds_gap_right == [0, 2, 3, 3, 3]


# =============================================================================
# Tests for get_quant_window_ranges_from_include_idxs function
# =============================================================================


def test_get_quant_window_ranges_from_include_idxs():
    """Test get_quant_window_ranges_from_include_idxs function."""
    include_idxs = [0, 1, 2, 10, 11, 12]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(0, 2), (10, 12)]


def test_get_quant_window_ranges_from_include_idxs_empty():
    """Test get_quant_window_ranges_from_include_idxs with empty list."""
    include_idxs = []
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == []


def test_get_quant_window_ranges_from_include_idxs_single():
    """Test get_quant_window_ranges_from_include_idxs with single range."""
    include_idxs = [50, 51, 52, 53]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53)]


def test_get_quant_window_ranges_from_include_idxs_single_gap():
    """Test get_quant_window_ranges_from_include_idxs with single gap."""
    include_idxs = [50, 51, 52, 53, 55]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53), (55, 55)]


def test_get_quant_window_ranges_from_include_idxs_multiple_gaps():
    """Test get_quant_window_ranges_from_include_idxs with multiple gaps."""
    include_idxs = [50, 51, 52, 53, 55, 56, 57, 58, 60]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53), (55, 58), (60, 60)]


# =============================================================================
# Tests for get_silent_edits function
# =============================================================================


def test_get_silent_edits():
    """Test get_silent_edits function."""
    ref_seq = "AGS"
    seq = "AGS"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGC")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGs"


def test_get_silent_edits_multiple_silent():
    """Test get_silent_edits with multiple silent edits."""
    ref_seq = "PHATKIDS"
    seq = "PHATKIDS"
    ref_codons = [("P", "CCT"), ("H", "CAT"), ("A", "GCT"), ("T", "ACT"), ("K", "AAA"), ("I", "ATT"), ("D", "GAT"), ("S", "AGC")]
    seq_codons = [("P", "CCC"), ("H", "CAC"), ("A", "GCT"), ("T", "ACT"), ("K", "AAA"), ("I", "ATT"), ("D", "GAC"), ("S", "AGT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "phATKIds"


def test_get_silent_edits_no_silent():
    """Test get_silent_edits with no silent edits."""
    ref_seq = "AGS"
    seq = "AGT"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("T", "ACT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGT"


def test_get_silent_edits_deletion():
    """Test get_silent_edits with deletion."""
    ref_seq = "AGS"
    seq = "AG-"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AG-"


def test_get_silent_edits_multiple_deletions():
    """Test get_silent_edits with multiple deletions."""
    ref_seq = "AGS"
    seq = "A--"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "A--"


def test_get_silent_edits_insertion():
    """Test get_silent_edits with insertion."""
    ref_seq = "AGS-"
    seq = "AGST"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT"), ("T", "ACT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGST"


def test_get_silent_edits_middle_insertion():
    """Test get_silent_edits with insertion in middle."""
    ref_seq = "AG-S"
    seq = "AGTS"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("T", "ACT"), ("S", "AGC")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGTs"


# =============================================================================
# Tests for reverse_complement function
# =============================================================================


def test_reverse_complement_basic():
    """Test reverse_complement with basic DNA sequence."""
    assert CRISPRessoShared.reverse_complement("ATCG") == "CGAT"


def test_reverse_complement_all_same():
    """Test reverse_complement with all same nucleotides."""
    assert CRISPRessoShared.reverse_complement("AAAA") == "TTTT"
    assert CRISPRessoShared.reverse_complement("TTTT") == "AAAA"
    assert CRISPRessoShared.reverse_complement("CCCC") == "GGGG"
    assert CRISPRessoShared.reverse_complement("GGGG") == "CCCC"


def test_reverse_complement_with_n():
    """Test reverse_complement with N (ambiguous) bases."""
    assert CRISPRessoShared.reverse_complement("ATNCG") == "CGNAT"
    assert CRISPRessoShared.reverse_complement("NNNN") == "NNNN"


def test_reverse_complement_single_base():
    """Test reverse_complement with single base."""
    assert CRISPRessoShared.reverse_complement("A") == "T"
    assert CRISPRessoShared.reverse_complement("T") == "A"
    assert CRISPRessoShared.reverse_complement("C") == "G"
    assert CRISPRessoShared.reverse_complement("G") == "C"


def test_reverse_complement_empty():
    """Test reverse_complement with empty string."""
    assert CRISPRessoShared.reverse_complement("") == ""


def test_reverse_complement_lowercase():
    """Test reverse_complement handles lowercase input."""
    assert CRISPRessoShared.reverse_complement("atcg") == "CGAT"


def test_reverse_complement_palindrome():
    """Test reverse_complement with palindromic sequence."""
    assert CRISPRessoShared.reverse_complement("ATAT") == "ATAT"
    assert CRISPRessoShared.reverse_complement("GCGC") == "GCGC"


# =============================================================================
# Tests for reverse function
# =============================================================================


def test_reverse_basic():
    """Test reverse with basic DNA sequence."""
    assert CRISPRessoShared.reverse("ATCG") == "GCTA"


def test_reverse_single():
    """Test reverse with single character."""
    assert CRISPRessoShared.reverse("A") == "A"


def test_reverse_empty():
    """Test reverse with empty string."""
    assert CRISPRessoShared.reverse("") == ""


def test_reverse_palindrome():
    """Test reverse with palindrome."""
    assert CRISPRessoShared.reverse("ATAT") == "TATA"


def test_reverse_lowercase():
    """Test reverse handles lowercase input."""
    assert CRISPRessoShared.reverse("atcg") == "GCTA"


# =============================================================================
# Tests for find_wrong_nt function
# =============================================================================


def test_find_wrong_nt_valid_sequence():
    """Test find_wrong_nt with valid DNA sequence."""
    assert CRISPRessoShared.find_wrong_nt("ATCG") == []
    assert CRISPRessoShared.find_wrong_nt("ATCGATCG") == []


def test_find_wrong_nt_with_n():
    """Test find_wrong_nt allows N as valid."""
    assert CRISPRessoShared.find_wrong_nt("ATCGN") == []
    assert CRISPRessoShared.find_wrong_nt("NNNN") == []


def test_find_wrong_nt_invalid_single():
    """Test find_wrong_nt with single invalid character."""
    result = CRISPRessoShared.find_wrong_nt("ATCGX")
    assert set(result) == {"X"}


def test_find_wrong_nt_invalid_multiple():
    """Test find_wrong_nt with multiple invalid characters."""
    result = CRISPRessoShared.find_wrong_nt("ATCGXYZ")
    assert set(result) == {"X", "Y", "Z"}


def test_find_wrong_nt_lowercase():
    """Test find_wrong_nt handles lowercase input."""
    assert CRISPRessoShared.find_wrong_nt("atcg") == []


def test_find_wrong_nt_empty():
    """Test find_wrong_nt with empty string."""
    assert CRISPRessoShared.find_wrong_nt("") == []


# =============================================================================
# Tests for capitalize_sequence function
# =============================================================================


def test_capitalize_sequence_lowercase():
    """Test capitalize_sequence with lowercase input."""
    assert CRISPRessoShared.capitalize_sequence("atcg") == "ATCG"


def test_capitalize_sequence_mixed():
    """Test capitalize_sequence with mixed case input."""
    assert CRISPRessoShared.capitalize_sequence("AtCg") == "ATCG"


def test_capitalize_sequence_already_upper():
    """Test capitalize_sequence with already uppercase input."""
    assert CRISPRessoShared.capitalize_sequence("ATCG") == "ATCG"


def test_capitalize_sequence_null():
    """Test capitalize_sequence with null value."""
    assert pd.isnull(CRISPRessoShared.capitalize_sequence(np.nan))


def test_capitalize_sequence_none():
    """Test capitalize_sequence with None value."""
    result = CRISPRessoShared.capitalize_sequence(None)
    # None is treated as null by pandas and returned as-is
    assert result is None


# =============================================================================
# Tests for slugify function
# =============================================================================


def test_slugify_spaces():
    """Test slugify replaces spaces with underscores."""
    assert CRISPRessoShared.slugify("hello world") == "hello_world"


def test_slugify_slashes():
    """Test slugify replaces slashes with underscores."""
    assert CRISPRessoShared.slugify("file/name") == "file_name"


def test_slugify_colons():
    """Test slugify replaces colons with underscores."""
    assert CRISPRessoShared.slugify("test:value") == "test_value"


def test_slugify_multiple_spaces():
    """Test slugify collapses multiple spaces/underscores."""
    assert CRISPRessoShared.slugify("multiple   spaces") == "multiple_spaces"


def test_slugify_special_chars():
    """Test slugify handles various special characters."""
    assert CRISPRessoShared.slugify("file*name") == "file_name"
    assert CRISPRessoShared.slugify("test[1]") == "test_1_"
    assert CRISPRessoShared.slugify("a<b>c") == "a_b_c"


def test_slugify_quotes():
    """Test slugify handles quotes."""
    assert CRISPRessoShared.slugify("it's") == "it_s"
    assert CRISPRessoShared.slugify('test"value') == "test_value"


def test_slugify_empty():
    """Test slugify with empty string."""
    assert CRISPRessoShared.slugify("") == ""


def test_slugify_already_clean():
    """Test slugify with already clean string."""
    assert CRISPRessoShared.slugify("clean_name") == "clean_name"


# =============================================================================
# Tests for get_amino_acids_and_codons function
# =============================================================================


def test_get_amino_acids_and_codons_basic():
    """Test get_amino_acids_and_codons with basic sequence."""
    result = CRISPRessoShared.get_amino_acids_and_codons("ATGGCT")
    assert result == [("M", "ATG"), ("A", "GCT")]


def test_get_amino_acids_and_codons_partial():
    """Test get_amino_acids_and_codons ignores incomplete codon at end."""
    result = CRISPRessoShared.get_amino_acids_and_codons("ATGGCTA")
    assert result == [("M", "ATG"), ("A", "GCT")]


def test_get_amino_acids_and_codons_stop():
    """Test get_amino_acids_and_codons with stop codon."""
    result = CRISPRessoShared.get_amino_acids_and_codons("ATGTAA")
    assert result == [("M", "ATG"), ("*", "TAA")]


def test_get_amino_acids_and_codons_single_codon():
    """Test get_amino_acids_and_codons with single codon."""
    result = CRISPRessoShared.get_amino_acids_and_codons("ATG")
    assert result == [("M", "ATG")]


def test_get_amino_acids_and_codons_too_short():
    """Test get_amino_acids_and_codons with sequence too short for codon."""
    result = CRISPRessoShared.get_amino_acids_and_codons("AT")
    assert result == []


def test_get_amino_acids_and_codons_empty():
    """Test get_amino_acids_and_codons with empty sequence."""
    result = CRISPRessoShared.get_amino_acids_and_codons("")
    assert result == []


# =============================================================================
# Tests for get_amino_acids_from_nucs function
# =============================================================================


def test_get_amino_acids_from_nucs_basic():
    """Test get_amino_acids_from_nucs with basic sequence."""
    assert CRISPRessoShared.get_amino_acids_from_nucs("ATGGCT") == "MA"


def test_get_amino_acids_from_nucs_with_deletion():
    """Test get_amino_acids_from_nucs with deletion (gap codon)."""
    assert CRISPRessoShared.get_amino_acids_from_nucs("ATG---GCT") == "M-A"


def test_get_amino_acids_from_nucs_stop():
    """Test get_amino_acids_from_nucs with stop codon."""
    assert CRISPRessoShared.get_amino_acids_from_nucs("ATGTAA") == "M*"


def test_get_amino_acids_from_nucs_partial():
    """Test get_amino_acids_from_nucs ignores incomplete codon."""
    assert CRISPRessoShared.get_amino_acids_from_nucs("ATGGCTA") == "MA"


def test_get_amino_acids_from_nucs_empty():
    """Test get_amino_acids_from_nucs with empty sequence."""
    assert CRISPRessoShared.get_amino_acids_from_nucs("") == ""


def test_get_amino_acids_from_nucs_multiple_deletions():
    """Test get_amino_acids_from_nucs with multiple deletions."""
    assert CRISPRessoShared.get_amino_acids_from_nucs("------") == "--"


# =============================================================================
# Tests for unexplode_cigar function
# =============================================================================


def test_unexplode_cigar_all_match():
    """Test unexplode_cigar with all matches."""
    assert CRISPRessoShared.unexplode_cigar("MMMMM") == ["5M"]


def test_unexplode_cigar_insertion_match():
    """Test unexplode_cigar with insertion and matches."""
    assert CRISPRessoShared.unexplode_cigar("IIMMMMMMM") == ["2I", "7M"]


def test_unexplode_cigar_match_deletion_match():
    """Test unexplode_cigar with deletion in middle."""
    assert CRISPRessoShared.unexplode_cigar("MMDDDMM") == ["2M", "3D", "2M"]


def test_unexplode_cigar_complex():
    """Test unexplode_cigar with complex pattern."""
    assert CRISPRessoShared.unexplode_cigar("MMMIIDDDDMMMM") == ["3M", "2I", "4D", "4M"]


def test_unexplode_cigar_single():
    """Test unexplode_cigar with single character."""
    assert CRISPRessoShared.unexplode_cigar("M") == ["1M"]
    assert CRISPRessoShared.unexplode_cigar("I") == ["1I"]
    assert CRISPRessoShared.unexplode_cigar("D") == ["1D"]


def test_unexplode_cigar_empty():
    """Test unexplode_cigar with empty string."""
    assert CRISPRessoShared.unexplode_cigar("") == []


# =============================================================================
# Tests for get_ref_length_from_cigar function
# =============================================================================


def test_get_ref_length_from_cigar_match_only():
    """Test get_ref_length_from_cigar with matches only."""
    assert CRISPRessoShared.get_ref_length_from_cigar("10M") == 10


def test_get_ref_length_from_cigar_with_insertion():
    """Test get_ref_length_from_cigar - insertions don't consume ref."""
    assert CRISPRessoShared.get_ref_length_from_cigar("5M3I5M") == 10


def test_get_ref_length_from_cigar_with_deletion():
    """Test get_ref_length_from_cigar - deletions consume ref."""
    assert CRISPRessoShared.get_ref_length_from_cigar("5M3D5M") == 13


def test_get_ref_length_from_cigar_complex():
    """Test get_ref_length_from_cigar with complex CIGAR."""
    assert CRISPRessoShared.get_ref_length_from_cigar("2M1I3M2D4M") == 11


def test_get_ref_length_from_cigar_single():
    """Test get_ref_length_from_cigar with single operation."""
    assert CRISPRessoShared.get_ref_length_from_cigar("1M") == 1
    assert CRISPRessoShared.get_ref_length_from_cigar("1D") == 1
    assert CRISPRessoShared.get_ref_length_from_cigar("1I") == 0


def test_get_ref_length_from_cigar_large_numbers():
    """Test get_ref_length_from_cigar with large numbers."""
    assert CRISPRessoShared.get_ref_length_from_cigar("100M50D25M") == 175


# =============================================================================
# Tests for clean_filename function
# =============================================================================


def test_clean_filename_spaces():
    """Test clean_filename replaces spaces."""
    assert CRISPRessoShared.clean_filename("hello world") == "hello_world"


def test_clean_filename_slashes():
    """Test clean_filename handles slashes."""
    assert CRISPRessoShared.clean_filename("test/file") == "test_file"


def test_clean_filename_special_chars():
    """Test clean_filename removes special characters."""
    result = CRISPRessoShared.clean_filename("file<>name")
    assert result == "file_name"


def test_clean_filename_colons():
    """Test clean_filename handles colons."""
    result = CRISPRessoShared.clean_filename("file:name:test")
    assert result == "file_name_test"


def test_clean_filename_already_clean():
    """Test clean_filename with already clean name."""
    assert CRISPRessoShared.clean_filename("clean_name") == "clean_name"


def test_clean_filename_numbers():
    """Test clean_filename preserves numbers."""
    assert CRISPRessoShared.clean_filename("test123") == "test123"


def test_clean_filename_empty():
    """Test clean_filename with empty string."""
    assert CRISPRessoShared.clean_filename("") == ""


def test_clean_filename_dots():
    """Test clean_filename preserves dots."""
    assert "." in CRISPRessoShared.clean_filename("file.txt")


def test_clean_filename_hyphens():
    """Test clean_filename preserves hyphens."""
    assert "-" in CRISPRessoShared.clean_filename("file-name")


# =============================================================================
# Tests for get_crispresso_logo function
# =============================================================================


def test_get_crispresso_logo():
    """Test get_crispresso_logo returns the expected ASCII art."""
    logo = CRISPRessoShared.get_crispresso_logo()
    assert "C)|" in logo
    assert "\\___/" in logo


# =============================================================================
# Tests for get_crispresso_footer function
# =============================================================================


def test_get_crispresso_footer():
    """Test get_crispresso_footer returns footer with expected ASCII art."""
    footer = CRISPRessoShared.get_crispresso_footer()
    assert "\\___/" in footer


# =============================================================================
# Tests for get_crispresso_header function
# =============================================================================


def test_get_crispresso_header():
    """Test get_crispresso_header contains the provided header string."""
    header = CRISPRessoShared.get_crispresso_header("Desc", "MyHeader")
    assert "MyHeader" in header


# =============================================================================
# Tests for format_cl_text function
# =============================================================================


def test_format_cl_text_basic():
    """Test format_cl_text returns text unchanged when short."""
    result = CRISPRessoShared.format_cl_text("test text")
    assert result == "test text"


def test_format_cl_text_with_max_chars():
    """Test format_cl_text wraps long text with newlines."""
    long_text = "a" * 200
    result = CRISPRessoShared.format_cl_text(long_text, max_chars=50)
    assert "\n" in result


def test_format_cl_text_converts_spaces_to_tabs():
    """Test format_cl_text handles spaces to tab conversion."""
    text = "    indented"
    result = CRISPRessoShared.format_cl_text(text, spaces_to_tab=4)
    assert result == "    indented"


def test_format_cl_text_empty():
    """Test format_cl_text with empty string."""
    result = CRISPRessoShared.format_cl_text("")
    assert result == ""


# =============================================================================
# Tests for is_C2Pro_installed function
# =============================================================================


def test_is_C2Pro_installed():
    """Test is_C2Pro_installed returns False when C2Pro is not installed."""
    result = CRISPRessoShared.is_C2Pro_installed()
    assert result is False


# =============================================================================
# Tests for check_file function
# =============================================================================


def test_check_file_valid_file():
    """Test check_file with valid file does not raise."""
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(b"test content")
        temp_path = f.name

    try:
        # check_file returns None on success, raises on failure
        result = CRISPRessoShared.check_file(temp_path)
        assert result is None  # Function returns None on success
    finally:
        os.remove(temp_path)


def test_check_file_nonexistent():
    """Test check_file with non-existent file raises exception."""
    with pytest.raises(Exception):
        CRISPRessoShared.check_file("/nonexistent/path/to/file.txt")


def test_check_file_empty_path():
    """Test check_file with empty path raises exception."""
    with pytest.raises(Exception):
        CRISPRessoShared.check_file("")


# =============================================================================
# Tests for insert_indels function
# =============================================================================


def test_insert_indels_no_indels():
    """Test insert_indels with no indels."""
    seq = "MAS"
    codons = [("M", "ATG"), ("A", "GCT"), ("S", "AGT")]
    result = CRISPRessoShared.insert_indels(seq, codons)
    assert result == codons


def test_insert_indels_with_gap():
    """Test insert_indels with gap in sequence."""
    seq = "M-S"
    codons = [("M", "ATG"), ("S", "AGT")]
    result = CRISPRessoShared.insert_indels(seq, codons)
    assert len(result) == 3
    assert result[1][0] == "-"  # Gap


def test_insert_indels_multiple_gaps():
    """Test insert_indels with multiple gaps."""
    seq = "M--S"
    codons = [("M", "ATG"), ("S", "AGT")]
    result = CRISPRessoShared.insert_indels(seq, codons)
    assert len(result) == 4


# =============================================================================
# Tests for check_output_folder function
# Note: check_output_folder requires a valid CRISPResso2 output folder with
# CRISPResso2_info.json, so we test error handling for invalid folders
# =============================================================================


def test_check_output_folder_invalid_raises():
    """Test check_output_folder raises for invalid folder.""" 
    with tempfile.TemporaryDirectory() as tmpdir:
        # An empty folder should raise as it's not a valid CRISPResso2 output
        with pytest.raises(CRISPRessoShared.OutputFolderIncompleteException):
            CRISPRessoShared.check_output_folder(tmpdir)


# =============================================================================
# Tests for force_symlink function
# =============================================================================


def test_force_symlink_creates_symlink():
    """Test force_symlink creates a symbolic link."""
    with tempfile.TemporaryDirectory() as tmpdir:
        src = os.path.join(tmpdir, "source.txt")
        dst = os.path.join(tmpdir, "link.txt")

        # Create source file
        with open(src, "w") as f:
            f.write("test content")

        CRISPRessoShared.force_symlink(src, dst)
        assert os.path.islink(dst)
        assert os.path.exists(dst)


def test_force_symlink_overwrites_existing():
    """Test force_symlink overwrites existing symlink."""
    with tempfile.TemporaryDirectory() as tmpdir:
        src1 = os.path.join(tmpdir, "source1.txt")
        src2 = os.path.join(tmpdir, "source2.txt")
        dst = os.path.join(tmpdir, "link.txt")

        # Create source files
        with open(src1, "w") as f:
            f.write("content 1")
        with open(src2, "w") as f:
            f.write("content 2")

        # Create first symlink
        CRISPRessoShared.force_symlink(src1, dst)
        # Overwrite with second
        CRISPRessoShared.force_symlink(src2, dst)

        assert os.path.islink(dst)
        # Should point to src2 now
        assert os.readlink(dst) == src2


# =============================================================================
# Tests for get_command_output function
# Note: get_command_output returns a generator
# =============================================================================


def test_get_command_output_echo():
    """Test get_command_output with echo command."""
    result = CRISPRessoShared.get_command_output("echo hello")
    # Function returns a generator, iterate to get output
    output = list(result)
    assert any("hello" in line for line in output)


def test_get_command_output_returns_generator():
    """Test get_command_output returns generator."""
    import types
    result = CRISPRessoShared.get_command_output("true")
    assert isinstance(result, types.GeneratorType)


def test_get_command_output_can_iterate():
    """Test get_command_output output can be iterated."""
    result = CRISPRessoShared.get_command_output("echo test")
    output = list(result)
    assert output[0] == "test\n"


# =============================================================================
# Tests for set_guide_array function
# =============================================================================


def test_set_guide_array_single_value():
    """Test set_guide_array with single value for multiple guides."""
    vals = "5"
    guides = ["guide1", "guide2", "guide3"]
    result = CRISPRessoShared.set_guide_array(vals, guides, "test_property")
    assert result == [5, 5, 5]


def test_set_guide_array_multiple_values():
    """Test set_guide_array with multiple values."""
    vals = "1,2,3"
    guides = ["guide1", "guide2", "guide3"]
    result = CRISPRessoShared.set_guide_array(vals, guides, "test_property")
    assert result == [1, 2, 3]


def test_set_guide_array_empty_guides():
    """Test set_guide_array with empty guides."""
    vals = "5"
    guides = []
    result = CRISPRessoShared.set_guide_array(vals, guides, "test_property")
    assert result == []


# =============================================================================
# Tests for parse_count_file function
# =============================================================================


def test_parse_count_file_valid():
    """Test parse_count_file with valid file."""
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        f.write("Amplicon\tA\tT\tC\tG\n")
        f.write("Ref\t10\t20\t30\t40\n")
        f.write("Mod\t5\t15\t25\t35\n")
        temp_path = f.name

    try:
        ampSeq, lab_freqs = CRISPRessoShared.parse_count_file(temp_path)
        assert ampSeq == "ATCG"
        assert "Ref" in lab_freqs
        assert "Mod" in lab_freqs
        assert lab_freqs["Ref"] == ["10", "20", "30", "40"]
    finally:
        os.remove(temp_path)


def test_parse_count_file_nonexistent():
    """Test parse_count_file with nonexistent file."""
    ampSeq, lab_freqs = CRISPRessoShared.parse_count_file("/nonexistent/file.txt")
    assert ampSeq is None
    assert lab_freqs is None


# =============================================================================
# Tests for parse_alignment_file function
# =============================================================================


def test_parse_alignment_file_valid():
    """Test parse_alignment_file with valid file."""
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        f.write("Amplicon\tA\tT\tC\tG\n")
        f.write("Read1\t1\t0\t0\t0\n")
        f.write("Read2\t0\t1\t1\t0\n")
        temp_path = f.name

    try:
        ampSeq, lab_freqs = CRISPRessoShared.parse_alignment_file(temp_path)
        assert ampSeq == "ATCG"
        assert "Read1" in lab_freqs
        assert "Read2" in lab_freqs
    finally:
        os.remove(temp_path)


def test_parse_alignment_file_nonexistent():
    """Test parse_alignment_file with nonexistent file."""
    ampSeq, lab_freqs = CRISPRessoShared.parse_alignment_file("/nonexistent/file.txt")
    assert ampSeq is None
    assert lab_freqs is None


# =============================================================================
# Tests for assert_fastq_format function
# =============================================================================


def test_assert_fastq_format_valid():
    """Test assert_fastq_format with valid FASTQ file."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("IIIIIII\n")
        temp_path = f.name

    try:
        result = CRISPRessoShared.assert_fastq_format(temp_path)
        assert result is True
    finally:
        os.remove(temp_path)


def test_assert_fastq_format_invalid_no_at():
    """Test assert_fastq_format with missing @ symbol."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fastq') as f:
        f.write("SEQ_ID\n")  # Missing @
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("IIIIIII\n")
        temp_path = f.name

    try:
        with pytest.raises(CRISPRessoShared.InputFileFormatException):
            CRISPRessoShared.assert_fastq_format(temp_path)
    finally:
        os.remove(temp_path)


def test_assert_fastq_format_invalid_no_plus():
    """Test assert_fastq_format with missing + symbol."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fastq') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("missing_plus\n")  # Should be +
        f.write("IIIIIII\n")
        temp_path = f.name

    try:
        with pytest.raises(CRISPRessoShared.InputFileFormatException):
            CRISPRessoShared.assert_fastq_format(temp_path)
    finally:
        os.remove(temp_path)


def test_assert_fastq_format_gzipped():
    """Test assert_fastq_format with gzipped FASTQ file."""
    import gzip

    with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as f:
        temp_path = f.name

    with gzip.open(temp_path, 'wt') as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("IIIIIII\n")

    try:
        result = CRISPRessoShared.assert_fastq_format(temp_path)
        assert result is True
    finally:
        os.remove(temp_path)


# =============================================================================
# Tests for getCRISPRessoArgParser function
# =============================================================================


def test_getCRISPRessoArgParser_core():
    """Test getCRISPRessoArgParser creates parser for Core tool."""
    parser = CRISPRessoShared.getCRISPRessoArgParser("Core")
    assert parser is not None
    # Should have version action
    assert "--version" in [a.option_strings[0] for a in parser._actions if a.option_strings]


def test_getCRISPRessoArgParser_batch():
    """Test getCRISPRessoArgParser creates parser for Batch tool."""
    parser = CRISPRessoShared.getCRISPRessoArgParser("Batch")
    assert parser is not None


def test_getCRISPRessoArgParser_pooled():
    """Test getCRISPRessoArgParser creates parser for Pooled tool."""
    parser = CRISPRessoShared.getCRISPRessoArgParser("Pooled")
    assert parser is not None


# =============================================================================
# Tests for get_core_crispresso_options function
# =============================================================================


def test_get_core_crispresso_options_returns_set():
    """Test get_core_crispresso_options returns a set of options."""
    options = CRISPRessoShared.get_core_crispresso_options()
    assert isinstance(options, set)
    assert len(options) > 0


def test_get_core_crispresso_options_contains_common():
    """Test get_core_crispresso_options contains common options."""
    options = CRISPRessoShared.get_core_crispresso_options()
    # Should contain some common options
    assert "fastq_r1" in options or "amplicon_seq" in options


# =============================================================================
# Tests for get_crispresso_options_lookup function
# =============================================================================


def test_get_crispresso_options_lookup_core():
    """Test get_crispresso_options_lookup for Core tool."""
    lookup = CRISPRessoShared.get_crispresso_options_lookup("Core")
    assert isinstance(lookup, dict)


def test_get_crispresso_options_lookup_contains_abbreviations():
    """Test get_crispresso_options_lookup contains abbreviations."""
    lookup = CRISPRessoShared.get_crispresso_options_lookup("Core")
    # r1 should map to fastq_r1
    if "r1" in lookup:
        assert lookup["r1"] == "fastq_r1"


# =============================================================================
# Tests for propagate_crispresso_options function
# =============================================================================


def test_propagate_crispresso_options_basic():
    """Test propagate_crispresso_options adds options to command."""
    cmd = "CRISPResso"
    options = ["name", "output_folder"]
    params = {"name": "test_sample", "output_folder": "/tmp/output"}

    result = CRISPRessoShared.propagate_crispresso_options(cmd, options, params)

    assert "--name test_sample" in result
    assert "--output_folder" in result


def test_propagate_crispresso_options_with_none():
    """Test propagate_crispresso_options handles None values."""
    cmd = "CRISPResso"
    options = ["name", "output_folder"]
    params = {"name": "test", "output_folder": None}

    result = CRISPRessoShared.propagate_crispresso_options(cmd, options, params)

    assert "--name test" in result
    assert "--output_folder" not in result


def test_propagate_crispresso_options_with_bool():
    """Test propagate_crispresso_options handles boolean values."""
    cmd = "CRISPResso"
    options = ["debug"]
    params = {"debug": True}

    result = CRISPRessoShared.propagate_crispresso_options(cmd, options, params)

    assert "--debug" in result


def test_propagate_crispresso_options_false_bool():
    """Test propagate_crispresso_options handles False boolean."""
    cmd = "CRISPResso"
    options = ["debug"]
    params = {"debug": False}

    result = CRISPRessoShared.propagate_crispresso_options(cmd, options, params)

    assert "--debug" not in result


def test_propagate_crispresso_options_with_spaces():
    """Test propagate_crispresso_options handles values with spaces."""
    cmd = "CRISPResso"
    options = ["name"]
    params = {"name": "sample with spaces"}

    result = CRISPRessoShared.propagate_crispresso_options(cmd, options, params)

    assert '"sample with spaces"' in result


# =============================================================================
# Tests for check_if_failed_run function
# =============================================================================


def test_check_if_failed_run_missing_info_file():
    """Test check_if_failed_run with missing info file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # No files created, should report as failed
        failed, message = CRISPRessoShared.check_if_failed_run(tmpdir, lambda x: None)
        assert failed is True


def test_check_if_failed_run_missing_status_file():
    """Test check_if_failed_run with missing status file."""
    import json

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create info file but not status file
        info_path = os.path.join(tmpdir, "CRISPResso2_info.json")
        with open(info_path, 'w') as f:
            json.dump({"test": "data"}, f)

        failed, message = CRISPRessoShared.check_if_failed_run(tmpdir, lambda x: None)
        assert failed is True


def test_check_if_failed_run_complete():
    """Test check_if_failed_run with complete run."""
    import json

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create both files with complete status
        info_path = os.path.join(tmpdir, "CRISPResso2_info.json")
        status_path = os.path.join(tmpdir, "CRISPResso_status.json")

        with open(info_path, 'w') as f:
            json.dump({"test": "data"}, f)

        with open(status_path, 'w') as f:
            json.dump({"percent_complete": 100.0, "status": "complete", "message": ""}, f)

        failed, message = CRISPRessoShared.check_if_failed_run(tmpdir, lambda x: None)
        assert failed is False


# =============================================================================
# Tests for set_console_log_level function
# =============================================================================


def test_set_console_log_level_basic():
    """Test set_console_log_level doesn't raise."""
    import logging
    logger = logging.getLogger("test_logger")
    # Should not raise
    CRISPRessoShared.set_console_log_level(logger, 3, debug=False)


def test_set_console_log_level_with_debug():
    """Test set_console_log_level with debug enabled."""
    import logging
    logger = logging.getLogger("test_logger_debug")
    # Should not raise
    CRISPRessoShared.set_console_log_level(logger, 4, debug=True)


# =============================================================================
# Tests for get_alignment_coordinates function
# =============================================================================


def test_get_alignment_coordinates_identical():
    """Test get_alignment_coordinates with identical sequences."""
    from CRISPResso2 import CRISPResso2Align
    aln_matrix = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")

    to_seq = "ATCGATCG"
    from_seq = "ATCGATCG"

    result = CRISPRessoShared.get_alignment_coordinates(
        to_seq, from_seq, aln_matrix, -20, -2
    )

    assert result is not None
    assert len(result) == 2  # Returns tuple of two lists


def test_get_alignment_coordinates_with_gap():
    """Test get_alignment_coordinates with gap in sequence."""
    from CRISPResso2 import CRISPResso2Align
    aln_matrix = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")

    to_seq = "ATCGATCG"
    from_seq = "ATCATCG"  # One base shorter

    result = CRISPRessoShared.get_alignment_coordinates(
        to_seq, from_seq, aln_matrix, -20, -2
    )

    assert result is not None


# =============================================================================
# Tests for zip_results function
# =============================================================================


def test_zip_results_creates_zip():
    """Test zip_results creates a zip file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some test files
        for i in range(3):
            with open(os.path.join(tmpdir, f"test_{i}.txt"), 'w') as f:
                f.write(f"content {i}")

        CRISPRessoShared.zip_results(tmpdir)

        # Check zip was created
        zip_path = tmpdir + ".zip"
        assert os.path.exists(zip_path)

        # Cleanup
        if os.path.exists(zip_path):
            os.remove(zip_path)


# =============================================================================
# Tests for additional edge cases
# =============================================================================


def test_reverse_complement_long_sequence():
    """Test reverse_complement with longer sequence."""
    seq = "ATCGATCGATCGATCGATCG"
    result = CRISPRessoShared.reverse_complement(seq)
    assert len(result) == len(seq)
    # Reverse complement of ATCG... should end with ...CGAT
    assert result.endswith("CGAT")


def test_find_wrong_nt_numbers():
    """Test find_wrong_nt with numbers in sequence."""
    result = CRISPRessoShared.find_wrong_nt("ATCG123")
    assert set(result) == {"1", "2", "3"}


def test_slugify_unicode_decomposable():
    """Test slugify decomposes accented characters to their base letter via NFKD normalization."""
    # NFKD decomposes ë -> e + combining diaeresis, then ascii encode drops the combining char
    assert CRISPRessoShared.slugify("café") == "cafe"
    assert CRISPRessoShared.slugify("naïve") == "naive"
    assert CRISPRessoShared.slugify("résumé") == "resume"


def test_slugify_unicode_non_decomposable():
    """Test slugify drops non-decomposable unicode characters (CJK, emoji)."""
    # CJK and emoji have no ASCII decomposition, so they are dropped entirely
    assert CRISPRessoShared.slugify("test_日本語") == "test_"
    assert CRISPRessoShared.slugify("emoji_😀") == "emoji_"


def test_clean_filename_pipe():
    """Test clean_filename handles pipe character."""
    result = CRISPRessoShared.clean_filename("file|name")
    assert result == "file_name"


def test_unexplode_cigar_soft_clip():
    """Test unexplode_cigar with soft clipping."""
    result = CRISPRessoShared.unexplode_cigar("SSMMMMS")
    assert result == ["2S", "4M", "1S"]


def test_get_ref_length_from_cigar_with_soft_clip():
    """Test get_ref_length_from_cigar ignores soft clips."""
    # S operations don't consume reference
    result = CRISPRessoShared.get_ref_length_from_cigar("5S10M5S")
    assert result == 10


# =============================================================================
# Tests for get_relative_coordinates function
# =============================================================================


def test_get_relative_coordinates_identical():
    """Test get_relative_coordinates with identical sequences."""
    ref = "ATCG"
    aln = "ATCG"
    s1inds, s2inds = CRISPRessoShared.get_relative_coordinates(ref, aln)

    assert list(s1inds) == [0, 1, 2, 3]
    assert list(s2inds) == [0, 1, 2, 3]


def test_get_relative_coordinates_with_gap_in_ref():
    """Test get_relative_coordinates with gap in reference."""
    ref = "AT--CG"
    aln = "ATGGCG"
    s1inds, s2inds = CRISPRessoShared.get_relative_coordinates(ref, aln)

    # Length should match alignment length
    assert len(s1inds) == 6
    assert len(s2inds) == 6


def test_get_relative_coordinates_with_gap_in_to_seq():
    """Test get_relative_coordinates with gap in to_sequence."""
    # Sequences must be same length for get_relative_coordinates
    # Both sequences need to be same length (aligned)
    to_seq = "AT--CG"  # Reference with gaps (insertions in from_seq)
    from_seq = "ATGGCG"  # Aligned sequence
    # When there are gaps, the function tracks positions
    s1inds_left, s1inds_right = CRISPRessoShared.get_relative_coordinates(to_seq, from_seq)

    # The function returns arrays for mapping positions
    assert len(s1inds_left) == len(to_seq)
    assert len(s1inds_right) == len(to_seq)


# =============================================================================
# Tests for codon translation
# =============================================================================


def test_codon_to_aa_lookup():
    """Test CODON_TO_AMINO_ACID lookup table exists and has entries."""
    assert hasattr(CRISPRessoShared, 'CODON_TO_AMINO_ACID')
    assert len(CRISPRessoShared.CODON_TO_AMINO_ACID) > 0
    assert CRISPRessoShared.CODON_TO_AMINO_ACID.get('ATG') == 'Met'


def test_codon_to_aa_single_char_lookup():
    """Test CODON_TO_AMINO_ACID_SINGLE_CHAR lookup table."""
    assert hasattr(CRISPRessoShared, 'CODON_TO_AMINO_ACID_SINGLE_CHAR')
    assert CRISPRessoShared.CODON_TO_AMINO_ACID_SINGLE_CHAR.get('ATG') == 'M'
    assert CRISPRessoShared.CODON_TO_AMINO_ACID_SINGLE_CHAR.get('TAA') == '*'


# =============================================================================
# Tests for get_quant_window_ranges_from_include_idxs function
# =============================================================================


def test_get_quant_window_ranges_continuous():
    """Test get_quant_window_ranges with continuous indices."""
    include_idxs = [0, 1, 2, 3, 4]
    result = CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs)

    assert result == [(0, 4)]


def test_get_quant_window_ranges_two_ranges():
    """Test get_quant_window_ranges with two separate ranges."""
    include_idxs = [0, 1, 2, 10, 11, 12]
    result = CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs)

    assert len(result) == 2
    assert result[0] == (0, 2)
    assert result[1] == (10, 12)


def test_get_quant_window_ranges_single_idx():
    """Test get_quant_window_ranges with single index."""
    include_idxs = [5]
    result = CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs)

    assert result == [(5, 5)]


def test_get_quant_window_ranges_empty():
    """Test get_quant_window_ranges with empty list."""
    result = CRISPRessoShared.get_quant_window_ranges_from_include_idxs([])
    assert result == []


# =============================================================================
# Tests for additional utility functions
# =============================================================================


def test_nt_complement_lookup():
    """Test nt_complement dictionary exists and is correct."""
    assert hasattr(CRISPRessoShared, 'nt_complement')
    assert CRISPRessoShared.nt_complement['A'] == 'T'
    assert CRISPRessoShared.nt_complement['T'] == 'A'
    assert CRISPRessoShared.nt_complement['C'] == 'G'
    assert CRISPRessoShared.nt_complement['G'] == 'C'


def test_cigar_lookup_exists():
    """Test CIGAR_LOOKUP dictionary exists."""
    assert hasattr(CRISPRessoShared, 'CIGAR_LOOKUP')
    # Check some key lookups
    assert CRISPRessoShared.CIGAR_LOOKUP[('A', 'A')] == 'M'
    assert CRISPRessoShared.CIGAR_LOOKUP[('A', '-')] == 'I'
    assert CRISPRessoShared.CIGAR_LOOKUP[('-', 'A')] == 'D'


def test_unexplode_cigar_basic():
    """Test unexplode_cigar with basic expanded cigar."""
    result = CRISPRessoShared.unexplode_cigar("MMMIID")
    assert result == ["3M", "2I", "1D"]


def test_unexplode_cigar_single_ops():
    """Test unexplode_cigar with single operations."""
    result = CRISPRessoShared.unexplode_cigar("MID")
    assert result == ["1M", "1I", "1D"]


# =============================================================================
# Tests for CIGAR operation extraction
# =============================================================================


def test_get_ref_length_from_cigar_basic():
    """Test get_ref_length_from_cigar with basic CIGAR."""
    result = CRISPRessoShared.get_ref_length_from_cigar("10M")
    assert result == 10


def test_get_ref_length_from_cigar_with_insertion():
    """Test get_ref_length_from_cigar with insertion (doesn't consume reference)."""
    result = CRISPRessoShared.get_ref_length_from_cigar("5M2I5M")
    assert result == 10  # Insertions don't consume reference


def test_get_ref_length_from_cigar_with_deletion():
    """Test get_ref_length_from_cigar with deletion (consumes reference)."""
    result = CRISPRessoShared.get_ref_length_from_cigar("5M2D5M")
    assert result == 12  # Deletions consume reference


def test_get_ref_length_from_cigar_complex():
    """Test get_ref_length_from_cigar with complex CIGAR."""
    result = CRISPRessoShared.get_ref_length_from_cigar("10M5I3D5M")
    # M=10, I=0 (doesn't consume ref), D=3, M=5 = 18
    assert result == 18


# =============================================================================
# Tests for exception classes
# =============================================================================


def test_bad_parameter_exception():
    """Test BadParameterException can be raised."""
    with pytest.raises(CRISPRessoShared.BadParameterException):
        raise CRISPRessoShared.BadParameterException("test error")


def test_auto_exception():
    """Test AutoException can be raised."""
    with pytest.raises(CRISPRessoShared.AutoException):
        raise CRISPRessoShared.AutoException("test auto error")


def test_output_folder_incomplete_exception():
    """Test OutputFolderIncompleteException can be raised."""
    with pytest.raises(CRISPRessoShared.OutputFolderIncompleteException):
        raise CRISPRessoShared.OutputFolderIncompleteException("test folder error")


def test_input_file_format_exception():
    """Test InputFileFormatException can be raised."""
    with pytest.raises(CRISPRessoShared.InputFileFormatException):
        raise CRISPRessoShared.InputFileFormatException("test format error")


# =============================================================================
# Tests for additional string manipulation
# =============================================================================


def test_reverse_identical():
    """Test reverse function with simple sequence."""
    result = CRISPRessoShared.reverse("ATCG")
    assert result == "GCTA"


def test_reverse_with_n():
    """Test reverse function handles N."""
    result = CRISPRessoShared.reverse("ATNCG")
    assert result == "GCNTA"


def test_find_wrong_nt_all_valid():
    """Test find_wrong_nt with all valid nucleotides."""
    result = CRISPRessoShared.find_wrong_nt("ATCGNATCGN")
    assert result == []


def test_find_wrong_nt_mixed():
    """Test find_wrong_nt with mixed valid and invalid."""
    result = CRISPRessoShared.find_wrong_nt("ATCGXYZ")
    assert set(result) == {"X", "Y", "Z"}


# =============================================================================
# Tests for CIGAR utility functions
# =============================================================================


def test_cigar_unexplode_pattern():
    """Test CIGAR_UNEXPLODE_PATTERN regex exists."""
    assert hasattr(CRISPRessoShared, 'cigarUnexplodePattern')
    # Test it can match CIGAR operations
    import re
    result = CRISPRessoShared.cigarUnexplodePattern.findall("MMMIID")
    assert len(result) > 0


def test_unexplode_cigar_all_same():
    """Test unexplode_cigar with all same operations."""
    result = CRISPRessoShared.unexplode_cigar("MMMMM")
    assert result == ["5M"]


def test_unexplode_cigar_alternating():
    """Test unexplode_cigar with alternating operations."""
    result = CRISPRessoShared.unexplode_cigar("MDMDMD")
    assert len(result) == 6
    assert result == ["1M", "1D", "1M", "1D", "1M", "1D"]


# =============================================================================
# Tests for check_custom_config function
# =============================================================================


def test_check_custom_config_none():
    """Test check_custom_config with None returns defaults."""
    result = CRISPRessoShared.check_custom_config(None)
    assert "colors" in result
    assert "guardrails" in result


def test_check_custom_config_returns_default_colors():
    """Test check_custom_config returns default color configuration."""
    result = CRISPRessoShared.check_custom_config(None)
    colors = result["colors"]
    assert colors == {
        "Substitution": "#0000FF",
        "Insertion": "#008000",
        "Deletion": "#FF0000",
        "A": "#7FC97F",
        "T": "#BEAED4",
        "C": "#FDC086",
        "G": "#FFFF99",
        "N": "#C8C8C8",
        "-": "#1E1E1E",
        "amino_acid_scheme": "unique",
    }
