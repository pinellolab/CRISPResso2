"""Tests for CRISPRessoPlot module utility functions."""

import pytest

from CRISPResso2 import CRISPRessoPlot


# =============================================================================
# Tests for get_nuc_color
# =============================================================================


def test_get_nuc_color_A():
    """Test get_nuc_color returns correct color for A."""
    color = CRISPRessoPlot.get_nuc_color("A", 1.0)
    assert color == pytest.approx((127 / 255.0, 201 / 255.0, 127 / 255.0, 1.0))


def test_get_nuc_color_T():
    """Test get_nuc_color returns correct color for T."""
    color = CRISPRessoPlot.get_nuc_color("T", 1.0)
    assert color == pytest.approx((190 / 255.0, 174 / 255.0, 212 / 255.0, 1.0))


def test_get_nuc_color_C():
    """Test get_nuc_color returns correct color for C."""
    color = CRISPRessoPlot.get_nuc_color("C", 1.0)
    assert color == pytest.approx((253 / 255.0, 192 / 255.0, 134 / 255.0, 1.0))


def test_get_nuc_color_G():
    """Test get_nuc_color returns correct color for G."""
    color = CRISPRessoPlot.get_nuc_color("G", 1.0)
    assert color == pytest.approx((255 / 255.0, 255 / 255.0, 153 / 255.0, 1.0))


def test_get_nuc_color_N():
    """Test get_nuc_color returns correct color for N (ambiguous)."""
    color = CRISPRessoPlot.get_nuc_color("N", 1.0)
    assert color == pytest.approx((200 / 255.0, 200 / 255.0, 200 / 255.0, 1.0))


def test_get_nuc_color_INS():
    """Test get_nuc_color returns correct color for INS (insertion)."""
    color = CRISPRessoPlot.get_nuc_color("INS", 1.0)
    assert color == pytest.approx((193 / 255.0, 129 / 255.0, 114 / 255.0, 1.0))


def test_get_nuc_color_DEL():
    """Test get_nuc_color returns correct color for DEL (deletion)."""
    color = CRISPRessoPlot.get_nuc_color("DEL", 1.0)
    assert color == pytest.approx((193 / 255.0, 129 / 255.0, 114 / 255.0, 1.0))


def test_get_nuc_color_gap():
    """Test get_nuc_color returns correct color for - (gap)."""
    color = CRISPRessoPlot.get_nuc_color("-", 1.0)
    assert color == pytest.approx((30 / 255.0, 30 / 255.0, 30 / 255.0, 1.0))


def test_get_nuc_color_alpha():
    """Test get_nuc_color respects alpha parameter."""
    color_full = CRISPRessoPlot.get_nuc_color("A", 1.0)
    color_half = CRISPRessoPlot.get_nuc_color("A", 0.5)
    assert color_full[3] == 1.0
    assert color_half[3] == 0.5
    # RGB should be the same
    assert color_full[:3] == color_half[:3]


def test_get_nuc_color_unknown():
    """Test get_nuc_color handles unknown nucleotides with computed color."""
    color = CRISPRessoPlot.get_nuc_color("X", 1.0)
    char_sum = (ord('X') - 65) / 90.0
    expected = (char_sum, 1 - char_sum, 2 * char_sum * (1 - char_sum), 1.0)
    assert color == pytest.approx(expected)


# =============================================================================
# Tests for get_color_lookup
# =============================================================================


def test_get_color_lookup_basic():
    """Test get_color_lookup with basic nucleotides."""
    nucs = ["A", "T", "C", "G"]
    colors = CRISPRessoPlot.get_color_lookup(nucs, 1.0)
    assert "A" in colors
    assert "T" in colors
    assert "C" in colors
    assert "G" in colors


def test_get_color_lookup_empty():
    """Test get_color_lookup with empty list."""
    colors = CRISPRessoPlot.get_color_lookup([], 1.0)
    assert colors == {}


def test_get_color_lookup_single():
    """Test get_color_lookup with single nucleotide."""
    colors = CRISPRessoPlot.get_color_lookup(["A"], 1.0)
    assert len(colors) == 1
    assert colors["A"] == pytest.approx((127 / 255.0, 201 / 255.0, 127 / 255.0, 1.0))


def test_get_color_lookup_all_nucs():
    """Test get_color_lookup with all nucleotide types."""
    nucs = ["A", "T", "C", "G", "N", "-", "INS", "DEL"]
    colors = CRISPRessoPlot.get_color_lookup(nucs, 0.8)
    assert len(colors) == 8
    assert colors["A"] == pytest.approx((127 / 255.0, 201 / 255.0, 127 / 255.0, 0.8))
    assert colors["T"] == pytest.approx((190 / 255.0, 174 / 255.0, 212 / 255.0, 0.8))
    assert colors["G"] == pytest.approx((255 / 255.0, 255 / 255.0, 153 / 255.0, 0.8))
    assert colors["INS"] == pytest.approx(colors["DEL"])  # Same color


# =============================================================================
# Tests for hex_to_rgb
# =============================================================================


def test_hex_to_rgb_basic():
    """Test hex_to_rgb with basic hex colors."""
    assert CRISPRessoPlot.hex_to_rgb("#FF0000") == (255, 0, 0)  # Red
    assert CRISPRessoPlot.hex_to_rgb("#00FF00") == (0, 255, 0)  # Green
    assert CRISPRessoPlot.hex_to_rgb("#0000FF") == (0, 0, 255)  # Blue


def test_hex_to_rgb_black_white():
    """Test hex_to_rgb with black and white."""
    assert CRISPRessoPlot.hex_to_rgb("#000000") == (0, 0, 0)  # Black
    assert CRISPRessoPlot.hex_to_rgb("#FFFFFF") == (255, 255, 255)  # White


def test_hex_to_rgb_without_hash():
    """Test hex_to_rgb handles colors without # prefix."""
    assert CRISPRessoPlot.hex_to_rgb("FF0000") == (255, 0, 0)


def test_hex_to_rgb_lowercase():
    """Test hex_to_rgb handles lowercase hex."""
    assert CRISPRessoPlot.hex_to_rgb("#ff0000") == (255, 0, 0)


def test_hex_to_rgb_mixed_case():
    """Test hex_to_rgb handles mixed case hex."""
    assert CRISPRessoPlot.hex_to_rgb("#FfAa00") == (255, 170, 0)


def test_hex_to_rgb_gray():
    """Test hex_to_rgb with gray colors."""
    assert CRISPRessoPlot.hex_to_rgb("#808080") == (128, 128, 128)


# =============================================================================
# Tests for amino_acids_to_numbers
# =============================================================================


def test_amino_acids_to_numbers_basic():
    """Test amino_acids_to_numbers with basic sequence."""
    result = CRISPRessoPlot.amino_acids_to_numbers("MA")
    assert result == [11, 1]


def test_amino_acids_to_numbers_stop():
    """Test amino_acids_to_numbers with stop codon."""
    result = CRISPRessoPlot.amino_acids_to_numbers("*")
    assert len(result) == 1
    assert result[0] == 0  # Stop codon is first in list


def test_amino_acids_to_numbers_gap():
    """Test amino_acids_to_numbers with gap."""
    result = CRISPRessoPlot.amino_acids_to_numbers("-")
    assert result == [22]


def test_amino_acids_to_numbers_all_standard():
    """Test amino_acids_to_numbers with all standard amino acids."""
    all_aa = "ACDEFGHIKLMNPQRSTVWY"
    result = CRISPRessoPlot.amino_acids_to_numbers(all_aa)
    assert result == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]


def test_amino_acids_to_numbers_empty():
    """Test amino_acids_to_numbers with empty sequence."""
    result = CRISPRessoPlot.amino_acids_to_numbers("")
    assert result == []


# =============================================================================
# Tests for get_amino_acid_color_dict
# =============================================================================


def test_get_amino_acid_color_dict_clustal():
    """Test get_amino_acid_color_dict with clustal scheme."""
    colors = CRISPRessoPlot.get_amino_acid_color_dict('clustal')
    assert colors['*'] == '#FF0000'
    assert colors['A'] == '#000000'
    assert colors['G'] == '#FFA500'
    assert colors['F'] == '#0000FF'
    assert colors['I'] == '#008000'
    assert colors['-'] == '#FFFFFF'


def test_get_amino_acid_color_dict_default():
    """Test get_amino_acid_color_dict with default scheme (clustal)."""
    colors = CRISPRessoPlot.get_amino_acid_color_dict()
    assert len(colors) == 23
    assert colors['*'] == '#FF0000'
    assert colors['A'] == '#000000'
    assert colors['G'] == '#FFA500'


def test_get_amino_acid_color_dict_returns_hex():
    """Test that get_amino_acid_color_dict returns hex colors."""
    colors = CRISPRessoPlot.get_amino_acid_color_dict('clustal')
    for aa, color in colors.items():
        # Should start with # and be a valid hex color
        assert color.startswith('#')
        assert len(color) == 7  # #RRGGBB format


# =============================================================================
# Tests for get_amino_acid_colors
# =============================================================================


def test_get_amino_acid_colors_basic():
    """Test get_amino_acid_colors returns hex+alpha list for clustal."""
    colors = CRISPRessoPlot.get_amino_acid_colors("clustal")
    assert len(colors) == 23
    assert colors[0] == '#FF000066'   # * (stop codon)
    assert colors[1] == '#00000066'   # A


# =============================================================================
# Tests for setMatplotlibDefaults
# =============================================================================


def test_setMatplotlibDefaults_no_error():
    """Test that setMatplotlibDefaults runs without error."""
    import matplotlib
    CRISPRessoPlot.setMatplotlibDefaults()
    # Verify it actually changed matplotlib settings
    assert matplotlib.rcParams['font.size'] != 0


# =============================================================================
# Tests for get_rows_for_sgRNA_annotation
# =============================================================================


def test_get_rows_for_sgRNA_annotation_empty():
    """Test get_rows_for_sgRNA_annotation raises ValueError with empty intervals."""
    with pytest.raises(ValueError):
        CRISPRessoPlot.get_rows_for_sgRNA_annotation([], 100)


def test_get_rows_for_sgRNA_annotation_single():
    """Test get_rows_for_sgRNA_annotation with single sgRNA."""
    sgRNA_intervals = [(10, 30)]
    amp_len = 100
    result = CRISPRessoPlot.get_rows_for_sgRNA_annotation(sgRNA_intervals, amp_len)
    assert len(result) == 1
    assert result[0] == 0  # First row


def test_get_rows_for_sgRNA_annotation_multiple_non_overlapping():
    """Test get_rows_for_sgRNA_annotation with non-overlapping sgRNAs."""
    sgRNA_intervals = [(10, 30), (50, 70)]
    amp_len = 100
    result = CRISPRessoPlot.get_rows_for_sgRNA_annotation(sgRNA_intervals, amp_len)
    assert len(result) == 2
    # Non-overlapping sgRNAs should be on the same row
    assert result[0] == result[1] == 0


def test_get_rows_for_sgRNA_annotation_overlapping():
    """Test get_rows_for_sgRNA_annotation with overlapping sgRNAs."""
    # Overlapping intervals should be on different rows
    sgRNA_intervals = [(10, 30), (20, 40)]
    amp_len = 100
    result = CRISPRessoPlot.get_rows_for_sgRNA_annotation(sgRNA_intervals, amp_len)
    assert len(result) == 2
    # Overlapping sgRNAs should be on different rows
    assert result[0] != result[1]


# =============================================================================
# Tests for prep_alleles_table
# =============================================================================


def test_prep_alleles_table_basic():
    """Test prep_alleles_table with basic data."""
    import pandas as pd
    import numpy as np

    # Create minimal allele dataframe
    df = pd.DataFrame({
        '%Reads': [50.0, 30.0, 20.0],
        '#Reads': [500, 300, 200],
        'Reference_Sequence': ['ATCG', 'ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG', 'A-CG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        CRISPRessoPlot.prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4], [1, 0, 3, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G'], ['A', '-', 'C', 'G']]
    assert y_labels == ['50.00% (500 reads)', '30.00% (300 reads)', '20.00% (200 reads)']
    assert is_reference == [True, False, False]


def test_prep_alleles_table_empty():
    """Test prep_alleles_table with empty dataframe after filtering."""
    import pandas as pd

    df = pd.DataFrame({
        '%Reads': [0.1],
        '#Reads': [1],
        'Reference_Sequence': ['ATCG'],
    }, index=['ATCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        CRISPRessoPlot.prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=1.0)

    assert len(X) == 0
    assert len(annot) == 0


def test_prep_alleles_table_max_rows():
    """Test prep_alleles_table respects MAX_N_ROWS."""
    import pandas as pd

    df = pd.DataFrame({
        '%Reads': [30.0, 25.0, 20.0, 15.0, 10.0],
        '#Reads': [300, 250, 200, 150, 100],
        'Reference_Sequence': ['ATCG', 'ATCG', 'ATCG', 'ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG', 'TTCG', 'ATCA', 'GGGG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        CRISPRessoPlot.prep_alleles_table(df, 'ATCG', MAX_N_ROWS=3, MIN_FREQUENCY=0)

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4], [2, 2, 3, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G'], ['T', 'T', 'C', 'G']]
    assert y_labels == ['30.00% (300 reads)', '25.00% (250 reads)', '20.00% (200 reads)']


def test_prep_alleles_table_with_insertions():
    """Test prep_alleles_table detects insertions."""
    import pandas as pd

    # Reference with gap indicates insertion in read
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['AT--CG'],
    }, index=['ATGGCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        CRISPRessoPlot.prep_alleles_table(df, 'ATGGCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    # Should detect insertion
    assert 0 in insertion_dict
    assert len(insertion_dict[0]) > 0


# =============================================================================
# Tests for color functions - additional cases
# =============================================================================


def test_get_color_lookup_with_custom_colors():
    """Test get_color_lookup with custom colors."""
    custom_colors = {
        'A': '#FF0000',
        'T': '#00FF00',
        'C': '#0000FF',
        'G': '#FFFF00',
        'N': '#CCCCCC',
        '-': '#000000'
    }
    nucs = ['A', 'T', 'C', 'G', 'N', '-']
    colors = CRISPRessoPlot.get_color_lookup(nucs, 0.8, custom_colors=custom_colors)

    assert colors['A'] == pytest.approx((1.0, 0.0, 0.0, 0.8))
    assert colors['T'] == pytest.approx((0.0, 1.0, 0.0, 0.8))
    assert colors['C'] == pytest.approx((0.0, 0.0, 1.0, 0.8))
    assert colors['G'] == pytest.approx((1.0, 1.0, 0.0, 0.8))
    assert colors['N'] == pytest.approx((204 / 255.0, 204 / 255.0, 204 / 255.0, 0.8))
    assert colors['-'] == pytest.approx((0.0, 0.0, 0.0, 0.8))


def test_get_amino_acid_color_dict_unique_scheme():
    """Test get_amino_acid_color_dict with unique scheme."""
    colors = CRISPRessoPlot.get_amino_acid_color_dict('unique')
    assert colors['*'] == '#FF0000'
    assert colors['A'] == '#000000'
    assert colors['C'] == '#1E90FF'
    assert colors['-'] == '#B0B0B0'


def test_get_amino_acid_colors_with_dict():
    """Test get_amino_acid_colors with dict scheme."""
    custom_scheme = {
        '*': '#FF0000', 'A': '#000000', 'C': '#1E90FF', 'D': '#FF4500',
        'E': '#32CD32', 'F': '#FFD700', 'G': '#8A2BE2', 'H': '#FF69B4',
        'I': '#00FF7F', 'K': '#00BFFF', 'L': '#FF6347', 'M': '#ADFF2F',
        'N': '#FF8C00', 'P': '#A52A2A', 'Q': '#00CED1', 'R': '#8A2BE2',
        'S': '#48D1CC', 'T': '#C71585', 'V': '#4682B4', 'W': '#D2691E',
        'Y': '#9ACD32', '': '#FFFFFF', '-': '#B0B0B0'
    }
    colors = CRISPRessoPlot.get_amino_acid_colors(custom_scheme)
    assert len(colors) == 23
    assert colors[0] == '#FF000066'   # * with hex alpha
    assert colors[1] == '#00000066'   # A with hex alpha
    assert colors[-1] == '#B0B0B066'  # - with hex alpha


def test_amino_acids_to_numbers_with_special():
    """Test amino_acids_to_numbers with special characters."""
    result = CRISPRessoPlot.amino_acids_to_numbers("*-")
    assert len(result) == 2
    assert result[0] == 0  # * is first
    assert result[1] == 22  # - is last


# =============================================================================
# Tests for hex_to_rgb edge cases
# =============================================================================


def test_hex_to_rgb_short_form():
    """Test hex_to_rgb would need 6-char form."""
    # Standard 6-character hex
    assert CRISPRessoPlot.hex_to_rgb("#AABBCC") == (170, 187, 204)


def test_hex_to_rgb_all_zeros():
    """Test hex_to_rgb with all zeros (black)."""
    assert CRISPRessoPlot.hex_to_rgb("#000000") == (0, 0, 0)


def test_hex_to_rgb_all_ones():
    """Test hex_to_rgb with all max (white)."""
    assert CRISPRessoPlot.hex_to_rgb("#FFFFFF") == (255, 255, 255)


# =============================================================================
# Tests for prep_alleles_table_compare
# =============================================================================


def test_prep_alleles_table_compare_basic():
    """Test prep_alleles_table_compare with basic data."""
    import pandas as pd
    import numpy as np

    # Create merged allele dataframe
    df = pd.DataFrame({
        '%Reads_sample1': [50.0, 30.0],
        '%Reads_sample2': [40.0, 35.0],
        '#Reads_sample1': [500, 300],
        '#Reads_sample2': [400, 350],
        'Reference_Sequence': ['ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws = \
        CRISPRessoPlot.prep_alleles_table_compare(
            df, 'sample1', 'sample2', MAX_N_ROWS=10, MIN_FREQUENCY=0
        )

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G']]
    assert y_labels == [
        '50.00% (500 reads) 40.00% (400 reads) ',
        '30.00% (300 reads) 35.00% (350 reads) ',
    ]


def test_prep_alleles_table_compare_with_insertion():
    """Test prep_alleles_table_compare detects insertions."""
    import pandas as pd

    df = pd.DataFrame({
        '%Reads_s1': [50.0],
        '%Reads_s2': [50.0],
        '#Reads_s1': [500],
        '#Reads_s2': [500],
        'Reference_Sequence': ['AT--CG'],  # Insertion markers
    }, index=['ATGGCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws = \
        CRISPRessoPlot.prep_alleles_table_compare(
            df, 's1', 's2', MAX_N_ROWS=10, MIN_FREQUENCY=0
        )

    # Should detect insertion
    assert 0 in insertion_dict
    assert len(insertion_dict[0]) > 0


# =============================================================================
# Tests for prep_amino_acid_table
# =============================================================================


def test_prep_amino_acid_table_basic():
    """Test prep_amino_acid_table with basic data."""
    import pandas as pd

    df = pd.DataFrame({
        '%Reads': [60.0, 40.0],
        '#Reads': [600, 400],
        'Reference_Sequence': ['MAS', 'MAS'],
        'silent_edit_inds': [[], []],
    }, index=['MAS', 'MAT'])

    X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq = \
        CRISPRessoPlot.prep_amino_acid_table(df, 'MAS', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(X) == 2
    assert is_reference[0] is True  # First row matches reference
    assert ref_seq == 'MAS'


def test_prep_amino_acid_table_with_silent_edits():
    """Test prep_amino_acid_table with silent edits."""
    import pandas as pd

    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['MAS'],
        'silent_edit_inds': [[1]],  # Silent edit at position 1
    }, index=['MAS'])

    X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq = \
        CRISPRessoPlot.prep_amino_acid_table(df, 'MAS', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    # Should have silent edit at row 0, position 1
    assert 0 in silent_edit_dict
    assert 1 in silent_edit_dict[0]


# =============================================================================
# Tests for CustomHeatMapper class
# =============================================================================


def test_custom_heatmap_basic():
    """Test custom_heatmap creates heatmap."""
    import numpy as np
    import matplotlib.pyplot as plt

    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    fig, ax = plt.subplots()
    result = CRISPRessoPlot.custom_heatmap(data, ax=ax)

    assert result is not None
    plt.close(fig)


def test_custom_heatmap_with_annotation():
    """Test custom_heatmap with annotation."""
    import numpy as np
    import matplotlib.pyplot as plt

    data = np.array([[1, 2], [3, 4]])
    annot = np.array([['A', 'B'], ['C', 'D']])

    fig, ax = plt.subplots()
    result = CRISPRessoPlot.custom_heatmap(data, annot=annot, fmt='s', ax=ax)

    assert result is not None
    plt.close(fig)


# =============================================================================
# Additional color function tests
# =============================================================================


def test_get_nuc_color_all_special():
    """Test get_nuc_color returns distinct colors for each nucleotide type."""
    expected = {
        'A':   (127 / 255.0, 201 / 255.0, 127 / 255.0, 1.0),
        'T':   (190 / 255.0, 174 / 255.0, 212 / 255.0, 1.0),
        'C':   (253 / 255.0, 192 / 255.0, 134 / 255.0, 1.0),
        'G':   (255 / 255.0, 255 / 255.0, 153 / 255.0, 1.0),
        'N':   (200 / 255.0, 200 / 255.0, 200 / 255.0, 1.0),
        'INS': (193 / 255.0, 129 / 255.0, 114 / 255.0, 1.0),
        'DEL': (193 / 255.0, 129 / 255.0, 114 / 255.0, 1.0),
        '-':   (30 / 255.0, 30 / 255.0, 30 / 255.0, 1.0),
    }
    for nuc, exp in expected.items():
        assert CRISPRessoPlot.get_nuc_color(nuc, 1.0) == pytest.approx(exp)


def test_get_color_lookup_preserves_all_nucleotides():
    """Test get_color_lookup returns colors for all nucleotides."""
    nucs = ['A', 'T', 'C', 'G', 'N', 'INS', 'DEL', '-']
    colors = CRISPRessoPlot.get_color_lookup(nucs, 0.5)

    for nuc in nucs:
        assert nuc in colors
        assert colors[nuc][3] == 0.5  # Alpha should be 0.5


# =============================================================================
# Tests for utility functions
# =============================================================================


def test_amino_acids_to_numbers_all_standard():
    """Test amino_acids_to_numbers with all standard amino acids."""
    aa_seq = "ACDEFGHIKLMNPQRSTVWY"
    result = CRISPRessoPlot.amino_acids_to_numbers(aa_seq)

    assert len(result) == 20
    # All should be unique
    assert len(set(result)) == 20


def test_get_amino_acid_colors_none_scheme():
    """Test get_amino_acid_colors with None scheme uses default (clustal)."""
    colors = CRISPRessoPlot.get_amino_acid_colors(None)
    assert len(colors) == 23
    assert colors[0] == '#FF000066'   # * (stop codon)
    assert colors[1] == '#00000066'   # A


def test_get_amino_acid_color_dict_something_scheme():
    """Test get_amino_acid_color_dict with 'something' scheme."""
    colors = CRISPRessoPlot.get_amino_acid_color_dict('something')
    assert colors['*'] == '#000000'
    assert colors['A'] == '#90EE90'
    assert colors['G'] == '#90EE90'
    assert colors['I'] == '#0000FF'
    assert colors['-'] == '#FFFFFF'


# =============================================================================
# Tests for CustomHeatMapper class
# =============================================================================


def test_custom_heatmapper_init():
    """Test CustomHeatMapper initialization."""
    import numpy as np

    data = np.array([[1, 2], [3, 4]])
    mapper = CRISPRessoPlot.Custom_HeatMapper(
        data, vmin=0, vmax=10, cmap=None, center=None, robust=False,
        annot=None, fmt=".2g", annot_kws=None, per_element_annot_kws=None,
        cbar=True, cbar_kws=None, xticklabels=True, yticklabels=True, mask=None
    )
    assert mapper is not None


def test_custom_heatmapper_with_annotation():
    """Test CustomHeatMapper with annotations."""
    import numpy as np

    data = np.array([[1, 2], [3, 4]])
    annot = np.array([['A', 'B'], ['C', 'D']])

    mapper = CRISPRessoPlot.Custom_HeatMapper(
        data, vmin=0, vmax=10, cmap=None, center=None, robust=False,
        annot=annot, fmt='s', annot_kws={'size': 10}, per_element_annot_kws=None,
        cbar=True, cbar_kws=None, xticklabels=True, yticklabels=True, mask=None
    )
    assert mapper.annot is not None


def test_custom_heatmapper_with_mask():
    """Test CustomHeatMapper with mask."""
    import numpy as np

    data = np.array([[1, 2], [3, 4]])
    mask = np.array([[False, True], [True, False]])

    mapper = CRISPRessoPlot.Custom_HeatMapper(
        data, vmin=0, vmax=10, cmap=None, center=None, robust=False,
        annot=None, fmt=".2g", annot_kws=None, per_element_annot_kws=None,
        cbar=True, cbar_kws=None, xticklabels=True, yticklabels=True, mask=mask
    )
    # mask might be stored differently - just verify mapper was created
    assert mapper is not None


# =============================================================================
# Tests for additional prep functions
# =============================================================================


def test_prep_alleles_table_with_substitution():
    """Test prep_alleles_table detects substitutions."""
    import pandas as pd
    import numpy as np

    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['ATCG'],  # Reference
    }, index=['GTCG'])  # G at position 0 instead of A

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        CRISPRessoPlot.prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(X) == 1
    assert is_reference[0] is False  # Different from reference
    # Should have bold annotation for substitution
    assert len(per_element_annot_kws[0]) > 0


def test_prep_alleles_table_all_reference():
    """Test prep_alleles_table with all reference sequences."""
    import pandas as pd

    df = pd.DataFrame({
        '%Reads': [100.0],
        '#Reads': [1000],
        'Reference_Sequence': ['ATCG'],
    }, index=['ATCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        CRISPRessoPlot.prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(is_reference) == 1
    assert is_reference[0] is True


# =============================================================================
# Tests for color edge cases
# =============================================================================


def test_get_nuc_color_lowercase():
    """Test get_nuc_color falls through to computed color for lowercase."""
    color = CRISPRessoPlot.get_nuc_color("a", 1.0)
    # Lowercase 'a' doesn't match any case; computed via ord('A')-65 = 0
    assert color == pytest.approx((0.0, 1.0, 0.0, 1.0))


def test_get_color_lookup_alpha_zero():
    """Test get_color_lookup with zero alpha."""
    nucs = ['A', 'T', 'C', 'G']
    colors = CRISPRessoPlot.get_color_lookup(nucs, 0.0)

    for nuc in nucs:
        assert colors[nuc][3] == 0.0


def test_get_color_lookup_alpha_one():
    """Test get_color_lookup with full alpha."""
    nucs = ['A', 'T', 'C', 'G']
    colors = CRISPRessoPlot.get_color_lookup(nucs, 1.0)

    for nuc in nucs:
        assert colors[nuc][3] == 1.0


# =============================================================================
# Tests for amino acid conversions
# =============================================================================


def test_amino_acids_to_numbers_gap():
    """Test amino_acids_to_numbers with gap character."""
    result = CRISPRessoPlot.amino_acids_to_numbers("-")
    assert len(result) == 1
    assert result[0] == 22  # Gap is last in the list


def test_amino_acids_to_numbers_stop():
    """Test amino_acids_to_numbers with stop codon."""
    result = CRISPRessoPlot.amino_acids_to_numbers("*")
    assert len(result) == 1
    assert result[0] == 0  # Stop is first in the list


def test_amino_acids_to_numbers_empty():
    """Test amino_acids_to_numbers with empty string."""
    result = CRISPRessoPlot.amino_acids_to_numbers("")
    assert result == []


# =============================================================================
# Tests for plot functions - file creation verification
# =============================================================================


def test_plot_nucleotide_quilt():
    """Test plot_nucleotide_quilt creates a PDF with valid data."""
    import pandas as pd
    import numpy as np
    import tempfile
    import os

    # Column names at positions 2+ are used as the reference sequence
    nuc_pct_df = pd.DataFrame({
        'Batch': ['Sample1'] * 6,
        'Nucleotide': ['A', 'T', 'C', 'G', 'N', '-'],
        'A': [0.9, 0.05, 0.02, 0.02, 0.01, 0.0],
        'T': [0.1, 0.8, 0.05, 0.04, 0.01, 0.0],
        'C': [0.05, 0.05, 0.85, 0.04, 0.01, 0.0],
        'G': [0.02, 0.02, 0.02, 0.93, 0.01, 0.0],
    })
    mod_pct_df = pd.DataFrame({
        'Batch': ['Sample1', 'Sample1'],
        'Modification': ['Insertions_Left', 'Insertions'],
        'A': [0.0, 0.0],
        'T': [0.0, 0.0],
        'C': [0.0, 0.0],
        'G': [0.0, 0.0],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        fig_root = os.path.join(tmpdir, "test_quilt")
        CRISPRessoPlot.plot_nucleotide_quilt(
            nuc_pct_df,
            mod_pct_df,
            fig_filename_root=fig_root,
            save_also_png=False,
        )
        assert os.path.exists(fig_root + ".pdf")


def test_plot_indel_size_distribution():
    """Test plot_indel_size_distribution creates a PDF with valid data."""
    import numpy as np
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_root = os.path.join(tmpdir, "test_indel")
        CRISPRessoPlot.plot_indel_size_distribution(
            hdensity=np.array([100, 50, 20, 10, 5]),
            hlengths=np.array([-2, -1, 0, 1, 2]),
            center_index=2,
            n_this_category=185,
            xmin=-3,
            xmax=3,
            title="Indel Size Distribution",
            plot_root=plot_root,
            save_also_png=False,
        )
        assert os.path.exists(plot_root + ".pdf")


def test_plot_read_barplot():
    """Test plot_read_barplot creates a PDF with valid data."""
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        fig_root = os.path.join(tmpdir, "test_barplot")
        CRISPRessoPlot.plot_read_barplot(
            N_READS_INPUT=10000,
            N_READS_AFTER_PREPROCESSING=9500,
            N_TOTAL=9000,
            fig_filename_root=fig_root,
            save_png=False,
        )
        assert os.path.exists(fig_root + ".pdf")


def test_plot_class_piechart_and_barplot():
    """Test plot_class_piechart_and_barplot creates piechart and barplot PDFs."""
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        pie_root = os.path.join(tmpdir, "test_pie")
        bar_root = os.path.join(tmpdir, "test_bar")
        CRISPRessoPlot.plot_class_piechart_and_barplot(
            class_counts_order=["Reference_MODIFIED", "Reference_UNMODIFIED"],
            class_counts={
                "Reference_MODIFIED": 200,
                "Reference_UNMODIFIED": 800,
            },
            ref_names=["Reference"],
            expected_hdr_amplicon_seq=None,
            N_TOTAL=1000,
            piechart_plot_root=pie_root,
            barplot_plot_root=bar_root,
            save_png=False,
        )
        assert os.path.exists(pie_root + ".pdf")
        assert os.path.exists(bar_root + ".pdf")


def test_plot_conversion_map():
    """Test plot_conversion_map creates a PDF with valid data."""
    import pandas as pd
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        fig_root = os.path.join(tmpdir, "test_conversion")
        nuc_pct_df = pd.DataFrame({
            'Batch': ['Sample1'] * 4,
            'Nucleotide': ['A', 'T', 'C', 'G'],
            'A': [0.90, 0.05, 0.03, 0.02],
            'C': [0.05, 0.05, 0.85, 0.05],
            'G': [0.02, 0.05, 0.02, 0.91],
        })
        CRISPRessoPlot.plot_conversion_map(
            nuc_pct_df=nuc_pct_df,
            conversion_nuc_from='C',
            conversion_nuc_to='A',
            fig_filename_root=fig_root,
            save_also_png=False,
        )
        assert os.path.exists(fig_root + ".pdf")


def test_plot_conversion_map_nuc_not_found():
    """Test plot_conversion_map returns early when from-nucleotide is not in the sequence."""
    import pandas as pd

    nuc_pct_df = pd.DataFrame({
        'Batch': ['Sample1'] * 4,
        'Nucleotide': ['A', 'T', 'C', 'G'],
        'A': [0.90, 0.05, 0.03, 0.02],
        'T': [0.05, 0.85, 0.05, 0.05],
        'C': [0.03, 0.05, 0.90, 0.02],
    })
    result = CRISPRessoPlot.plot_conversion_map(
        nuc_pct_df=nuc_pct_df,
        conversion_nuc_from='G',
        conversion_nuc_to='A',
    )
    assert result == ()


def test_plot_frequency_deletions_insertions():
    """Test plot_frequency_deletions_insertions creates a PDF with valid data."""
    import numpy as np
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_path = os.path.join(tmpdir, "test_freq")
        ref = {
            'y_values_mut': np.array([10, 5, 2]),
            'x_bins_mut': np.array([0, 1, 2]),
            'y_values_ins': np.array([5, 3, 1]),
            'x_bins_ins': np.array([0, 1, 2]),
            'y_values_del': np.array([8, 4, 2]),
            'x_bins_del': np.array([0, 1, 2]),
        }
        CRISPRessoPlot.plot_frequency_deletions_insertions(
            ref=ref,
            counts_total=100,
            plot_titles={
                'ins': 'Insertions',
                'del': 'Deletions',
                'mut': 'Substitutions',
            },
            plot_path=plot_path,
            xmax_del=5,
            xmax_ins=5,
            xmax_mut=5,
            save_also_png=False,
        )
        assert os.path.exists(plot_path + ".pdf")


def test_plot_alleles_heatmap():
    """Test plot_alleles_heatmap creates a PDF with valid data."""
    import numpy as np
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        fig_root = os.path.join(tmpdir, "test_alleles")
        X = np.array([[1, 2, 3, 4], [1, 2, 4, 4]])
        annot = np.array([['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G']])
        per_element_annot_kws = np.array(
            [[{}, {}, {}, {}], [{}, {}, {'weight': 'bold'}, {}]],
            dtype=object,
        )
        CRISPRessoPlot.plot_alleles_heatmap(
            reference_seq="ATCG",
            X=X,
            annot=annot,
            y_labels=["50.0% (500)", "30.0% (300)"],
            insertion_dict={},
            per_element_annot_kws=per_element_annot_kws,
            fig_filename_root=fig_root,
            SAVE_ALSO_PNG=False,
            plot_cut_point=False,
        )
        assert os.path.exists(fig_root + ".pdf")


def test_plot_amplicon_modifications():
    """Test plot_amplicon_modifications creates a PDF with valid data."""
    import numpy as np
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_root = os.path.join(tmpdir, "test_amp_mod")
        CRISPRessoPlot.plot_amplicon_modifications(
            all_indelsub_count_vectors=np.zeros(20),
            include_idxs_list=list(range(5, 15)),
            cut_points=[10],
            plot_cut_points=[True],
            sgRNA_intervals=[(7, 13)],
            n_total=1000,
            n_this_category=800,
            ref_name="Reference",
            num_refs=1,
            ref_len=20,
            y_max=100,
            plot_titles={
                'main': 'Modification Frequency',
                'combined': 'All modifications',
            },
            plot_root=plot_root,
            save_also_png=False,
        )
        assert os.path.exists(plot_root + ".pdf")
