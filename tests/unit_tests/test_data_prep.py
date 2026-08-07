"""Tests for CRISPResso2.plots.data_prep — extracted data preparation functions.

Each prep function takes a CorePlotContext as its only argument and returns
the kwargs dict that the corresponding CRISPRessoPlot function expects.
"""

import os

import numpy as np
import pandas as pd
from collections import Counter
from types import SimpleNamespace

from CRISPResso2.plots.data_prep import (
    _prep_windowed_alleles,
    amino_acids_to_numbers,
    get_base_edit_target_sequence,
    get_bp_substitutions,
    get_refpos_values,
    get_upset_plot_counts,
    prep_alleles_around_cut,
    prep_alleles_table,
    prep_alleles_table_compare,
    prep_alternate_allele_counts,
    prep_amino_acid_table,
    prep_amino_acid_table_for_plot,
    prep_amplicon_modifications,
    prep_class_piechart_and_barplot,
    prep_base_edit_quilt,
    prep_conversion_at_sel_nucs,
    prep_dsODN_piechart,
    prep_frequency_deletions_insertions,
    prep_global_frameshift_data,
    prep_global_modifications_reference,
    prep_hdr_nucleotide_quilt,
    prep_indel_size_distribution,
    prep_modification_frequency,
    prep_nucleotide_quilt,
    prep_pe_nucleotide_quilt,
    prep_pe_nucleotide_quilt_around_sgRNA,
    write_base_edit_counts,
    write_core_amplicon_data_files,
    _to_numeric_ignore_columns,
)
from CRISPResso2.plots.plot_context import CorePlotContext


# =============================================================================
# Helper: Build minimal CorePlotContext for tests
# =============================================================================


def _make_ctx(**overrides):
    """Build a minimal CorePlotContext, merging *overrides* into defaults."""
    defaults = dict(
        args=SimpleNamespace(
            plot_histogram_outliers=False,
            plot_window_size=20,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
            conversion_nuc_from='C',
            expected_hdr_amplicon_seq='',
            base_editor_output=False,
            coding_seq='',
        ),
        run_data={'running_info': {}},
        output_directory='',
        save_png=False,
        _jp=lambda f: f,
        custom_config={},
        refs={},
        ref_names=[],
        counts_total={},
        counts_modified={},
        counts_unmodified={},
        counts_discarded={},
        counts_insertion={},
        counts_deletion={},
        counts_substitution={},
        class_counts={},
        N_TOTAL=0,
        df_alleles=None,
        all_insertion_count_vectors={},
        all_insertion_left_count_vectors={},
        all_deletion_count_vectors={},
        all_substitution_count_vectors={},
        all_indelsub_count_vectors={},
        all_substitution_base_vectors={},
        all_base_count_vectors={},
        insertion_count_vectors={},
        deletion_count_vectors={},
        substitution_count_vectors={},
        insertion_length_vectors={},
        deletion_length_vectors={},
        hists_frameshift={},
        hists_inframe={},
        counts_modified_frameshift={},
        counts_modified_non_frameshift={},
        counts_non_modified_non_frameshift={},
        counts_splicing_sites_modified={},
    )
    defaults.update(overrides)
    return CorePlotContext(**defaults)


def _ref_dict(**overrides):
    """Build a minimal per-reference dict, merging *overrides* into defaults."""
    defaults = dict(
        sequence='ACGT',
        sequence_length=4,
        ref_plot_name='',
        include_idxs=[0, 1, 2, 3],
        sgRNA_cut_points=[2],
        sgRNA_plot_cut_points=[2],
        sgRNA_intervals=[(1, 3)],
        sgRNA_names=[''],
        sgRNA_mismatches=[],
        sgRNA_sequences=['ACGT'],
        sgRNA_orig_sequences=['ACGT'],
        sgRNA_plot_idxs=[[0, 1, 2, 3]],
        hdensity=np.array([100]),
        hlengths=np.array([0]),
        center_index=0,
        x_bins_ins=np.array([0, 1]),
        y_values_ins=np.array([90, 10]),
        x_bins_del=np.array([0, 1]),
        y_values_del=np.array([80, 20]),
        x_bins_mut=np.array([0, 1]),
        y_values_mut=np.array([85, 15]),
        contains_coding_seq=False,
        exon_intervals=[],
        exon_positions=[],
        exon_len_mods=[],
    )
    defaults.update(overrides)
    return defaults


def _to_serializable(obj):
    """Recursively convert numpy types to plain Python for snapshotting."""
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in sorted(obj.items())}
    if isinstance(obj, (list, tuple)):
        return type(obj)(_to_serializable(v) for v in obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


# =============================================================================
# Tests for _to_numeric_ignore_columns
# =============================================================================


def test_to_numeric_ignore_columns_basic():
    """Test _to_numeric_ignore_columns with basic dataframe."""
    df = pd.DataFrame({
        "name": ["a", "b", "c"],
        "value": ["1", "2", "3"],
        "count": ["10", "20", "30"]
    })
    result = _to_numeric_ignore_columns(df, {"name"})
    assert result["value"].dtype in [int, float, "int64", "float64"]
    assert result["count"].dtype in [int, float, "int64", "float64"]
    assert result["name"].dtype == object


def test_to_numeric_ignore_columns_multiple_ignore():
    """Test _to_numeric_ignore_columns ignoring multiple columns."""
    df = pd.DataFrame({
        "name": ["a", "b"],
        "label": ["x", "y"],
        "value": ["1", "2"]
    })
    result = _to_numeric_ignore_columns(df, {"name", "label"})
    assert result["name"].dtype == object
    assert result["label"].dtype == object
    assert result["value"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_empty_ignore():
    """Test _to_numeric_ignore_columns with empty ignore set."""
    df = pd.DataFrame({
        "a": ["1", "2"],
        "b": ["3", "4"]
    })
    result = _to_numeric_ignore_columns(df, set())
    assert result["a"].dtype in [int, float, "int64", "float64"]
    assert result["b"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_all_numeric():
    """Test _to_numeric_ignore_columns with all numeric columns."""
    df = pd.DataFrame({
        "a": ["1", "2", "3"],
        "b": ["4.5", "5.5", "6.5"]
    })

    result = _to_numeric_ignore_columns(df, set())

    assert result["a"].dtype in [int, float, "int64", "float64"]
    assert result["b"].dtype in [int, float, "int64", "float64"]


def test_to_numeric_ignore_columns_preserve_strings():
    """Test _to_numeric_ignore_columns preserves string columns in ignore list."""
    df = pd.DataFrame({
        "name": ["abc", "def", "ghi"],
        "value": ["1", "2", "3"]
    })

    result = _to_numeric_ignore_columns(df, {"name"})

    assert result["name"].dtype == object
    assert list(result["name"]) == ["abc", "def", "ghi"]


def test_to_numeric_ignore_columns_float_strings():
    """Test _to_numeric_ignore_columns with float strings."""
    df = pd.DataFrame({
        "val": ["1.5", "2.5", "3.5"]
    })

    result = _to_numeric_ignore_columns(df, set())

    assert result["val"].dtype == float


def test_to_numeric_ignore_columns_mixed_types():
    """Test _to_numeric_ignore_columns with mixed data types."""
    df = pd.DataFrame({
        "name": ["a", "b", "c"],
        "int_val": ["1", "2", "3"],
        "float_val": ["1.1", "2.2", "3.3"]
    })

    result = _to_numeric_ignore_columns(df, {"name"})

    assert result["name"].dtype == object
    assert result["int_val"].dtype in [int, float, "int64", "float64"]
    assert result["float_val"].dtype in [float, "float64"]


# =============================================================================
# Tests: prep_amplicon_modifications (plot 4a)
# =============================================================================


class TestPrepAmpliconModifications:

    def test_single_ref(self):
        ctx = _make_ctx(
            ref_names=['ref1'],
            refs={'ref1': _ref_dict(
                sequence='ACGTX',
                sequence_length=5,
                include_idxs=[1, 2, 3],
                sgRNA_cut_points=[2],
                sgRNA_plot_cut_points=[2],
                sgRNA_intervals=[(1, 3)],
            )},
            all_indelsub_count_vectors={'ref1': np.array([0, 1, 2, 3, 4])},
            counts_total={'ref1': 80},
            N_TOTAL=100,
        )
        ctx.ref_name = 'ref1'
        result = prep_amplicon_modifications(ctx)

        assert result['y_max'] == 4 * 1.1
        assert result['n_total'] == 100
        assert result['n_this_category'] == 80
        assert result['ref_len'] == 5
        assert result['plot_titles']['main'] == 'Mutation position distribution'

    def test_multi_ref_title(self):
        ctx = _make_ctx(
            ref_names=['myref', 'other'],
            refs={
                'myref': _ref_dict(),
                'other': _ref_dict(),
            },
            all_indelsub_count_vectors={'myref': np.array([1])},
            counts_total={'myref': 10},
            N_TOTAL=10,
        )
        ctx.ref_name = 'myref'
        result = prep_amplicon_modifications(ctx)
        assert 'myref' in result['plot_titles']['main']


# =============================================================================
# Tests: prep_modification_frequency (plot 4b)
# =============================================================================


class TestPrepModificationFrequency:

    def test_y_max_from_indelsub(self):
        """y_max is computed internally from all_indelsub_count_vectors."""
        ctx = _make_ctx(
            ref_names=['ref1'],
            refs={'ref1': _ref_dict(
                sequence='ACG',
                sequence_length=3,
                include_idxs=[0, 1, 2],
                sgRNA_intervals=[(0, 2)],
            )},
            all_insertion_count_vectors={'ref1': np.array([0, 1, 0])},
            all_deletion_count_vectors={'ref1': np.array([1, 0, 0])},
            all_substitution_count_vectors={'ref1': np.array([0, 0, 1])},
            all_indelsub_count_vectors={'ref1': np.array([1, 1, 1])},
            counts_total={'ref1': 40},
            N_TOTAL=50,
        )
        ctx.ref_name = 'ref1'
        result = prep_modification_frequency(ctx)
        assert result['y_max'] == 1 * 1.1
        assert result['n_total'] == 50
        assert result['n_this_category'] == 40


# =============================================================================
# Tests: prep_indel_size_distribution (plot 3a)
# =============================================================================


class TestPrepIndelSizeDistribution:

    def test_clipping_without_outliers(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                hdensity=np.array([1, 5, 94, 1]),
                hlengths=np.array([-30, -5, 0, 30]),
                center_index=2,
            )},
            counts_total={'r': 101},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['xmin'] == -15
        assert result['xmax'] == 15
        assert result['title'] == 'Indel size distribution'

    def test_no_clipping_with_outliers(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                hdensity=np.array([1, 5, 90, 3, 1]),
                hlengths=np.array([-50, -10, 0, 10, 50]),
                center_index=2,
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=True),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['xmin'] == -50
        assert result['xmax'] == 50

    def test_xmin_xmax_clamped_to_15(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                hdensity=np.array([50, 50]),
                hlengths=np.array([-2, 2]),
                center_index=0,
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=True),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['xmin'] == -15
        assert result['xmax'] == 15

    def test_multi_ref_title(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(hdensity=np.array([100]), hlengths=np.array([0])),
                'HDR': _ref_dict(),
            },
            counts_total={'FANC': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'FANC'
        result = prep_indel_size_distribution(ctx)
        assert result['title'] == 'Indel size distribution: FANC'


# =============================================================================
# Tests: prep_frequency_deletions_insertions (plot 3b)
# =============================================================================


class TestPrepFrequencyDeletionsInsertions:

    def test_clipping(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                x_bins_ins=np.arange(101),
                y_values_ins=np.zeros(101),
                x_bins_del=np.arange(101),
                y_values_del=np.zeros(101),
                x_bins_mut=np.arange(101),
                y_values_mut=np.zeros(101),
                hdensity=np.zeros(100),
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=True),
        )
        ctx.ref_name = 'r'
        result = prep_frequency_deletions_insertions(ctx)
        assert result['xmax_ins'] == 100
        assert result['xmax_del'] == 100

    def test_xmax_clamped_to_15(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                x_bins_ins=np.arange(5),
                y_values_ins=np.array([10, 5, 3, 1, 1]),
                x_bins_del=np.arange(5),
                y_values_del=np.array([10, 5, 3, 1, 1]),
                x_bins_mut=np.arange(5),
                y_values_mut=np.array([10, 5, 3, 1, 1]),
                hdensity=np.ones(20),
            )},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'r'
        result = prep_frequency_deletions_insertions(ctx)
        assert result['xmax_ins'] >= 15
        assert result['xmax_del'] >= 15


# =============================================================================
# Tests: prep_dsODN_piechart (plot 1d)
# =============================================================================


class TestPrepDsODNPiechart:

    def test_basic(self):
        df = pd.DataFrame({
            '#Reads': [60, 40],
            'contains dsODN': [True, False],
        })
        ctx = _make_ctx(df_alleles=df, N_TOTAL=100)
        result = prep_dsODN_piechart(ctx)
        assert len(result['labels']) == 2
        assert abs(result['sizes'][0] - 60.0) < 0.01
        assert abs(result['sizes'][1] - 40.0) < 0.01


# =============================================================================
# Tests: prep_nucleotide_quilt (plot 2a)
# =============================================================================


class TestPrepNucleotideQuilt:

    def test_basic(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence='AC',
                sequence_length=2,
                include_idxs=[0, 1],
            )},
            all_base_count_vectors={
                'r_A': np.array([80, 10]),
                'r_C': np.array([10, 80]),
                'r_G': np.array([5, 5]),
                'r_T': np.array([5, 5]),
                'r_N': np.array([0, 0]),
                'r_-': np.array([0, 0]),
            },
            all_insertion_count_vectors={'r': np.array([2, 3])},
            all_insertion_left_count_vectors={'r': np.array([1, 2])},
            all_deletion_count_vectors={'r': np.array([1, 1])},
            all_substitution_count_vectors={'r': np.array([3, 3])},
            all_indelsub_count_vectors={'r': np.array([5, 5])},
            counts_total={'r': 100},
        )
        ctx.ref_name = 'r'
        result = prep_nucleotide_quilt(ctx)

        assert 'nuc_pct_df' in result
        assert 'mod_pct_df' in result
        assert result['nuc_pct_df'].shape[0] == 6  # A,C,G,T,N,-
        assert result['mod_pct_df'].shape[0] == 6  # Ins, InsL, Del, Sub, All, Total


# =============================================================================
# Tests: prep_global_frameshift_data
# =============================================================================


class TestPrepGlobalFrameshiftData:

    def test_aggregation(self):
        ctx = _make_ctx(
            ref_names=['ref1', 'ref2'],
            refs={
                'ref1': _ref_dict(contains_coding_seq=True),
                'ref2': _ref_dict(contains_coding_seq=False),
            },
            counts_modified_frameshift={'ref1': 10, 'ref2': 5},
            counts_modified_non_frameshift={'ref1': 20, 'ref2': 3},
            counts_non_modified_non_frameshift={'ref1': 60, 'ref2': 80},
            counts_splicing_sites_modified={'ref1': 2, 'ref2': 1},
            counts_total={'ref1': 100, 'ref2': 100},
            counts_modified={'ref1': 30, 'ref2': 8},
            counts_unmodified={'ref1': 70, 'ref2': 92},
            hists_frameshift={'ref1': Counter({0: 5, 3: 10}), 'ref2': Counter({0: 80})},
            hists_inframe={'ref1': Counter({0: 60, 3: 20}), 'ref2': Counter({0: 92})},
        )
        result = prep_global_frameshift_data(ctx)
        # Only ref1 has coding seq
        assert result['global_modified_frameshift'] == 10
        assert result['global_modified_non_frameshift'] == 20
        assert result['global_count_total'] == 100

    def test_hdr_adds_unmodified(self):
        ctx = _make_ctx(
            ref_names=['HDR'],
            refs={'HDR': _ref_dict(contains_coding_seq=True)},
            counts_modified_frameshift={'HDR': 5},
            counts_modified_non_frameshift={'HDR': 10},
            counts_non_modified_non_frameshift={'HDR': 20},
            counts_splicing_sites_modified={'HDR': 1},
            counts_total={'HDR': 100},
            counts_modified={'HDR': 15},
            counts_unmodified={'HDR': 85},
            hists_frameshift={'HDR': Counter({0: 10})},
            hists_inframe={'HDR': Counter({0: 20})},
        )
        result = prep_global_frameshift_data(ctx)
        # HDR unmodified reads added to non_modified_non_frameshift
        assert result['global_non_modified_non_frameshift'] == 20 + 85


# =============================================================================
# Tests: prep_hdr_nucleotide_quilt (plot 4g)
# =============================================================================


class TestPrepHdrNucleotideQuilt:

    def test_basic(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(sequence='AC', sequence_length=2),
                'HDR': _ref_dict(sequence='AC', sequence_length=2),
            },
            counts_total={'FANC': 100, 'HDR': 50},
            ref1_all_base_count_vectors={
                'FANC_A': np.array([80, 10]),
                'FANC_C': np.array([10, 80]),
                'FANC_G': np.array([5, 5]),
                'FANC_T': np.array([5, 5]),
                'FANC_N': np.array([0, 0]),
                'FANC_-': np.array([0, 0]),
                'HDR_A': np.array([40, 5]),
                'HDR_C': np.array([5, 40]),
                'HDR_G': np.array([2, 2]),
                'HDR_T': np.array([3, 3]),
                'HDR_N': np.array([0, 0]),
                'HDR_-': np.array([0, 0]),
            },
            ref1_all_insertion_count_vectors={'FANC': np.array([1, 2]), 'HDR': np.array([0, 1])},
            ref1_all_insertion_left_count_vectors={'FANC': np.array([0, 1]), 'HDR': np.array([0, 0])},
            ref1_all_deletion_count_vectors={'FANC': np.array([1, 0]), 'HDR': np.array([0, 0])},
            ref1_all_substitution_count_vectors={'FANC': np.array([2, 1]), 'HDR': np.array([1, 0])},
            ref1_all_indelsub_count_vectors={'FANC': np.array([3, 3]), 'HDR': np.array([1, 1])},
        )
        result = prep_hdr_nucleotide_quilt(ctx)
        assert 'nuc_pct_df' in result
        assert 'mod_pct_df' in result
        assert result['quantification_window_idxs'] == []  # HDR default


# =============================================================================
# Tests: prep_pe_nucleotide_quilt (plot 11a)
# =============================================================================


class TestPrepPeNucleotideQuilt:

    def test_includes_quantification_window(self):
        ctx = _make_ctx(
            ref_names=['WT', 'PE'],
            refs={
                'WT': _ref_dict(sequence='AC', sequence_length=2, include_idxs=[0, 1]),
                'PE': _ref_dict(sequence='AC', sequence_length=2),
            },
            counts_total={'WT': 100, 'PE': 50},
            ref1_all_base_count_vectors={
                'WT_A': np.array([80, 10]), 'WT_C': np.array([10, 80]),
                'WT_G': np.array([5, 5]), 'WT_T': np.array([5, 5]),
                'WT_N': np.array([0, 0]), 'WT_-': np.array([0, 0]),
                'PE_A': np.array([40, 5]), 'PE_C': np.array([5, 40]),
                'PE_G': np.array([2, 2]), 'PE_T': np.array([3, 3]),
                'PE_N': np.array([0, 0]), 'PE_-': np.array([0, 0]),
            },
            ref1_all_insertion_count_vectors={'WT': np.array([1, 2]), 'PE': np.array([0, 1])},
            ref1_all_insertion_left_count_vectors={'WT': np.array([0, 1]), 'PE': np.array([0, 0])},
            ref1_all_deletion_count_vectors={'WT': np.array([1, 0]), 'PE': np.array([0, 0])},
            ref1_all_substitution_count_vectors={'WT': np.array([2, 1]), 'PE': np.array([1, 0])},
            ref1_all_indelsub_count_vectors={'WT': np.array([3, 3]), 'PE': np.array([1, 1])},
        )
        result = prep_pe_nucleotide_quilt(ctx)
        # PE uses first ref's include_idxs
        assert result['quantification_window_idxs'] == [0, 1]


# =============================================================================
# Tests: _prep_windowed_alleles (private helper, unchanged)
# =============================================================================


class TestPrepWindowedAlleles:

    def test_basic_windowing(self):
        df = pd.DataFrame({
            'Aligned_Sequence': ['ACG', 'AXG'],
            'Reference_Sequence': ['ACG', 'ACG'],
            '#Reads': [80, 20],
            '%Reads': [80.0, 20.0],
        }).set_index('Aligned_Sequence')

        df_out, df_plot, ref_seq, intervals, start = _prep_windowed_alleles(
            df_alleles_around_cut=df,
            cut_point=1,
            window_left=1,
            window_right=1,
            ref_sequence='AACGG',
            sgRNA_intervals=[(0, 4)],
            count_total=100,
            allele_plot_pcts_only_for_assigned_reference=False,
            expand_allele_plots_by_quantification=True,
        )
        assert ref_seq == 'AC'
        assert len(intervals) == 1


# =============================================================================
# Tests: prep_class_piechart_and_barplot (utility, not CorePlotContext-based)
# =============================================================================


class TestPrepClassPiechartAndBarplot:

    def test_basic_single_ref(self):
        result = prep_class_piechart_and_barplot(
            class_counts_order=['Reference_UNMODIFIED', 'Reference_MODIFIED'],
            class_counts={'Reference_UNMODIFIED': 80, 'Reference_MODIFIED': 20},
            ref_names=['Reference'],
            expected_hdr_amplicon_seq='',
            N_TOTAL=100,
        )
        assert 'UNMODIFIED' in result['labels'][0]
        assert abs(result['sizes'][0] - 80.0) < 0.01

    def test_hdr_mode_labels(self):
        result = prep_class_piechart_and_barplot(
            class_counts_order=['FANC_MODIFIED', 'HDR_MODIFIED', 'HDR_UNMODIFIED'],
            class_counts={'FANC_MODIFIED': 10, 'HDR_MODIFIED': 5, 'HDR_UNMODIFIED': 85},
            ref_names=['FANC', 'HDR'],
            expected_hdr_amplicon_seq='ACGT',
            N_TOTAL=100,
        )
        assert 'NHEJ' in result['labels'][0]
        assert 'Imperfect HDR' in result['labels'][1]
        assert 'HDR' in result['labels'][2]

    def test_return_keys(self):
        result = prep_class_piechart_and_barplot(
            class_counts_order=['R_UNMODIFIED'],
            class_counts={'R_UNMODIFIED': 50},
            ref_names=['R'],
            expected_hdr_amplicon_seq='',
            N_TOTAL=50,
        )
        assert sorted(result.keys()) == ['N_TOTAL', 'labels', 'sizes']


# =============================================================================
# Tests: prep_alternate_allele_counts (utility, not CorePlotContext-based)
# =============================================================================


class TestPrepAlternateAlleleCounts:

    def test_basic(self):
        sub_base_vectors = {
            'R_A': np.array([0, 0, 0]),
            'R_C': np.array([0, 90, 0]),
            'R_G': np.array([0, 0, 95]),
            'R_T': np.array([0, 5, 0]),
            'R_N': np.array([0, 0, 0]),
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', 'ACG')
        assert result['C']['T'] == 5
        assert result['C']['C'] == 90
        assert result['G']['G'] == 95

    def test_all_zeros(self):
        sub_base_vectors = {
            'R_A': np.array([0, 0]),
            'R_C': np.array([0, 0]),
            'R_G': np.array([0, 0]),
            'R_T': np.array([0, 0]),
            'R_N': np.array([0, 0]),
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', 'AC')
        for a in 'ACGTN':
            for b in 'ACGTN':
                assert result[a][b] == 0

    def test_multiple_positions_same_base(self):
        sub_base_vectors = {
            'R_A': np.array([3, 2]),
            'R_C': np.array([85, 90]),
            'R_G': np.array([0, 0]),
            'R_T': np.array([12, 8]),
            'R_N': np.array([0, 0]),
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', 'CC')
        assert result['C']['A'] == 5   # 3 + 2
        assert result['C']['T'] == 20  # 12 + 8
        assert result['C']['C'] == 175 # 85 + 90

    def test_list_ref_sequence(self):
        sub_base_vectors = {
            'R_A': [10, 0],
            'R_C': [0, 5],
            'R_G': [0, 0],
            'R_T': [0, 0],
            'R_N': [0, 0],
        }
        result = prep_alternate_allele_counts(sub_base_vectors, 'R', ['A', 'C'])
        assert result['A']['A'] == 10
        assert result['C']['C'] == 5


# =============================================================================
# Tests: prep_global_modifications_reference (plot 4e/4f)
# =============================================================================


class TestPrepGlobalModificationsReference:

    def test_ref0_is_4e(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(sequence_length=4, include_idxs=[0, 1, 2, 3]),
                'HDR': _ref_dict(),
            },
            N_TOTAL=100,
            ref1_all_insertion_count_vectors={'FANC': np.array([1, 0, 0, 0])},
            ref1_all_deletion_count_vectors={'FANC': np.array([0, 1, 0, 0])},
            ref1_all_substitution_count_vectors={'FANC': np.array([0, 0, 1, 0])},
        )
        ctx.ref_name = 'FANC'
        result = prep_global_modifications_reference(ctx)
        assert '4e' in result['fig_filename_root']
        assert 'all reads' in result['plot_title']

    def test_hdr_is_4f(self):
        ctx = _make_ctx(
            ref_names=['FANC', 'HDR'],
            refs={
                'FANC': _ref_dict(sequence_length=4, include_idxs=[0, 1, 2, 3]),
                'HDR': _ref_dict(sequence_length=4, include_idxs=[0, 1, 2, 3]),
            },
            N_TOTAL=100,
            ref1_all_insertion_count_vectors={'HDR': np.array([1, 0, 0, 0])},
            ref1_all_deletion_count_vectors={'HDR': np.array([0, 1, 0, 0])},
            ref1_all_substitution_count_vectors={'HDR': np.array([0, 0, 1, 0])},
        )
        ctx.ref_name = 'HDR'
        result = prep_global_modifications_reference(ctx)
        assert '4f' in result['fig_filename_root']
        assert 'HDR' in result['plot_title']


# =============================================================================
# Tests: prep_conversion_at_sel_nucs (plots 10e/10f/10g)
# =============================================================================


class TestPrepConversionAtSelNucs:

    def test_basic(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence='ACCA',
                sequence_length=4,
                sgRNA_plot_idxs=[[0, 1, 2, 3]],
            )},
            all_base_count_vectors={
                'r_A': np.array([80, 5, 5, 80]),
                'r_C': np.array([10, 85, 85, 10]),
                'r_G': np.array([5, 5, 5, 5]),
                'r_T': np.array([5, 5, 5, 5]),
                'r_N': np.array([0, 0, 0, 0]),
                'r_-': np.array([0, 0, 0, 0]),
            },
            counts_total={'r': 100},
            args=SimpleNamespace(conversion_nuc_from='C'),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        result = prep_conversion_at_sel_nucs(ctx)

        # C appears at positions 1, 2 in 'ACCA'
        assert len(result['from_nuc_indices']) == 2
        assert result['just_sel_nuc_pcts'].shape[1] == 2
        assert result['just_sel_nuc_freqs'].shape[1] == 2


# =============================================================================
# Tests: plot root generation
# =============================================================================


class TestPlotRootGeneration:
    """Verify that fig_filename_root/fig_filename_root is generated from ctx._jp."""

    def test_indel_size_uses_jp(self):
        def fake_jp(name):
            return '/output/' + name

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(ref_plot_name='R.')},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
            _jp=fake_jp,
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['fig_filename_root'] == '/output/3a.R.Indel_size_distribution'

    def test_no_jp_falls_back_to_filename(self):
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(ref_plot_name='')},
            counts_total={'r': 100},
            args=SimpleNamespace(plot_histogram_outliers=False),
        )
        ctx.ref_name = 'r'
        result = prep_indel_size_distribution(ctx)
        assert result['fig_filename_root'] == '3a.Indel_size_distribution'


# =============================================================================
# Tests for amino_acids_to_numbers
# =============================================================================


def test_amino_acids_to_numbers_basic():
    """Test amino_acids_to_numbers with basic sequence."""
    result = amino_acids_to_numbers("MA")
    assert result == [11, 1]


def test_amino_acids_to_numbers_stop():
    """Test amino_acids_to_numbers with stop codon."""
    result = amino_acids_to_numbers("*")
    assert len(result) == 1
    assert result[0] == 0  # Stop codon is first in list


def test_amino_acids_to_numbers_gap():
    """Test amino_acids_to_numbers with gap."""
    result = amino_acids_to_numbers("-")
    assert result == [22]


def test_amino_acids_to_numbers_all_standard():
    """Test amino_acids_to_numbers with all standard amino acids."""
    all_aa = "ACDEFGHIKLMNPQRSTVWY"
    result = amino_acids_to_numbers(all_aa)
    assert result == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    # All should be unique
    assert len(set(result)) == 20


def test_amino_acids_to_numbers_empty():
    """Test amino_acids_to_numbers with empty sequence."""
    result = amino_acids_to_numbers("")
    assert result == []


def test_amino_acids_to_numbers_with_special():
    """Test amino_acids_to_numbers with special characters."""
    result = amino_acids_to_numbers("*-")
    assert len(result) == 2
    assert result[0] == 0  # * is first
    assert result[1] == 22  # - is last


# =============================================================================
# Tests for prep_alleles_table
# =============================================================================


def test_prep_alleles_table_basic():
    """Test prep_alleles_table with basic data."""
    df = pd.DataFrame({
        '%Reads': [50.0, 30.0, 20.0],
        '#Reads': [500, 300, 200],
        'Reference_Sequence': ['ATCG', 'ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG', 'A-CG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4], [1, 0, 3, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G'], ['A', '-', 'C', 'G']]
    assert y_labels == ['50.00% (500 reads)', '30.00% (300 reads)', '20.00% (200 reads)']
    assert is_reference == [True, False, False]


def test_prep_alleles_table_empty():
    """Test prep_alleles_table with empty dataframe after filtering."""
    df = pd.DataFrame({
        '%Reads': [0.1],
        '#Reads': [1],
        'Reference_Sequence': ['ATCG'],
    }, index=['ATCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=1.0)

    assert len(X) == 0
    assert len(annot) == 0


def test_prep_alleles_table_max_rows():
    """Test prep_alleles_table respects MAX_N_ROWS."""
    df = pd.DataFrame({
        '%Reads': [30.0, 25.0, 20.0, 15.0, 10.0],
        '#Reads': [300, 250, 200, 150, 100],
        'Reference_Sequence': ['ATCG', 'ATCG', 'ATCG', 'ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG', 'TTCG', 'ATCA', 'GGGG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=3, MIN_FREQUENCY=0)

    assert X == [[1, 2, 3, 4], [1, 2, 4, 4], [2, 2, 3, 4]]
    assert annot == [['A', 'T', 'C', 'G'], ['A', 'T', 'G', 'G'], ['T', 'T', 'C', 'G']]
    assert y_labels == ['30.00% (300 reads)', '25.00% (250 reads)', '20.00% (200 reads)']


def test_prep_alleles_table_with_insertions():
    """Test prep_alleles_table detects insertions."""
    # Reference with gap indicates insertion in read
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['AT--CG'],
    }, index=['ATGGCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATGGCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    # Should detect insertion
    assert 0 in insertion_dict
    assert len(insertion_dict[0]) > 0


def test_prep_alleles_table_with_substitution():
    """Test prep_alleles_table detects substitutions."""
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['ATCG'],  # Reference
    }, index=['GTCG'])  # G at position 0 instead of A

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(X) == 1
    assert is_reference[0] is False  # Different from reference
    # Should have bold annotation for substitution
    assert len(per_element_annot_kws[0]) > 0


def test_prep_alleles_table_all_reference():
    """Test prep_alleles_table with all reference sequences."""
    df = pd.DataFrame({
        '%Reads': [100.0],
        '#Reads': [1000],
        'Reference_Sequence': ['ATCG'],
    }, index=['ATCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = \
        prep_alleles_table(df, 'ATCG', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(is_reference) == 1
    assert is_reference[0] is True


# =============================================================================
# Tests for prep_alleles_table_compare
# =============================================================================


def test_prep_alleles_table_compare_basic():
    """Test prep_alleles_table_compare with basic data."""
    # Create merged allele dataframe
    df = pd.DataFrame({
        '%Reads_sample1': [50.0, 30.0],
        '%Reads_sample2': [40.0, 35.0],
        '#Reads_sample1': [500, 300],
        '#Reads_sample2': [400, 350],
        'Reference_Sequence': ['ATCG', 'ATCG'],
    }, index=['ATCG', 'ATGG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws = \
        prep_alleles_table_compare(
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
    df = pd.DataFrame({
        '%Reads_s1': [50.0],
        '%Reads_s2': [50.0],
        '#Reads_s1': [500],
        '#Reads_s2': [500],
        'Reference_Sequence': ['AT--CG'],  # Insertion markers
    }, index=['ATGGCG'])

    X, annot, y_labels, insertion_dict, per_element_annot_kws = \
        prep_alleles_table_compare(
            df, 's1', 's2', MAX_N_ROWS=10, MIN_FREQUENCY=0
        )

    # Should detect insertion
    assert 0 in insertion_dict
    assert len(insertion_dict[0]) > 0


# =============================================================================
# Tests for prep_amino_acid_table_for_plot
# =============================================================================


def test_prep_amino_acid_table_for_plot_basic():
    """Test prep_amino_acid_table_for_plot with basic data."""
    df = pd.DataFrame({
        '%Reads': [60.0, 40.0],
        '#Reads': [600, 400],
        'Reference_Sequence': ['MAS', 'MAS'],
        'silent_edit_inds': [[], []],
    }, index=['MAS', 'MAT'])

    X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq = \
        prep_amino_acid_table_for_plot(df, 'MAS', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    assert len(X) == 2
    assert is_reference[0] is True  # First row matches reference
    assert ref_seq == 'MAS'


def test_prep_amino_acid_table_for_plot_with_silent_edits():
    """Test prep_amino_acid_table_for_plot with silent edits."""
    df = pd.DataFrame({
        '%Reads': [50.0],
        '#Reads': [500],
        'Reference_Sequence': ['MAS'],
        'silent_edit_inds': [[1]],  # Silent edit at position 1
    }, index=['MAS'])

    X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq = \
        prep_amino_acid_table_for_plot(df, 'MAS', MAX_N_ROWS=10, MIN_FREQUENCY=0)

    # Should have silent edit at row 0, position 1
    assert 0 in silent_edit_dict
    assert 1 in silent_edit_dict[0]


# =============================================================================
# Helper: Build a minimal df_alleles for integration-style prep tests
# =============================================================================


def _make_df_alleles(ref_name, aligned_seqs, ref_seq, reads, read_status=None):
    """Build a df_alleles-like DataFrame suitable for CRISPRessoShared cut functions.

    Parameters
    ----------
    ref_name : str
        Reference name to populate ``Reference_Name``.
    aligned_seqs : list[str]
        Aligned read sequences (may contain ``-`` for deletions).
    ref_seq : str
        The aligned reference sequence (same length as each aligned_seq).
    reads : list[int]
        ``#Reads`` for each allele.
    read_status : list[str] or None
        ``Read_Status`` values; defaults to ``'MODIFIED'`` for all.

    Returns a DataFrame with the columns expected by
    ``CRISPRessoShared.get_dataframe_around_cut_asymmetrical`` et al.
    """
    total_reads = sum(reads)
    n = len(aligned_seqs)
    if read_status is None:
        read_status = ['MODIFIED'] * n

    # Build ref_positions: for each position in the alignment, what is the
    # reference coordinate (0-indexed)?  Gaps in the reference get -1.
    ref_positions = []
    pos = 0
    for c in ref_seq:
        if c == '-':
            ref_positions.append(-1)
        else:
            ref_positions.append(pos)
            pos += 1

    # Compute simple n_deleted / n_inserted / n_mutated per allele
    rows = []
    for i, (aln, status) in enumerate(zip(aligned_seqs, read_status)):
        n_del = sum(1 for a, r in zip(aln, ref_seq) if a == '-' and r != '-')
        n_ins = sum(1 for a, r in zip(aln, ref_seq) if a != '-' and r == '-')
        n_mut = sum(1 for a, r in zip(aln, ref_seq) if a != '-' and r != '-' and a != r)
        rows.append({
            'Aligned_Sequence': aln,
            'Reference_Sequence': ref_seq,
            'ref_positions': ref_positions,
            'Read_Status': status,
            'n_deleted': n_del,
            'n_inserted': n_ins,
            'n_mutated': n_mut,
            '#Reads': reads[i],
            '%Reads': 100.0 * reads[i] / total_reads,
            'Reference_Name': ref_name,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Tests: prep_alleles_around_cut (plot 9)
# =============================================================================


class TestPrepAllelesAroundCut:

    def test_basic(self):
        """Smoke test: single allele, unmodified, returns expected keys."""
        #                       0123456789
        ref_sequence =         'AACCGGTTAA'
        aligned_ref =          'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence],
            ref_seq=aligned_ref,
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)

        assert 'df_alleles_around_cut' in result
        assert 'ref_seq_around_cut' in result
        assert 'new_sgRNA_intervals' in result
        assert 'new_cut_point' in result
        assert 'window_truncated' in result
        assert 'plot_input' in result
        # Window of size 3 around cut_point=4 → positions 2..7 → 6 chars
        assert len(result['ref_seq_around_cut']) <= 7

    def test_returns_sgRNA_legend_for_downstream_captions(self):
        """CRISPRessoPro writes its own caption for the interactive table and
        needs the same legend string, rather than re-deriving the format rule."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence],
            ref_seq=ref_sequence,
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['guide_a'],
                sgRNA_orig_sequences=['ACGT'],
                sgRNA_mismatches=[],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)

        assert result['sgRNA_legend'] == 'guide_a (ACGT)'
        # and the shared caption is built from the same string
        assert 'guide_a (ACGT)' in result['caption']

    def test_quantification_window_idxs_are_remapped_to_window_coordinates(self):
        """The window is a slice of the amplicon, so the quantification window
        indices must be rebased onto it -- the same remapping new_sgRNA_intervals
        gets. Without it the widget would draw the band at amplicon coordinates."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence],
            ref_seq=ref_sequence,
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
                include_idxs=[3, 4, 5],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)

        # new_sel_cols_start = cut_point(4) - window_left(3) = 1, and this path
        # carries an extra -1 (see new_cut_point), so the shift is 2
        assert result['quantification_window_idxs'] == [1, 2, 3]

        # the invariant that actually matters: the quantification window and the
        # sgRNA intervals must be rebased identically, or the band and the guide
        # would disagree about where they are
        shift = 2 - result['new_sgRNA_intervals'][0][0]  # sgRNA_intervals=[(2, 6)]
        assert result['quantification_window_idxs'] == [
            x - shift for x in [3, 4, 5]
        ]

    def test_with_substitution(self):
        """Allele with a substitution produces plot_input with expected content."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA', 'AATCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)

        assert result['plot_input'] is not None
        pi = result['plot_input']
        assert 'reference_seq' in pi
        assert 'prepped_df_alleles' in pi
        assert 'fig_filename_root' in pi
        assert '9.' in pi['fig_filename_root']
        # Two alleles in the window
        assert len(pi['y_labels']) >= 1

    def test_plot_input_none_below_threshold(self):
        """When all alleles are below min_frequency, plot_input is None."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=99999,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        assert result['plot_input'] is None

    def test_window_truncation_near_edge(self):
        """Window near amplicon edge sets window_truncated=True."""
        ref_sequence = 'ACGT'
        df = _make_df_alleles(
            'r', aligned_seqs=['ACGT'], ref_seq='ACGT',
            reads=[100], read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=4,
                sgRNA_cut_points=[0],
                sgRNA_plot_cut_points=[0],
                sgRNA_intervals=[(0, 3)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=20,  # much larger than sequence
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        assert result['window_truncated'] is True


    def test_pct_adjustment_for_assigned_reference(self):
        """allele_plot_pcts_only_for_assigned_reference mutates %Reads in place."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA', 'AATCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
            )},
            counts_total={'r': 50},  # different from sum of reads to test recalc
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=True,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        df_out = result['df_alleles_around_cut']

        # Should have %AllReads (original) and recalculated %Reads
        assert '%AllReads' in df_out.columns
        # %Reads should be recalculated: #Reads / counts_total * 100
        for _, row in df_out.iterrows():
            expected = row['#Reads'] / 50 * 100
            assert abs(row['%Reads'] - expected) < 0.01

    def test_new_cut_point_inside_sgRNA(self):
        """new_cut_point is set when cut_point falls inside an sgRNA interval."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[100],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],  # cut_point 4 is inside this
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        assert result['new_cut_point'] is not None

    def test_new_cut_point_outside_sgRNA(self):
        """new_cut_point is None when cut_point is outside all sgRNA intervals."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[100],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(7, 9)],  # cut_point 4 is NOT inside this
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        assert result['new_cut_point'] is None

    def test_groupby_collapse(self):
        """expand_allele_plots_by_quantification=False collapses alleles by sequence."""
        ref_sequence = 'AACCGGTTAA'
        # Two identical alleles that should collapse
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA', 'AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[60, 40],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=False,
                annotate_wildtype_allele='',
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_alleles_around_cut(ctx)
        # The two identical alleles should be collapsed into one row in plot_input
        pi = result['plot_input']
        if pi is not None:
            assert len(pi['y_labels']) == 1


# =============================================================================
# Tests: prep_base_edit_quilt (plot 10h)
# =============================================================================


class TestPrepBaseEditQuilt:

    def test_basic(self):
        """Smoke test: base edit quilt with one C position."""
        #                 0123456789
        ref_sequence =   'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA', 'AATCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                conversion_nuc_from='C',
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_base_edit_quilt(ctx)

        assert 'df_alleles_around_cut' in result
        assert 'ref_seq_around_cut' in result
        assert 'plot_input' in result

    def test_plot_input_none_below_threshold(self):
        """When all alleles are below min_frequency, plot_input is None."""
        ref_sequence = 'AACCGGTTAA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['AACCGGTTAA'],
            ref_seq='AACCGGTTAA',
            reads=[100],
            read_status=['UNMODIFIED'],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                conversion_nuc_from='C',
                min_frequency_alleles_around_cut_to_plot=99999,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_base_edit_quilt(ctx)
        assert result['plot_input'] is None

    def test_x_labels_are_conversion_nuc_positions(self):
        """x_labels should be 1-indexed positions of conversion_nuc_from."""
        ref_sequence = 'ACCA'
        df = _make_df_alleles(
            'r',
            aligned_seqs=['ACCA', 'ATCA'],
            ref_seq='ACCA',
            reads=[70, 30],
        )
        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=4,
                sgRNA_cut_points=[1],
                sgRNA_plot_cut_points=[1],
                sgRNA_intervals=[(0, 3)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                plot_window_size=3,
                conversion_nuc_from='C',
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
            ),
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0

        result = prep_base_edit_quilt(ctx)
        # 'ACCA' — C at positions 1 and 2 (0-indexed) → 1-indexed: [2, 3]
        pi = result['plot_input']
        if pi is not None:
            assert pi['x_labels'] == [2, 3]


# =============================================================================
# Tests: prep_amino_acid_table (plot 9a)
# =============================================================================


class TestPrepAminoAcidTable:

    def test_basic(self):
        """Smoke test: amino acid table with a simple coding seq."""
        # 9-base coding seq → 3 amino acids
        coding_seq = 'ATGGCTTAA'  # M, A, *
        # Build a ref that contains this coding seq starting at position 0
        ref_sequence = coding_seq
        ref_len = len(ref_sequence)

        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence, 'ATGACTTAA'],
            ref_seq=ref_sequence,
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=ref_len,
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
                contains_coding_seq=True,
                exon_positions=list(range(ref_len)),
                exon_intervals=[(0, ref_len - 1)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='',
            ),
            run_data={'running_info': {'coding_seqs': [coding_seq]}},
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        ctx.coding_seq_ind = 0

        result = prep_amino_acid_table(ctx)

        assert 'coding_seq_amino_acids' in result
        assert 'amino_acid_cut_point' in result
        assert 'df_to_plot' in result
        assert 'plot_input' in result
        # Coding seq ATG GCT TAA → M A *
        assert result['coding_seq_amino_acids'] == 'MA*'

    def test_plot_input_none_below_threshold(self):
        """When no rows pass threshold, plot_input is None."""
        coding_seq = 'ATGGCTTAA'
        ref_sequence = coding_seq

        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence],
            ref_seq=ref_sequence,
            reads=[100],
            read_status=['UNMODIFIED'],
        )

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                contains_coding_seq=True,
                exon_positions=list(range(len(ref_sequence))),
                exon_intervals=[(0, len(ref_sequence) - 1)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                min_frequency_alleles_around_cut_to_plot=99999,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='',
            ),
            run_data={'running_info': {'coding_seqs': [coding_seq]}},
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        ctx.coding_seq_ind = 0

        result = prep_amino_acid_table(ctx)
        assert result['plot_input'] is None

    def test_plot_input_populated_above_threshold(self):
        """When rows pass threshold, plot_input contains expected keys."""
        coding_seq = 'ATGGCTTAA'  # M, A, *
        ref_sequence = coding_seq

        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence, 'ATGACTTAA'],
            ref_seq=ref_sequence,
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
                contains_coding_seq=True,
                exon_positions=list(range(len(ref_sequence))),
                exon_intervals=[(0, len(ref_sequence) - 1)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='',
            ),
            run_data={'running_info': {'coding_seqs': [coding_seq]}},
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        ctx.coding_seq_ind = 0

        result = prep_amino_acid_table(ctx)
        pi = result['plot_input']
        assert pi is not None
        assert 'reference_seq_amino_acids' in pi
        assert 'fig_filename_root' in pi
        assert 'X' in pi
        assert 'annot' in pi
        assert 'y_labels' in pi
        assert pi['reference_seq_amino_acids'] == 'MA*'
        assert '9a.' in pi['fig_filename_root']

    def test_annotate_wildtype_allele(self):
        """annotate_wildtype_allele appends to y_labels of reference rows."""
        coding_seq = 'ATGGCTTAA'  # M, A, *
        ref_sequence = coding_seq

        df = _make_df_alleles(
            'r',
            aligned_seqs=[ref_sequence, 'ATGACTTAA'],
            ref_seq=ref_sequence,
            reads=[80, 20],
            read_status=['UNMODIFIED', 'MODIFIED'],
        )

        ctx = _make_ctx(
            ref_names=['r'],
            refs={'r': _ref_dict(
                sequence=ref_sequence,
                sequence_length=len(ref_sequence),
                sgRNA_cut_points=[4],
                sgRNA_plot_cut_points=[4],
                sgRNA_intervals=[(2, 6)],
                sgRNA_names=['sgRNA1'],
                sgRNA_mismatches=[],
                contains_coding_seq=True,
                exon_positions=list(range(len(ref_sequence))),
                exon_intervals=[(0, len(ref_sequence) - 1)],
            )},
            counts_total={'r': 100},
            df_alleles=df,
            args=SimpleNamespace(
                min_frequency_alleles_around_cut_to_plot=0,
                max_rows_alleles_around_cut_to_plot=10,
                annotate_wildtype_allele='****',
            ),
            run_data={'running_info': {'coding_seqs': [coding_seq]}},
        )
        ctx.ref_name = 'r'
        ctx.sgRNA_ind = 0
        ctx.coding_seq_ind = 0

        result = prep_amino_acid_table(ctx)
        pi = result['plot_input']
        if pi is not None:
            # At least one y_label should end with '****'
            wt_labels = [l for l in pi['y_labels'] if l.endswith('****')]
            assert len(wt_labels) >= 1


# =============================================================================
# Tests: base editing utility functions (get_refpos_values, get_bp_substitutions,
#        get_base_edit_target_sequence, get_upset_plot_counts, write_base_edit_counts)
# =============================================================================

class TestGetRefposValues:

    def test_no_gaps(self):
        """No gaps: each ref position maps directly to its read base."""
        result = get_refpos_values("ATCG", "ATCG")
        assert result[0] == "A"
        assert result[1] == "T"
        assert result[2] == "C"
        assert result[3] == "G"

    def test_gap_in_ref(self):
        """Leading gaps in ref: insertions accumulate on position 0."""
        result = get_refpos_values("--ATGC", "GGATGC")
        assert result[0] == "GGA"
        assert result[1] == "T"
        assert result[2] == "G"
        assert result[3] == "C"

    def test_gap_in_read(self):
        """Deletion in read is represented as '-'."""
        result = get_refpos_values("ATGC", "A-GC")
        assert result[0] == "A"
        assert result[1] == "-"
        assert result[2] == "G"
        assert result[3] == "C"

    def test_example_from_docstring(self):
        """Example straight from the function docstring."""
        result = get_refpos_values("--A-TGC-", "GGAGTCGA")
        assert result[0] == "GGAG"
        assert result[1] == "T"
        assert result[2] == "C"
        assert result[3] == "GA"

    def test_insertion_middle(self):
        """Insertion in middle of ref attaches to preceding ref position."""
        result = get_refpos_values("AT-CG", "ATGCG")
        assert result[0] == "A"
        assert result[1] == "TG"
        assert result[2] == "C"
        assert result[3] == "G"

    def test_deletion(self):
        """Two consecutive deletions each become '-'."""
        result = get_refpos_values("ATCG", "A--G")
        assert result[0] == "A"
        assert result[1] == "-"
        assert result[2] == "-"
        assert result[3] == "G"

    def test_complex(self):
        """Mixed insertion and deletion in same alignment."""
        result = get_refpos_values("A-TC--G", "AGTCAAG")
        assert result[0] == "AG"
        assert result[1] == "T"
        assert result[2] == "CAA"
        assert result[3] == "G"

    def test_all_matches(self):
        """All-match alignment: every position is the identity base."""
        ref = "ATCGATCG"
        result = get_refpos_values(ref, ref)
        for i, base in enumerate(ref):
            assert result[i] == base

    def test_insertion_at_start(self):
        """Insertions before the first ref base accumulate on position 0."""
        result = get_refpos_values("--ATCG", "GGATCG")
        assert result[0] == "GGA"

    def test_insertion_at_end(self):
        """Insertions after the last ref base attach to the last position."""
        result = get_refpos_values("ATCG--", "ATCGGG")
        assert result[3] == "GGG"

    def test_deletions(self):
        """Two deleted positions in the middle are both '-'."""
        result = get_refpos_values("ATCGATCG", "AT--ATCG")
        assert result[2] == "-"
        assert result[3] == "-"

    def test_mixed_indels(self):
        """Mixed insertions and deletions produce a non-empty dict."""
        result = get_refpos_values("AT-CGATCG", "ATG--ATCG")
        assert len(result) > 0


class TestGetBpSubstitutions:

    def test_one_sub(self):
        """Single substitution at position 0."""
        ref_changes_dict = get_refpos_values("AAA", "CAA")
        result = get_bp_substitutions(ref_changes_dict, "AAA", [0, 1, 2])
        assert len(result) == 1
        assert result[0] == (0, 'A', 'C')

    def test_two_subs(self):
        """Two substitutions at positions 0 and 1."""
        ref_changes_dict = get_refpos_values("AAA", "CCA")
        result = get_bp_substitutions(ref_changes_dict, "AAA", [0, 1, 2])
        assert result == [(0, 'A', 'C'), (1, 'A', 'C')]

    def test_insertions(self):
        """Insertions adjacent to ref positions are captured."""
        ref_changes_dict = get_refpos_values("AAA-AAA-----", "AAACAAACCCCC")
        result = get_bp_substitutions(ref_changes_dict, "AAAAAA", [0, 1, 2, 3, 4, 5])
        assert len(result) == 2
        assert result[0] == (2, 'A', 'AC')
        assert result[1] == (5, 'A', 'ACCCCC')

    def test_no_changes(self):
        """Identical sequences produce empty list."""
        ref_changes_dict = get_refpos_values("ATCG", "ATCG")
        assert get_bp_substitutions(ref_changes_dict, "ATCG", [0, 1, 2, 3]) == []

    def test_single_sub(self):
        """A→G at position 0."""
        ref_changes_dict = get_refpos_values("ATCG", "GTCG")
        assert get_bp_substitutions(ref_changes_dict, "ATCG", [0, 1, 2, 3]) == [(0, 'A', 'G')]

    def test_multiple_subs(self):
        """Three substitutions across the sequence."""
        ref_changes_dict = get_refpos_values("ATCG", "GACT")
        assert get_bp_substitutions(ref_changes_dict, "ATCG", [0, 1, 2, 3]) == [
            (0, 'A', 'G'), (1, 'T', 'A'), (3, 'G', 'T'),
        ]

    def test_partial_positions(self):
        """Only positions in ref_positions_to_include are checked."""
        ref_changes_dict = get_refpos_values("ATCG", "GTCA")
        assert get_bp_substitutions(ref_changes_dict, "ATCG", [0, 3]) == [
            (0, 'A', 'G'), (3, 'G', 'A'),
        ]

    def test_all_match(self):
        """All identical: empty result."""
        ref_changes_dict = get_refpos_values("ATCG", "ATCG")
        assert get_bp_substitutions(ref_changes_dict, "ATCG", [0, 1, 2, 3]) == []

    def test_all_different(self):
        """All four positions substituted."""
        ref_changes_dict = get_refpos_values("ATCG", "TAGC")
        assert get_bp_substitutions(ref_changes_dict, "ATCG", [0, 1, 2, 3]) == [
            (0, 'A', 'T'), (1, 'T', 'A'), (2, 'C', 'G'), (3, 'G', 'C'),
        ]

    def test_with_insertion(self):
        """Pre-built dict with an insertion at position 1 is flagged."""
        ref_changes_dict = {0: 'A', 1: 'TG', 2: 'C', 3: 'G'}
        result = get_bp_substitutions(ref_changes_dict, "ATCG", [0, 1, 2, 3])
        assert any(sub[0] == 1 for sub in result)


class TestGetBaseEditTargetSequence:

    def test_df_alleles(self):
        """Integration test with the real df_alleles fixture."""
        df_alleles = pd.read_csv('tests/df_alleles.txt')
        ref_seq = (
            'CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGG'
            'AATCCCTTCTGCAGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGC'
            'CGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG'
        )
        target_seq = get_base_edit_target_sequence(ref_seq, df_alleles, 0)
        assert target_seq == (
            'AATACGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGG'
            'AATCCCTTCTGCAGCCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTG'
            'CGCCGGGCCTTGCCGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGGCATGGCCCCATTCGCACGGCTCTG'
            'GAGCGGC'
        )

    def test_small_df(self):
        """Small fixture: most-common non-reference allele is the target."""
        df_alleles = pd.read_csv('tests/test_be_df.txt')
        target_seq = get_base_edit_target_sequence('AAAA', df_alleles, 0)
        assert target_seq == 'AAGA'


class TestGetUpsetPlotCounts:

    def test_basic(self):
        """Counts dict has expected keys and totals."""
        df_alleles = pd.read_csv('tests/test_be_df.txt')
        counts_dict = get_upset_plot_counts(df_alleles, [(3, 'A', 'G')], 'TEST')
        assert len(counts_dict) == 19
        assert counts_dict['total_alleles'] == 3
        assert counts_dict['total_alleles_reads'] == 100
        assert sum(counts_dict['binary_allele_counts'].values()) == 100


class TestWriteBaseEditCounts:

    def test_files_written_and_cleaned(self, tmp_path):
        """write_base_edit_counts writes exactly the five expected files."""
        df_alleles = pd.read_csv('tests/test_be_df.txt')
        counts_dict = get_upset_plot_counts(df_alleles, [(3, 'A', 'G')], 'TEST')

        def _jp(filename):
            return str(tmp_path / filename)

        write_base_edit_counts('TEST', counts_dict, [(3, 'A', 'G')], _jp)

        expected = [
            '10i.TEST.arrays.txt',
            '10i.TEST.binary_allele_counts.txt',
            '10i.TEST.category_allele_counts.txt',
            '10i.TEST.counts.txt',
            '10i.TEST.precise_allele_counts.txt',
        ]
        for filename in expected:
            assert (tmp_path / filename).exists(), f"Missing: {filename}"


class TestWriteCoreAmpliconDataFiles:

    def test_writes_expected_files_and_metadata(self, tmp_path):
        """Per-amplicon CORE exports are written and recorded in run metadata."""
        ref_name = 'amp1'
        ref_plot_name = 'Amp1_'
        ctx = _make_ctx(
            args=SimpleNamespace(
                plot_histogram_outliers=False,
                plot_window_size=20,
                allele_plot_pcts_only_for_assigned_reference=False,
                expand_allele_plots_by_quantification=True,
                conversion_nuc_from='C',
                expected_hdr_amplicon_seq='',
                base_editor_output=True,
                coding_seq='',
                crispresso1_mode=False,
            ),
            output_directory=str(tmp_path),
            _jp=lambda f: os.path.join(str(tmp_path), f),
            refs={ref_name: _ref_dict(sequence='ACGT', ref_plot_name=ref_plot_name, include_idxs=[1, 3])},
            ref_names=[ref_name],
            counts_total={ref_name: 10},
            all_base_count_vectors={
                ref_name + '_A': np.array([1, 2, 3, 4]),
                ref_name + '_C': np.array([5, 6, 7, 8]),
                ref_name + '_G': np.array([9, 10, 11, 12]),
                ref_name + '_T': np.array([13, 14, 15, 16]),
                ref_name + '_N': np.array([17, 18, 19, 20]),
                ref_name + '_-': np.array([21, 22, 23, 24]),
            },
            substitution_base_vectors={
                ref_name + '_A': np.array([31, 32]),
                ref_name + '_C': np.array([33, 34]),
                ref_name + '_G': np.array([35, 36]),
                ref_name + '_T': np.array([37, 38]),
                ref_name + '_N': np.array([39, 40]),
            },
            all_substitution_base_vectors={
                ref_name + '_A': np.array([41, 42, 43, 44]),
                ref_name + '_C': np.array([45, 46, 47, 48]),
                ref_name + '_G': np.array([49, 50, 51, 52]),
                ref_name + '_T': np.array([53, 54, 55, 56]),
                ref_name + '_N': np.array([57, 58, 59, 60]),
            },
        )
        crispresso2_info = {'results': {'refs': {ref_name: {}}}}

        write_core_amplicon_data_files(ctx, crispresso2_info)

        expected = {
            'quant_window_nuc_freq_filename': ref_plot_name + 'Quantification_window_nucleotide_frequency_table.txt',
            'quant_window_nuc_pct_filename': ref_plot_name + 'Quantification_window_nucleotide_percentage_table.txt',
            'nuc_freq_filename': ref_plot_name + 'Nucleotide_frequency_table.txt',
            'nuc_pct_filename': ref_plot_name + 'Nucleotide_percentage_table.txt',
            'quant_window_sub_freq_filename': ref_plot_name + 'Quantification_window_substitution_frequency_table.txt',
            'sub_freq_table_filename': ref_plot_name + 'Substitution_frequency_table.txt',
        }
        info_ref = crispresso2_info['results']['refs'][ref_name]
        for key, filename in expected.items():
            assert info_ref[key] == filename
            assert (tmp_path / filename).exists(), f"Missing: {filename}"

        quant_window_nuc = pd.read_csv(tmp_path / expected['quant_window_nuc_freq_filename'], sep='\t', index_col=0)
        assert list(quant_window_nuc.columns) == ['C', 'T']
        assert quant_window_nuc.loc['A'].tolist() == [2, 4]
        assert quant_window_nuc.loc['-'].tolist() == [22, 24]


# =============================================================================
# Batch / Aggregate prep functions
# =============================================================================


def _make_batch_ctx(tmp_path, **overrides):
    """Build a minimal BatchPlotContext for testing."""
    from CRISPResso2.plots.plot_context import BatchPlotContext

    seq = 'ACGTACGTAC'
    seq_cols = list(seq)
    n_rows_nuc = 4  # 2 batches × 2 nucleotides
    n_rows_mod = 2  # 2 batches
    n_rows_mod_freq = 6  # 2 batches × 3 mod types

    nuc_data = {'Batch': ['s1', 's1', 's2', 's2'], 'Nucleotide': ['A', 'C', 'A', 'C']}
    for i, base in enumerate(seq_cols):
        nuc_data[i] = [0.25] * n_rows_nuc
    nuc_pct_df = pd.DataFrame(nuc_data)
    nuc_pct_df.columns = ['Batch', 'Nucleotide'] + seq_cols

    mod_pct_data = {'Batch': ['s1', 's2'], 'Modification': ['Insertions', 'Insertions']}
    for i, base in enumerate(seq_cols):
        mod_pct_data[i] = [0.05] * n_rows_mod
    mod_pct_df = pd.DataFrame(mod_pct_data)
    mod_pct_df.columns = ['Batch', 'Modification'] + seq_cols

    mod_freq_data = {'Batch': ['s1', 's1', 's1', 's2', 's2', 's2'],
                     'Modification': ['Insertions', 'Deletions', 'Substitutions'] * 2}
    for i, base in enumerate(seq_cols):
        mod_freq_data[i] = [5] * n_rows_mod_freq
    mod_freq_df = pd.DataFrame(mod_freq_data)
    mod_freq_df.columns = ['Batch', 'Modification'] + seq_cols

    defaults = dict(
        args=SimpleNamespace(
            conversion_nuc_from='C',
            conversion_nuc_to='T',
            base_editor_output=False,
            use_matplotlib=True,
        ),
        run_data={'results': {'general_plots': {}}},
        output_directory=str(tmp_path),
        save_png=False,
        _jp=lambda f: os.path.join(str(tmp_path), f),
        custom_config={'colors': {}},
        amplicon_names=['Amp1'],
        nucleotide_frequency_summary_dfs={'Amp1': nuc_pct_df.copy()},
        nucleotide_percentage_summary_dfs={'Amp1': nuc_pct_df},
        modification_frequency_summary_dfs={'Amp1': mod_freq_df},
        modification_percentage_summary_dfs={'Amp1': mod_pct_df},
        consensus_guides={'Amp1': ['ACGT']},
        consensus_include_idxs={'Amp1': np.array([3, 4, 5, 6])},
        consensus_sgRNA_intervals={'Amp1': [(3, 7)]},
        consensus_sgRNA_plot_idxs={'Amp1': [np.array([2, 3, 4, 5, 6, 7])]},
        guides_all_same={'Amp1': True},
        amplicon_name='Amp1',
    )
    defaults.update(overrides)
    return BatchPlotContext(**defaults)


class TestComputeSubSgRNAIntervals:
    """Tests for _compute_sub_sgRNA_intervals."""

    def test_basic_mapping(self):
        from CRISPResso2.plots.data_prep import _compute_sub_sgRNA_intervals
        sgRNA_plot_idxs = np.array([2, 3, 4, 5, 6, 7])
        consensus_sgRNA_intervals = [(3, 7)]
        include_idxs = np.array([3, 4, 5, 6])
        sub_intervals, sub_include = _compute_sub_sgRNA_intervals(
            sgRNA_plot_idxs, consensus_sgRNA_intervals, include_idxs,
        )
        assert len(sub_intervals) == 1
        assert sub_intervals[0][0] >= 0
        assert sub_include[0] == 3 - 2  # offset by plot_idxs[0]


class TestPrepBatchNucQuiltAroundSgRNA:
    """Tests for prep_batch_nuc_quilt_around_sgRNA."""

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_nuc_quilt_around_sgRNA
        ctx = _make_batch_ctx(tmp_path)
        ctx.sgRNA_ind = 0
        result = prep_batch_nuc_quilt_around_sgRNA(ctx)
        assert 'nuc_pct_df' in result
        assert 'mod_pct_df' in result
        assert 'fig_filename_root' in result
        assert 'sgRNA_intervals' in result

    def test_slices_dataframe(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_nuc_quilt_around_sgRNA
        ctx = _make_batch_ctx(tmp_path)
        ctx.sgRNA_ind = 0
        result = prep_batch_nuc_quilt_around_sgRNA(ctx)
        # sliced DataFrame should have fewer columns than original
        orig_cols = ctx.nucleotide_percentage_summary_dfs['Amp1'].shape[1]
        result_cols = result['nuc_pct_df'].shape[1]
        assert result_cols <= orig_cols


class TestPrepBatchNucQuilt:
    """Tests for prep_batch_nuc_quilt."""

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_nuc_quilt
        ctx = _make_batch_ctx(tmp_path)
        result = prep_batch_nuc_quilt(ctx)
        assert 'nuc_pct_df' in result
        assert 'fig_filename_root' in result

    def test_includes_sgRNA_data_when_guides_same(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_nuc_quilt
        ctx = _make_batch_ctx(tmp_path)
        result = prep_batch_nuc_quilt(ctx)
        assert 'sgRNA_intervals' in result
        assert 'sgRNA_sequences' in result

    def test_no_sgRNA_when_guides_differ(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_nuc_quilt
        ctx = _make_batch_ctx(tmp_path, guides_all_same={'Amp1': False})
        result = prep_batch_nuc_quilt(ctx)
        assert 'sgRNA_intervals' not in result


class TestPrepBatchConversionMap:
    """Tests for prep_batch_conversion_map."""

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_conversion_map
        ctx = _make_batch_ctx(tmp_path)
        result = prep_batch_conversion_map(ctx)
        assert 'nuc_pct_df' in result
        assert 'conversion_nuc_from' in result
        assert 'conversion_nuc_to' in result


class TestPrepBatchAlleleModificationHeatmap:
    """Tests for prep_batch_allele_modification_heatmap."""

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_allele_modification_heatmap
        ctx = _make_batch_ctx(tmp_path)
        ctx.mod_type = 'Insertions'
        result = prep_batch_allele_modification_heatmap(ctx)
        assert 'sample_values' in result
        assert 'sample_sgRNA_intervals' in result
        assert 'plot_path' in result
        assert 'div_id' in result

    def test_filters_by_mod_type(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_allele_modification_heatmap
        ctx = _make_batch_ctx(tmp_path)
        ctx.mod_type = 'Insertions'
        result = prep_batch_allele_modification_heatmap(ctx)
        # Should only have rows for 'Insertions' (2 batches)
        assert result['sample_values'].shape[0] == 2


class TestPrepBatchAlleleModificationLine:
    """Tests for prep_batch_allele_modification_line."""

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_batch_allele_modification_line
        ctx = _make_batch_ctx(tmp_path)
        ctx.mod_type = 'Deletions'
        result = prep_batch_allele_modification_line(ctx)
        assert 'sample_values' in result
        assert 'plot_path' in result
        assert 'div_id' in result


# =============================================================================
# Multi-mode summary plot prep functions (Pooled, WGS, Aggregate)
# =============================================================================


class TestPrepReadsTotal:
    """Tests for prep_reads_total."""

    def _make_summary_ctx(self, tmp_path):
        """Build a minimal context with df_summary_quantification."""
        from CRISPResso2.plots.plot_context import PooledPlotContext
        df = pd.DataFrame({
            'Name': ['region1', 'region2'],
            'Reads_total': [1000, 500],
            'Reads_aligned': [900, 400],
            'Unmodified': [800, 350],
            'Modified': [100, 50],
        })
        return PooledPlotContext(
            args=SimpleNamespace(min_reads_to_use_region=10),
            run_data={},
            output_directory=str(tmp_path),
            save_png=True,
            _jp=lambda f: os.path.join(str(tmp_path), f),
            custom_config={},
            df_summary_quantification=df,
        )

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_reads_total
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert 'df_summary_quantification' in result
        assert 'fig_filename_root' in result
        assert 'save_png' in result
        assert 'cutoff' in result

    def test_filename_uses_prefix(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_reads_total
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert 'CRISPRessoPooled_reads_summary' in result['fig_filename_root']

    def test_different_prefix(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_reads_total
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_reads_total(ctx, prefix='CRISPRessoWGS')
        assert 'CRISPRessoWGS_reads_summary' in result['fig_filename_root']

    def test_cutoff_from_args(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_reads_total
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert result['cutoff'] == 10

    def test_save_png_passthrough(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_reads_total
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert result['save_png'] is True

    def test_dataframe_is_a_copy(self, tmp_path):
        """The returned df must be a copy so sanitization doesn't
        mutate ``ctx.df_summary_quantification``."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert result['df_summary_quantification'] is not ctx.df_summary_quantification


class TestSanitizeSummaryNames:
    """Tests for the Name column sanitization in prep_reads_total /
    prep_unmod_mod_pcts."""

    def _make_summary_ctx(self, tmp_path, df):
        from CRISPResso2.plots.plot_context import PooledPlotContext
        return PooledPlotContext(
            args=SimpleNamespace(min_reads_to_use_region=10),
            run_data={},
            output_directory=str(tmp_path),
            save_png=True,
            _jp=lambda f: os.path.join(str(tmp_path), f),
            custom_config={},
            df_summary_quantification=df,
        )

    def test_ascii_names_unchanged(self, tmp_path):
        """ASCII-only names should pass through slugify unchanged
        (or nearly so — slugify is a no-op for clean ASCII)."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        df = pd.DataFrame({
            'Name': ['HEK3', 'sample1'],
            'Reads_total': [100, 200],
        })
        ctx = self._make_summary_ctx(tmp_path, df)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert list(result['df_summary_quantification']['Name']) == ['HEK3', 'sample1']

    def test_emoji_stripped(self, tmp_path):
        """Emoji must be stripped — this is the core reason the
        sanitization exists."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        df = pd.DataFrame({
            'Name': ['FAN /+-C%_:@#%a😀sdf    -'],
            'Reads_total': [100],
        })
        ctx = self._make_summary_ctx(tmp_path, df)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        sanitized = list(result['df_summary_quantification']['Name'])[0]
        assert '😀' not in sanitized
        assert sanitized == 'FAN_+-C%_@#%asdf_-'

    def test_source_dataframe_not_mutated(self, tmp_path):
        """``ctx.df_summary_quantification`` must be unchanged by
        the prep call — sanitization only affects the returned copy.
        CSVs / crispresso2_info / log lines still see raw names."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        original = 'FAN /+-C%_:@#%a😀sdf    -'
        df = pd.DataFrame({'Name': [original], 'Reads_total': [100]})
        ctx = self._make_summary_ctx(tmp_path, df)
        prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert list(ctx.df_summary_quantification['Name']) == [original]

    def test_non_string_entries_coerced(self, tmp_path):
        """Integer / None / NaN entries in Name are coerced via
        str() before slugify — must not crash."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        df = pd.DataFrame({
            'Name': [1, 'HEK3', None, np.nan],
            'Reads_total': [100, 200, 300, 400],
        })
        ctx = self._make_summary_ctx(tmp_path, df)
        # Must not raise
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        sanitized = list(result['df_summary_quantification']['Name'])
        assert sanitized[0] == '1'
        assert sanitized[1] == 'HEK3'
        # None and NaN become the string 'None' / 'nan' via str() coercion.
        # The exact text doesn't matter; the point is no crash.
        assert len(sanitized) == 4

    def test_missing_name_column_ok(self, tmp_path):
        """If the df has no 'Name' column (unusual), sanitization
        is a no-op — must not crash."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        df = pd.DataFrame({'Reads_total': [100]})
        ctx = self._make_summary_ctx(tmp_path, df)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        # Copy should still be returned; just without a Name column.
        assert 'Name' not in result['df_summary_quantification'].columns

    def test_empty_dataframe(self, tmp_path):
        """Empty DataFrame with a Name column — sanitization is a
        no-op, returns empty df."""
        from CRISPResso2.plots.data_prep import prep_reads_total
        df = pd.DataFrame({'Name': [], 'Reads_total': []})
        ctx = self._make_summary_ctx(tmp_path, df)
        result = prep_reads_total(ctx, prefix='CRISPRessoPooled')
        assert len(result['df_summary_quantification']) == 0

    def test_applies_to_prep_unmod_mod_pcts_too(self, tmp_path):
        """The sanitization must also apply via prep_unmod_mod_pcts,
        not just prep_reads_total."""
        from CRISPResso2.plots.data_prep import prep_unmod_mod_pcts
        df = pd.DataFrame({
            'Name': ['clean', 'emoji😀here'],
            'Reads_total': [100, 200],
        })
        ctx = self._make_summary_ctx(tmp_path, df)
        result = prep_unmod_mod_pcts(ctx, prefix='CRISPRessoPooled')
        sanitized = list(result['df_summary_quantification']['Name'])
        assert '😀' not in sanitized[1]


class TestPrepUnmodModPcts:
    """Tests for prep_unmod_mod_pcts."""

    def _make_summary_ctx(self, tmp_path):
        """Build a minimal context with df_summary_quantification."""
        from CRISPResso2.plots.plot_context import PooledPlotContext
        df = pd.DataFrame({
            'Name': ['region1', 'region2'],
            'Reads_total': [1000, 500],
            'Reads_aligned': [900, 400],
            'Unmodified': [800, 350],
            'Modified': [100, 50],
        })
        return PooledPlotContext(
            args=SimpleNamespace(min_reads_to_use_region=10),
            run_data={},
            output_directory=str(tmp_path),
            save_png=False,
            _jp=lambda f: os.path.join(str(tmp_path), f),
            custom_config={},
            df_summary_quantification=df,
        )

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_unmod_mod_pcts
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_unmod_mod_pcts(ctx, prefix='CRISPRessoPooled')
        assert 'df_summary_quantification' in result
        assert 'fig_filename_root' in result
        assert 'save_png' in result
        assert 'cutoff' in result

    def test_filename_uses_prefix(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_unmod_mod_pcts
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_unmod_mod_pcts(ctx, prefix='CRISPRessoPooled')
        assert 'CRISPRessoPooled_modification_summary' in result['fig_filename_root']

    def test_save_png_passthrough(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_unmod_mod_pcts
        ctx = self._make_summary_ctx(tmp_path)
        result = prep_unmod_mod_pcts(ctx, prefix='CRISPRessoPooled')
        assert result['save_png'] is False

    def test_works_with_wgs_context(self, tmp_path):
        """Verify the function works with WGSPlotContext too."""
        from CRISPResso2.plots.plot_context import WGSPlotContext
        from CRISPResso2.plots.data_prep import prep_unmod_mod_pcts
        df = pd.DataFrame({
            'Name': ['region1'],
            'Reads_total': [1000],
            'Reads_aligned': [900],
            'Unmodified': [800],
            'Modified': [100],
        })
        ctx = WGSPlotContext(
            args=SimpleNamespace(min_reads_to_use_region=50),
            run_data={},
            output_directory=str(tmp_path),
            save_png=True,
            _jp=lambda f: os.path.join(str(tmp_path), f),
            custom_config={},
            df_summary_quantification=df,
        )
        result = prep_unmod_mod_pcts(ctx, prefix='CRISPRessoWGS')
        assert 'CRISPRessoWGS_modification_summary' in result['fig_filename_root']
        assert result['cutoff'] == 50


# =============================================================================
# Compare prep functions
# =============================================================================


def _make_compare_ctx(tmp_path, **overrides):
    """Build a minimal ComparePlotContext for testing."""
    from CRISPResso2.plots.plot_context import ComparePlotContext
    defaults = dict(
        args=SimpleNamespace(
            reported_qvalue_cutoff=0.05,
            min_frequency_alleles_around_cut_to_plot=0.05,
            max_rows_alleles_around_cut_to_plot=50,
            offset_around_cut_to_plot=20,
            crispresso_output_folder_1='/tmp/s1',
            crispresso_output_folder_2='/tmp/s2',
        ),
        run_data={'results': {'general_plots': {}}},
        output_directory=str(tmp_path),
        save_png=False,
        _jp=lambda f: os.path.join(str(tmp_path), f),
        custom_config={},
        amplicon_names=['Amp1'],
        sample_1_name='WT',
        sample_2_name='KO',
        run_info_1={'running_info': {'report_filename': 'report.html'}, 'results': {'refs': {'Amp1': {'sgRNA_cut_points': [10], 'sgRNA_intervals': [(5, 15)], 'include_idxs': np.array([8, 9, 10, 11, 12])}}}},
        run_info_2={'running_info': {'report_filename': 'report.html'}, 'results': {'refs': {'Amp1': {'sgRNA_cut_points': [10], 'sgRNA_intervals': [(5, 15)], 'include_idxs': np.array([8, 9, 10, 11, 12])}}}},
        amplicon_info_1={'Amp1': {'Reads_aligned': '1000', 'Unmodified': '800', 'Modified': '200'}},
        amplicon_info_2={'Amp1': {'Reads_aligned': '1000', 'Unmodified': '900', 'Modified': '100'}},
        amplicon_name='Amp1',
    )
    defaults.update(overrides)
    return ComparePlotContext(**defaults)


class TestPrepCompareEditingBarchart:
    """Tests for prep_compare_editing_barchart."""

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_editing_barchart
        ctx = _make_compare_ctx(tmp_path)
        result = prep_compare_editing_barchart(ctx)
        assert 'n_total_1' in result
        assert 'n_total_2' in result
        assert 'sample_1_name' in result
        assert 'plot_titles' in result
        assert 'plot_path' in result

    def test_reads_from_amplicon_info(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_editing_barchart
        ctx = _make_compare_ctx(tmp_path)
        result = prep_compare_editing_barchart(ctx)
        assert result['n_total_1'] == 1000.0
        assert result['n_modified_1'] == 200.0
        assert result['n_unmodified_2'] == 900.0

    def test_sample_names_passthrough(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_editing_barchart
        ctx = _make_compare_ctx(tmp_path)
        result = prep_compare_editing_barchart(ctx)
        assert result['sample_1_name'] == 'WT'
        assert result['sample_2_name'] == 'KO'


class TestPrepCompareModificationPositions:
    """Tests for prep_compare_modification_positions."""

    def _make_ctx_with_mod_freqs(self, tmp_path):
        ctx = _make_compare_ctx(tmp_path)
        seq = 'ACGTACGTACGTACGTACGT'
        ctx.mod_freqs_1 = {
            'Insertions': [10] * len(seq),
            'Deletions': [5] * len(seq),
            'Substitutions': [3] * len(seq),
            'All_modifications': [15] * len(seq),
            'Total': [100] * len(seq),
        }
        ctx.mod_freqs_2 = {
            'Insertions': [5] * len(seq),
            'Deletions': [10] * len(seq),
            'Substitutions': [2] * len(seq),
            'All_modifications': [12] * len(seq),
            'Total': [100] * len(seq),
        }
        ctx.consensus_sequence = seq
        ctx.cut_points = [10]
        ctx.sgRNA_intervals = [(5, 15)]
        ctx.quant_windows_1 = np.array([8, 9, 10, 11, 12])
        ctx.quant_windows_2 = np.array([8, 9, 10, 11, 12])
        return ctx

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_modification_positions
        ctx = self._make_ctx_with_mod_freqs(tmp_path)
        ctx.mod_type = 'Insertions'
        result = prep_compare_modification_positions(ctx)
        assert 'plot_kwargs' in result
        assert 'mod_df' in result
        assert 'sig_count' in result
        assert 'sig_count_quant_window' in result

    def test_fisher_test_produces_pvalues(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_modification_positions
        ctx = self._make_ctx_with_mod_freqs(tmp_path)
        ctx.mod_type = 'Insertions'
        result = prep_compare_modification_positions(ctx)
        pvalues = result['plot_kwargs']['pvalues']
        assert len(pvalues) == len(ctx.consensus_sequence)
        # All positions have same counts so p-values should all be equal
        assert all(p == pvalues[0] for p in pvalues)

    def test_sig_count_is_int(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_modification_positions
        ctx = self._make_ctx_with_mod_freqs(tmp_path)
        ctx.mod_type = 'Insertions'
        result = prep_compare_modification_positions(ctx)
        assert isinstance(result['sig_count'], (int, np.integer))

    def test_mod_df_has_correct_shape(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_modification_positions
        ctx = self._make_ctx_with_mod_freqs(tmp_path)
        ctx.mod_type = 'Deletions'
        result = prep_compare_modification_positions(ctx)
        # 7 rows: sample1_mod, sample1_total, sample2_mod, sample2_total, odds_ratios, pvalues, qval_bonferroni
        assert result['mod_df'].shape[0] == 7
        # columns: Reference + one per base
        assert result['mod_df'].shape[1] == len(ctx.consensus_sequence) + 1


class TestPrepCompareAlleleTable:
    """Tests for prep_compare_allele_table."""

    def _make_ctx_with_allele_pairs(self, tmp_path):
        ctx = _make_compare_ctx(tmp_path)
        ctx.consensus_sequence = 'ACGTACGTACGTACGTACGT'
        ctx.sgRNA_intervals = [(5, 15)]
        ctx.cut_points = [10]

        df1 = pd.DataFrame({
            'Aligned_Sequence': ['ACGT', 'ACGA'],
            'Reference_Sequence': ['ACGT', 'ACGT'],
            'Unedited': [True, False],
            'n_deleted': [0, 0],
            'n_inserted': [0, 0],
            'n_mutated': [0, 1],
            '#Reads': [80, 20],
            '%Reads': [80.0, 20.0],
        })
        df2 = pd.DataFrame({
            'Aligned_Sequence': ['ACGT', 'ACGC'],
            'Reference_Sequence': ['ACGT', 'ACGT'],
            'Unedited': [True, False],
            'n_deleted': [0, 0],
            'n_inserted': [0, 0],
            'n_mutated': [0, 1],
            '#Reads': [60, 40],
            '%Reads': [60.0, 40.0],
        })
        ctx.allele_pairs = [('alleles_around_sgRNA.txt', 'alleles_around_sgRNA.txt', df1, df2)]
        return ctx

    def test_returns_expected_keys(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_allele_table
        ctx = self._make_ctx_with_allele_pairs(tmp_path)
        result = prep_compare_allele_table(ctx, 0)
        assert 'merged_df' in result
        assert 'plot_top_kwargs' in result
        assert 'plot_bottom_kwargs' in result
        assert 'is_base_edit' in result
        assert 'file_root' in result

    def test_merged_df_has_lfc(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_allele_table
        ctx = self._make_ctx_with_allele_pairs(tmp_path)
        result = prep_compare_allele_table(ctx, 0)
        assert 'each_LFC' in result['merged_df'].columns

    def test_top_bottom_sorted_differently(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_allele_table
        ctx = self._make_ctx_with_allele_pairs(tmp_path)
        result = prep_compare_allele_table(ctx, 0)
        top_lfc = result['plot_top_kwargs']['df_alleles']['each_LFC'].iloc[0]
        bottom_lfc = result['plot_bottom_kwargs']['df_alleles']['each_LFC'].iloc[0]
        # Top should be enriched in sample_1 (higher LFC), bottom in sample_2 (lower LFC)
        assert top_lfc >= bottom_lfc

    def test_is_base_edit_detection(self, tmp_path):
        from CRISPResso2.plots.data_prep import prep_compare_allele_table
        ctx = self._make_ctx_with_allele_pairs(tmp_path)
        result = prep_compare_allele_table(ctx, 0)
        assert result['is_base_edit'] is False

        # Test with base_edit in filename
        ctx.allele_pairs = [('base_edit_alleles.txt', 'base_edit_alleles.txt',
                             ctx.allele_pairs[0][2], ctx.allele_pairs[0][3])]
        result = prep_compare_allele_table(ctx, 0)
        assert result['is_base_edit'] is True
