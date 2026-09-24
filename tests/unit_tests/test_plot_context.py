"""Unit tests for PlotContext hierarchy."""

import pytest
from types import SimpleNamespace


# =============================================================================
# Helpers
# =============================================================================


def _make_base_kwargs(**overrides):
    """Build the 6 required base PlotContext kwargs."""
    defaults = dict(
        args=SimpleNamespace(),
        run_data={},
        output_directory="",
        save_png=False,
        _jp=lambda f: f,
        custom_config={},
    )
    defaults.update(overrides)
    return defaults


def _make_core_required_kwargs(**overrides):
    """Build the minimal required kwargs for CorePlotContext."""
    defaults = _make_base_kwargs()
    defaults.update(dict(
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
    ))
    defaults.update(overrides)
    return defaults


# =============================================================================
# Base PlotContext
# =============================================================================


class TestPlotContextBase:
    """Test the PlotContext base class."""

    def test_construction(self):
        from CRISPResso2.plots.plot_context import PlotContext

        ctx = PlotContext(**_make_base_kwargs(
            args=SimpleNamespace(name='test'),
            run_data={'results': {}},
            output_directory='/tmp/out',
            save_png=True,
        ))
        assert ctx.args.name == 'test'
        assert ctx.run_data == {'results': {}}
        assert ctx.output_directory == '/tmp/out'
        assert ctx.save_png is True
        assert ctx._jp is not None
        assert ctx.custom_config == {}

    def test_all_fields_required(self):
        """PlotContext has no defaults — omitting any field raises TypeError."""
        from CRISPResso2.plots.plot_context import PlotContext

        with pytest.raises(TypeError):
            PlotContext(args=SimpleNamespace())  # missing other required fields

    def test_jp_callable(self):
        from CRISPResso2.plots.plot_context import PlotContext

        jp = lambda f: '/out/' + f
        ctx = PlotContext(**_make_base_kwargs(_jp=jp))
        assert ctx._jp('file.txt') == '/out/file.txt'


# =============================================================================
# CorePlotContext
# =============================================================================


class TestCorePlotContextConstruction:
    """Test CorePlotContext can be constructed and fields are accessible."""

    def test_minimal_construction(self):
        """CorePlotContext can be constructed with required fields only."""
        from CRISPResso2.plots.plot_context import CorePlotContext

        ctx = CorePlotContext(**_make_core_required_kwargs(
            args=SimpleNamespace(name='test'),
            run_data={'results': {}},
            refs={'ref1': {}},
            ref_names=['ref1'],
            counts_total={'ref1': 100},
            counts_modified={'ref1': 50},
            counts_unmodified={'ref1': 50},
            counts_discarded={'ref1': 0},
            counts_insertion={'ref1': 10},
            counts_deletion={'ref1': 20},
            counts_substitution={'ref1': 5},
            class_counts={'Reference': 50, 'NHEJ': 50},
            N_TOTAL=100,
            all_insertion_count_vectors={'ref1': [0, 1, 2]},
            all_insertion_left_count_vectors={'ref1': [0, 0, 1]},
            all_deletion_count_vectors={'ref1': [0, 2, 0]},
            all_substitution_count_vectors={'ref1': [1, 0, 0]},
            all_indelsub_count_vectors={'ref1': [1, 2, 2]},
            all_substitution_base_vectors={'ref1': {}},
            all_base_count_vectors={'ref1': {}},
            insertion_count_vectors={'ref1': [0, 1]},
            deletion_count_vectors={'ref1': [0, 1]},
            substitution_count_vectors={'ref1': [1, 0]},
            insertion_length_vectors={'ref1': [0, 1]},
            deletion_length_vectors={'ref1': [0, 1]},
            hists_frameshift={'ref1': {}},
            hists_inframe={'ref1': {}},
            counts_modified_frameshift={'ref1': 0},
            counts_modified_non_frameshift={'ref1': 0},
            counts_non_modified_non_frameshift={'ref1': 0},
            counts_splicing_sites_modified={'ref1': 0},
        ))
        assert ctx.args.name == 'test'
        assert ctx.ref_names == ['ref1']
        assert ctx.N_TOTAL == 100

    def test_defaults(self):
        """Optional fields have correct defaults."""
        from CRISPResso2.plots.plot_context import CorePlotContext

        ctx = CorePlotContext(**_make_core_required_kwargs())
        # HDR vectors default to empty dicts
        assert ctx.ref1_all_insertion_count_vectors == {}
        assert ctx.ref1_all_deletion_count_vectors == {}
        assert ctx.ref1_all_substitution_count_vectors == {}
        assert ctx.ref1_all_indelsub_count_vectors == {}
        assert ctx.ref1_all_insertion_left_count_vectors == {}
        assert ctx.ref1_all_base_count_vectors == {}
        # New data fields default to empty
        assert ctx.class_counts_order == []
        assert ctx.homology_scores == []
        assert ctx.homology_counts == []
        assert ctx.insertion_count_vectors_noncoding == {}
        assert ctx.deletion_count_vectors_noncoding == {}
        assert ctx.substitution_count_vectors_noncoding == {}
        assert ctx.substitution_base_vectors == {}
        assert ctx.df_scaffold_insertion_sizes is None
        # Scope fields default to None
        assert ctx.ref_name is None
        assert ctx.sgRNA_ind is None
        assert ctx.coding_seq_ind is None

    def test_scope_fields_mutable(self):
        """ref_name and sgRNA_ind can be set after construction."""
        from CRISPResso2.plots.plot_context import CorePlotContext

        ctx = CorePlotContext(**_make_core_required_kwargs(
            refs={'amp1': {}, 'amp2': {}},
            ref_names=['amp1', 'amp2'],
        ))
        ctx.ref_name = 'amp1'
        ctx.sgRNA_ind = 0
        assert ctx.ref_name == 'amp1'
        assert ctx.sgRNA_ind == 0

        ctx.ref_name = 'amp2'
        ctx.sgRNA_ind = 1
        assert ctx.ref_name == 'amp2'
        assert ctx.sgRNA_ind == 1

    def test_new_fields_construction(self):
        """CorePlotContext can be constructed with the additional data fields."""
        from CRISPResso2.plots.plot_context import CorePlotContext
        import numpy as np

        ctx = CorePlotContext(**_make_core_required_kwargs(
            class_counts_order=['Ref_UNMODIFIED', 'Ref_MODIFIED'],
            homology_scores=[95.0, 80.0, 60.0],
            homology_counts=[100, 50, 10],
            insertion_count_vectors_noncoding={'ref1': np.array([0, 1])},
            deletion_count_vectors_noncoding={'ref1': np.array([1, 0])},
            substitution_count_vectors_noncoding={'ref1': np.array([0, 0])},
            substitution_base_vectors={'ref1_A': [1, 2], 'ref1_C': [0, 0]},
            df_scaffold_insertion_sizes='fake_df',
        ))
        assert ctx.class_counts_order == ['Ref_UNMODIFIED', 'Ref_MODIFIED']
        assert ctx.homology_scores == [95.0, 80.0, 60.0]
        assert ctx.homology_counts == [100, 50, 10]
        assert len(ctx.insertion_count_vectors_noncoding) == 1
        assert ctx.substitution_base_vectors == {'ref1_A': [1, 2], 'ref1_C': [0, 0]}
        assert ctx.df_scaffold_insertion_sizes == 'fake_df'

    def test_references_not_copied(self):
        """CorePlotContext holds references, not copies, of data structures."""
        from CRISPResso2.plots.plot_context import CorePlotContext

        counts = {'ref1': 100}
        ctx = CorePlotContext(**_make_core_required_kwargs(counts_total=counts))
        # Mutate original — CorePlotContext sees it (zero-copy)
        counts['ref1'] = 200
        assert ctx.counts_total['ref1'] == 200

    def test_base_fields_accessible(self):
        """Base class fields (args, _jp, etc.) are accessible on CorePlotContext."""
        from CRISPResso2.plots.plot_context import CorePlotContext

        jp = lambda f: '/out/' + f
        ctx = CorePlotContext(**_make_core_required_kwargs(
            output_directory='/out',
            save_png=True,
            _jp=jp,
            custom_config={'colors': {}},
        ))
        assert ctx.output_directory == '/out'
        assert ctx.save_png is True
        assert ctx._jp('file.txt') == '/out/file.txt'
        assert ctx.custom_config == {'colors': {}}


# =============================================================================
# Inheritance contract
# =============================================================================


class TestInheritance:
    """Test that all context classes inherit from PlotContext."""

    def test_core_isinstance(self):
        from CRISPResso2.plots.plot_context import PlotContext, CorePlotContext

        ctx = CorePlotContext(**_make_core_required_kwargs())
        assert isinstance(ctx, PlotContext)

    def test_batch_isinstance(self):
        from CRISPResso2.plots.plot_context import PlotContext, BatchPlotContext

        ctx = BatchPlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1'],
            nucleotide_frequency_summary_dfs={},
            nucleotide_percentage_summary_dfs={},
            modification_frequency_summary_dfs={},
            modification_percentage_summary_dfs={},
            consensus_guides={},
            consensus_include_idxs={},
            consensus_sgRNA_intervals={},
            consensus_sgRNA_plot_idxs={},
            guides_all_same={},
        )
        assert isinstance(ctx, PlotContext)

    def test_aggregate_isinstance(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import PlotContext, AggregatePlotContext

        ctx = AggregatePlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1'],
            nucleotide_frequency_summary_dfs={},
            nucleotide_percentage_summary_dfs={},
            modification_frequency_summary_dfs={},
            modification_percentage_summary_dfs={},
            consensus_guides={},
            consensus_include_idxs={},
            consensus_sgRNA_intervals={},
            consensus_sgRNA_plot_idxs={},
            guides_all_same={},
            df_summary_quantification=pd.DataFrame(),
            sample_count={},
        )
        assert isinstance(ctx, PlotContext)

    def test_pooled_isinstance(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import PlotContext, PooledPlotContext

        ctx = PooledPlotContext(**_make_base_kwargs(), df_summary_quantification=pd.DataFrame())
        assert isinstance(ctx, PlotContext)

    def test_wgs_isinstance(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import PlotContext, WGSPlotContext

        ctx = WGSPlotContext(**_make_base_kwargs(), df_summary_quantification=pd.DataFrame())
        assert isinstance(ctx, PlotContext)

    def test_compare_isinstance(self):
        from CRISPResso2.plots.plot_context import PlotContext, ComparePlotContext

        ctx = ComparePlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1'],
            sample_1_name='Sample 1',
            sample_2_name='Sample 2',
            run_info_1={},
            run_info_2={},
            amplicon_info_1={},
            amplicon_info_2={},
        )
        assert isinstance(ctx, PlotContext)


# =============================================================================
# BatchPlotContext
# =============================================================================


class TestBatchPlotContext:
    """Test BatchPlotContext construction and scope fields."""

    def test_construction(self):
        from CRISPResso2.plots.plot_context import BatchPlotContext

        ctx = BatchPlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1', 'amp2'],
            nucleotide_frequency_summary_dfs={'amp1': 'df1'},
            nucleotide_percentage_summary_dfs={'amp1': 'df1'},
            modification_frequency_summary_dfs={'amp1': 'df1'},
            modification_percentage_summary_dfs={'amp1': 'df1'},
            consensus_guides={'amp1': ['ATCG']},
            consensus_include_idxs={'amp1': [0, 1, 2]},
            consensus_sgRNA_intervals={'amp1': [(5, 25)]},
            consensus_sgRNA_plot_idxs={'amp1': [[0, 1, 2]]},
            guides_all_same={'amp1': True, 'amp2': False},
        )
        assert ctx.amplicon_names == ['amp1', 'amp2']
        assert ctx.guides_all_same['amp1'] is True

    def test_scope_fields_default_none(self):
        from CRISPResso2.plots.plot_context import BatchPlotContext

        ctx = BatchPlotContext(
            **_make_base_kwargs(),
            amplicon_names=[],
            nucleotide_frequency_summary_dfs={},
            nucleotide_percentage_summary_dfs={},
            modification_frequency_summary_dfs={},
            modification_percentage_summary_dfs={},
            consensus_guides={},
            consensus_include_idxs={},
            consensus_sgRNA_intervals={},
            consensus_sgRNA_plot_idxs={},
            guides_all_same={},
        )
        assert ctx.amplicon_name is None
        assert ctx.sgRNA_ind is None

    def test_scope_fields_mutable(self):
        from CRISPResso2.plots.plot_context import BatchPlotContext

        ctx = BatchPlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1'],
            nucleotide_frequency_summary_dfs={},
            nucleotide_percentage_summary_dfs={},
            modification_frequency_summary_dfs={},
            modification_percentage_summary_dfs={},
            consensus_guides={},
            consensus_include_idxs={},
            consensus_sgRNA_intervals={},
            consensus_sgRNA_plot_idxs={},
            guides_all_same={},
        )
        ctx.amplicon_name = 'amp1'
        ctx.sgRNA_ind = 0
        assert ctx.amplicon_name == 'amp1'
        assert ctx.sgRNA_ind == 0


# =============================================================================
# AggregatePlotContext
# =============================================================================


class TestAggregatePlotContext:
    """Test AggregatePlotContext construction and extra fields."""

    def test_construction(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import AggregatePlotContext

        df = pd.DataFrame({'Name': ['s1'], 'Reads_total': [100]})
        ctx = AggregatePlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1'],
            nucleotide_frequency_summary_dfs={},
            nucleotide_percentage_summary_dfs={},
            modification_frequency_summary_dfs={},
            modification_percentage_summary_dfs={},
            consensus_guides={},
            consensus_include_idxs={},
            consensus_sgRNA_intervals={},
            consensus_sgRNA_plot_idxs={},
            guides_all_same={},
            df_summary_quantification=df,
            sample_count={'amp1': 5},
        )
        assert ctx.sample_count['amp1'] == 5
        assert len(ctx.df_summary_quantification) == 1

    def test_scope_fields_default_none(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import AggregatePlotContext

        ctx = AggregatePlotContext(
            **_make_base_kwargs(),
            amplicon_names=[],
            nucleotide_frequency_summary_dfs={},
            nucleotide_percentage_summary_dfs={},
            modification_frequency_summary_dfs={},
            modification_percentage_summary_dfs={},
            consensus_guides={},
            consensus_include_idxs={},
            consensus_sgRNA_intervals={},
            consensus_sgRNA_plot_idxs={},
            guides_all_same={},
            df_summary_quantification=pd.DataFrame(),
            sample_count={},
        )
        assert ctx.amplicon_name is None
        assert ctx.sgRNA_ind is None


# =============================================================================
# PooledPlotContext
# =============================================================================


class TestPooledPlotContext:
    """Test PooledPlotContext construction."""

    def test_construction(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import PooledPlotContext

        df = pd.DataFrame({'Name': ['r1'], 'Reads_total': [500]})
        ctx = PooledPlotContext(**_make_base_kwargs(), df_summary_quantification=df)
        assert len(ctx.df_summary_quantification) == 1
        assert ctx.args is not None


# =============================================================================
# WGSPlotContext
# =============================================================================


class TestWGSPlotContext:
    """Test WGSPlotContext construction."""

    def test_construction(self):
        import pandas as pd
        from CRISPResso2.plots.plot_context import WGSPlotContext

        df = pd.DataFrame({'Name': ['r1'], 'Reads_total': [500]})
        ctx = WGSPlotContext(**_make_base_kwargs(), df_summary_quantification=df)
        assert len(ctx.df_summary_quantification) == 1
        assert ctx.args is not None


# =============================================================================
# ComparePlotContext
# =============================================================================


class TestComparePlotContext:
    """Test ComparePlotContext construction and scope fields."""

    def test_construction(self):
        from CRISPResso2.plots.plot_context import ComparePlotContext

        ctx = ComparePlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1', 'amp2'],
            sample_1_name='WT',
            sample_2_name='KO',
            run_info_1={'results': {}},
            run_info_2={'results': {}},
            amplicon_info_1={'amp1': {}},
            amplicon_info_2={'amp1': {}},
        )
        assert ctx.sample_1_name == 'WT'
        assert ctx.sample_2_name == 'KO'
        assert ctx.amplicon_names == ['amp1', 'amp2']

    def test_scope_field_default_none(self):
        from CRISPResso2.plots.plot_context import ComparePlotContext

        ctx = ComparePlotContext(
            **_make_base_kwargs(),
            amplicon_names=[],
            sample_1_name='S1',
            sample_2_name='S2',
            run_info_1={},
            run_info_2={},
            amplicon_info_1={},
            amplicon_info_2={},
        )
        assert ctx.amplicon_name is None

    def test_scope_field_mutable(self):
        from CRISPResso2.plots.plot_context import ComparePlotContext

        ctx = ComparePlotContext(
            **_make_base_kwargs(),
            amplicon_names=['amp1'],
            sample_1_name='S1',
            sample_2_name='S2',
            run_info_1={},
            run_info_2={},
            amplicon_info_1={},
            amplicon_info_2={},
        )
        ctx.amplicon_name = 'amp1'
        assert ctx.amplicon_name == 'amp1'


# =============================================================================
# Import
# =============================================================================


class TestImports:
    """Test that all context classes are importable."""

    def test_import_all(self):
        from CRISPResso2.plots.plot_context import (
            PlotContext,
            CorePlotContext,
            BatchPlotContext,
            AggregatePlotContext,
            PooledPlotContext,
            WGSPlotContext,
            ComparePlotContext,
        )
        assert PlotContext is not None
        assert CorePlotContext is not None
        assert BatchPlotContext is not None
        assert AggregatePlotContext is not None
        assert PooledPlotContext is not None
        assert WGSPlotContext is not None
        assert ComparePlotContext is not None
