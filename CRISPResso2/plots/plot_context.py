"""PlotContext hierarchy: typed data views for CRISPResso2 plot plugins.

This module defines the base :class:`PlotContext` and mode-specific
subclasses that CRISPRessoPro and plugins consume to generate custom
plots.  Each CORE module constructs its corresponding context once per
run, after analysis completes but before the plotting section begins.

Hierarchy::

    PlotContext (base — 6 required fields, no defaults)
    ├── CorePlotContext      (single-sample CRISPR analysis)
    ├── BatchPlotContext     (cross-batch summary data)
    ├── AggregatePlotContext (cross-folder summary data)
    ├── PooledPlotContext    (region summary)
    ├── WGSPlotContext       (region summary)
    └── ComparePlotContext   (paired-sample comparison)

"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import pandas as pd


# =============================================================================
# Base class — shared configuration / plumbing
# =============================================================================


@dataclass
class PlotContext:
    """Base context shared by every CRISPResso2 mode.

    Contains only configuration and plumbing fields that are common to
    all modes.  All fields are **required** (no defaults) so that
    dataclass inheritance works on Python < 3.10 without ``kw_only``.

    .. warning::

        Contexts should be treated as read-only (except for scope
        fields).  Mutating the underlying dicts, arrays, or DataFrames
        will corrupt the analysis results and the output report.
    """

    args: argparse.Namespace            # All CLI arguments
    run_data: dict                      # crispresso2_info dict
    output_directory: str               # Output path
    save_png: bool                      # Whether to save PNG alongside HTML/SVG
    _jp: Callable[[str], str]           # Joins filename with output directory
    custom_config: dict                 # Color/style/report configuration


# =============================================================================
# Core (single-sample) context
# =============================================================================


@dataclass
class CorePlotContext(PlotContext):
    """View into CRISPRessoCORE's computed data for a single-sample run.

    Constructed once per run by CORE after the analysis phase completes.
    Wraps references to CORE's local variables — no data is copied.

    .. warning::

        All fields except ``ref_name``, ``sgRNA_ind``, and
        ``coding_seq_ind`` should be treated as read-only.  Mutating the
        underlying dicts, arrays, or DataFrames will corrupt the analysis
        results and the output report.
    """

    # === Run-level data (always available) ===

    refs: dict                      # Per-amplicon reference data
    ref_names: list[str]            # Ordered list of amplicon/reference names

    # Alignment counts (ref_name -> int)
    counts_total: dict[str, int]
    counts_modified: dict[str, int]
    counts_unmodified: dict[str, int]
    counts_discarded: dict[str, int]
    counts_insertion: dict[str, int]
    counts_deletion: dict[str, int]
    counts_substitution: dict[str, int]

    # Class counts (class_name -> count)
    class_counts: dict[str, int]
    N_TOTAL: int

    # Allele data
    df_alleles: pd.DataFrame

    # Per-position count vectors (ref_name -> numpy array)
    all_insertion_count_vectors: dict[str, np.ndarray]
    all_insertion_left_count_vectors: dict[str, np.ndarray]
    all_deletion_count_vectors: dict[str, np.ndarray]
    all_substitution_count_vectors: dict[str, np.ndarray]
    all_indelsub_count_vectors: dict[str, np.ndarray]
    all_substitution_base_vectors: dict[str, np.ndarray]
    all_base_count_vectors: dict[str, np.ndarray]

    # Quantification window vectors
    insertion_count_vectors: dict[str, np.ndarray]
    deletion_count_vectors: dict[str, np.ndarray]
    substitution_count_vectors: dict[str, np.ndarray]
    insertion_length_vectors: dict[str, np.ndarray]
    deletion_length_vectors: dict[str, np.ndarray]

    # Histograms (ref_name -> Counter)
    hists_frameshift: dict
    hists_inframe: dict

    # Frameshift / coding counts (ref_name -> int)
    counts_modified_frameshift: dict[str, int]
    counts_modified_non_frameshift: dict[str, int]
    counts_non_modified_non_frameshift: dict[str, int]
    counts_splicing_sites_modified: dict[str, int]

    # === Fields with defaults (must follow all required fields) ===

    # HDR / ref1-aligned vectors (populated when expected_hdr_amplicon_seq is set)
    ref1_all_insertion_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_deletion_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_substitution_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_indelsub_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_insertion_left_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)
    ref1_all_base_count_vectors: dict[str, np.ndarray] = field(default_factory=dict)

    # === Additional data fields (populated by CORE when available) ===

    # Allele classification ordering (for plot_1b/1c)
    class_counts_order: list[str] = field(default_factory=list)

    # Allele homology data (for plot_1e)
    homology_scores: list[float] = field(default_factory=list)
    homology_counts: list[int] = field(default_factory=list)

    # Noncoding modification vectors (for plot_7; ref_name → numpy array)
    insertion_count_vectors_noncoding: dict[str, np.ndarray] = field(default_factory=dict)
    deletion_count_vectors_noncoding: dict[str, np.ndarray] = field(default_factory=dict)
    substitution_count_vectors_noncoding: dict[str, np.ndarray] = field(default_factory=dict)

    # Quantification-window substitution base vectors (for plot_10c)
    # Same shape as all_substitution_base_vectors but window-only
    substitution_base_vectors: dict[str, np.ndarray] = field(default_factory=dict)

    # Scaffold insertion sizes DataFrame (for plot_11c)
    df_scaffold_insertion_sizes: Optional[pd.DataFrame] = None

    # Alternate allele counts (for plots 10b/10c; populated when base_editor_output)
    # ref_name → {ref_nuc: {obs_nuc: count}}
    alt_nuc_counts: dict[str, dict] = field(default_factory=dict)       # quantification window
    alt_nuc_counts_all: dict[str, dict] = field(default_factory=dict)   # full amplicon

    # === Scope fields (set by Pro during iteration) ===

    ref_name: Optional[str] = None
    sgRNA_ind: Optional[int] = None
    coding_seq_ind: Optional[int] = None


# =============================================================================
# Batch context
# =============================================================================


@dataclass
class BatchPlotContext(PlotContext):
    """View into CRISPRessoBatchCORE's aggregated data for cross-batch plots.

    Constructed once per run after per-amplicon summary DataFrames are
    built.  All per-amplicon data is stored in dicts keyed by amplicon
    name; the ``amplicon_name`` scope field selects the current slice
    during iteration.
    """

    # === Run-level data ===

    amplicon_names: list[str]                                   # All amplicon names in this batch

    # Per-amplicon summary DataFrames (amplicon_name -> DataFrame)
    nucleotide_frequency_summary_dfs: dict[str, pd.DataFrame]
    nucleotide_percentage_summary_dfs: dict[str, pd.DataFrame]
    modification_frequency_summary_dfs: dict[str, pd.DataFrame]
    modification_percentage_summary_dfs: dict[str, pd.DataFrame]

    # Per-amplicon consensus guide data (amplicon_name -> ...)
    consensus_guides: dict[str, list[str]]
    consensus_include_idxs: dict[str, np.ndarray]
    consensus_sgRNA_intervals: dict[str, list[tuple]]
    consensus_sgRNA_plot_idxs: dict[str, list[np.ndarray]]
    guides_all_same: dict[str, bool]

    # Per-amplicon summary file paths (amplicon_name -> {kind: filename})
    # Populated by CRISPRessoBatchCORE; used by built-in plot runner.
    all_summary_filenames: dict[str, dict[str, str]] = field(default_factory=dict)

    # Quantification-window summary filenames (populated when base_editor_output).
    sub_nucleotide_frequency_summary_filename: str = ''
    sub_nucleotide_percentage_summary_filename: str = ''

    # === Scope fields ===

    amplicon_name: Optional[str] = None
    sgRNA_ind: Optional[int] = None
    mod_type: Optional[str] = None


# =============================================================================
# Aggregate context
# =============================================================================


@dataclass
class AggregatePlotContext(PlotContext):
    """View into CRISPRessoAggregateCORE's data for cross-folder plots.

    Same per-amplicon dict structure as :class:`BatchPlotContext`, plus
    ``df_summary_quantification`` for the reads-total / unmod-mod-pcts
    summary plots and ``sample_count`` for pagination.
    """

    # === Run-level data ===

    amplicon_names: list[str]

    # Per-amplicon summary DataFrames (amplicon_name -> DataFrame)
    nucleotide_frequency_summary_dfs: dict[str, pd.DataFrame]
    nucleotide_percentage_summary_dfs: dict[str, pd.DataFrame]
    modification_frequency_summary_dfs: dict[str, pd.DataFrame]
    modification_percentage_summary_dfs: dict[str, pd.DataFrame]

    # Per-amplicon consensus guide data (amplicon_name -> ...)
    consensus_guides: dict[str, list[str]]
    consensus_include_idxs: dict[str, np.ndarray]
    consensus_sgRNA_intervals: dict[str, list[tuple]]
    consensus_sgRNA_plot_idxs: dict[str, list[np.ndarray]]
    guides_all_same: dict[str, bool]

    # Summary quantification (for reads_total / unmod_mod_pcts plots)
    df_summary_quantification: pd.DataFrame

    # Per-amplicon sample count (for pagination)
    sample_count: dict[str, int]

    # Per-amplicon summary file paths (amplicon_name -> {kind: filename}).
    # Populated by CRISPRessoAggregateCORE; used by the shared
    # ``run_builtin_batch_plots`` runner.
    all_summary_filenames: dict[str, dict[str, str]] = field(default_factory=dict)

    # Quantification-window summary filenames (populated when base_editor_output).
    sub_nucleotide_frequency_summary_filename: str = ''
    sub_nucleotide_percentage_summary_filename: str = ''

    # === Scope fields ===

    amplicon_name: Optional[str] = None
    sgRNA_ind: Optional[int] = None
    mod_type: Optional[str] = None


# =============================================================================
# Pooled context
# =============================================================================


@dataclass
class PooledPlotContext(PlotContext):
    """View into CRISPRessoPooledCORE's data for region summary plots.

    Only carries ``df_summary_quantification`` — the two plots
    (reads_total, unmod_mod_pcts) read everything else from ``ctx.args``.
    """

    df_summary_quantification: pd.DataFrame


# =============================================================================
# WGS context
# =============================================================================


@dataclass
class WGSPlotContext(PlotContext):
    """View into CRISPRessoWGSCORE's data for region summary plots.

    Same shape as :class:`PooledPlotContext`; separate class for
    independent evolution and plugin dispatch.
    """

    df_summary_quantification: pd.DataFrame


# =============================================================================
# Compare context
# =============================================================================


@dataclass
class ComparePlotContext(PlotContext):
    """View into CRISPRessoCompareCORE's data for paired-sample plots.

    Holds raw inputs from both samples.  Per-amplicon derived data
    (profiles, Fisher test results, merged allele tables) are computed
    by prep functions, not stored on the context.
    """

    # === Run-level data ===

    amplicon_names: list[str]           # Amplicon names present in both samples
    sample_1_name: str
    sample_2_name: str
    run_info_1: dict                    # Full crispresso2_info from sample 1
    run_info_2: dict                    # Full crispresso2_info from sample 2
    amplicon_info_1: dict[str, dict]    # Per-amplicon quantification info, sample 1
    amplicon_info_2: dict[str, dict]    # Per-amplicon quantification info, sample 2

    # === Scope fields (set by CORE per amplicon iteration) ===

    amplicon_name: Optional[str] = None

    # Per-amplicon loaded data
    profile_1: Optional[np.ndarray] = None
    profile_2: Optional[np.ndarray] = None
    mod_freqs_1: Optional[dict] = None
    mod_freqs_2: Optional[dict] = None
    consensus_sequence: Optional[str] = None
    quant_windows_1: Optional[np.ndarray] = None
    quant_windows_2: Optional[np.ndarray] = None
    cut_points: Optional[list] = None
    sgRNA_intervals: Optional[list] = None

    # Per-mod-type scope
    mod_type: Optional[str] = None

    # Allele comparison data (list of (df1, df2, allele_file_1, allele_file_2) tuples)
    allele_pairs: Optional[list] = None
