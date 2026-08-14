"""Data preparation functions for CRISPResso2 plots.

**Every built-in plot has a corresponding prep function here — no exceptions.**

Each CorePlotContext-based prep takes a single :class:`CorePlotContext`
argument and returns the kwargs dict that the corresponding CRISPRessoPlot
function expects, giving every plot a uniform ``prep(ctx) → csv write →
plot → metadata`` shape in CORE.

These are pure functions — no file I/O, no side effects — with the
exception of ``prep_alleles_around_cut`` and ``prep_base_edit_quilt``
which intentionally mutate their input DataFrames (see their docstrings).

CORE calls these to build plot inputs. CRISPRessoPro can call them
from CorePlotContext to generate plots independently.

This module also contains serialization-boundary helper functions
(``prep_alleles_table``, ``prep_alleles_table_compare``,
``prep_amino_acid_table_for_plot``) that reduce DataFrame size before
dispatching to plot worker threads. These are not CorePlotContext-based —
they take explicit arguments and are called by the CorePlotContext-based
prep functions or directly from CRISPRessoPlot's combined plot functions.

Scope fields on CorePlotContext
-------------------------------
Most functions require ``ctx.ref_name`` to be set (which reference amplicon
is being processed). Functions that operate per-sgRNA additionally require
``ctx.sgRNA_ind``. ``prep_amino_acid_table`` also requires
``ctx.coding_seq_ind``.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)
import os
import re
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

from CRISPResso2 import CRISPRessoShared
from CRISPResso2.plots.plot_context import CorePlotContext


# =============================================================================
# Private helpers (unchanged)
# =============================================================================


def _get_run_info(ctx, key, default=''):
    """Safely get a value from ctx.run_data['running_info']."""
    return ctx.run_data.get('running_info', {}).get(key, default)


def _get_ref_info(ctx, ref_name, key, default=''):
    """Safely get a value from ctx.run_data['results']['refs'][ref_name]."""
    return (ctx.run_data
            .get('results', {})
            .get('refs', {})
            .get(ref_name, {})
            .get(key, default))


def _to_numeric_ignore_columns(df, ignore_columns):
    """Convert DataFrame columns to numeric, ignoring specified columns.

    Uses positional indexing so duplicate column labels (for example the
    repeated base letters used by nucleotide quilt DataFrames) are converted
    exactly once per column rather than repeatedly selecting all columns with
    the same label.
    """
    converted_columns = []
    for idx, col in enumerate(df.columns):
        series = df.iloc[:, idx]
        if col not in ignore_columns:
            series = pd.to_numeric(series, errors='coerce')
        converted_columns.append(series)

    result = pd.concat(converted_columns, axis=1)
    result.columns = df.columns
    return result


def plot_title_with_ref_name(title, ref_name, num_refs):
    """Add ref_name suffix to plot title when there are multiple references."""
    if num_refs > 1:
        return title + ": " + ref_name
    return title


def _clip_to_percentile(values, counts, percentile=0.99, total=None):
    """Find the value at which cumulative counts reach the given percentile.

    Scans ``values`` and ``counts`` in order, returning the value where
    the running sum of counts first exceeds ``percentile * total``.
    If the cutoff is never reached, returns the last value.

    Parameters
    ----------
    total : float, optional
        Override the total used for cutoff calculation. If None, uses
        ``counts.sum()``.

    """
    if total is None:
        total = counts.sum()
    if total == 0:
        return values[-1] if len(values) > 0 else 0
    cutoff = percentile * total
    running = 0
    for val, cnt in zip(values, counts):
        running += cnt
        if running > cutoff:
            return val
    return values[-1]


# =============================================================================
# New private helpers for CorePlotContext extraction
# =============================================================================


def _make_fig_filename_root(ctx: CorePlotContext, filename: str) -> str:
    """Construct fig_filename_root path, or return *filename* if _jp is unavailable."""
    if ctx._jp is not None:
        return ctx._jp(filename)
    return filename


def _data_file_basename(ctx, filename: str) -> str:
    """Return the basename of a data file, including any file_prefix.

    When a ``file_prefix`` is set (e.g. batch runs), ``_jp`` prepends
    the slugified prefix to the filename.  This helper uses ``_jp`` to
    reconstruct the correct on-disk basename so that data-file links in
    the HTML report resolve correctly.
    """
    if ctx._jp is not None:
        return os.path.basename(ctx._jp(filename))
    return filename


def _ref(ctx: CorePlotContext) -> dict:
    """Shortcut for the current reference's data dict."""
    return ctx.refs[ctx.ref_name]


def _ref_plot_name(ctx: CorePlotContext) -> str:
    """Return the current reference's plot-name prefix."""
    return _ref(ctx)['ref_plot_name']


def _sgRNA_label(ctx: CorePlotContext) -> str:
    """Compute the file-name label for the current sgRNA."""
    ref = _ref(ctx)
    sgRNA = ref['sgRNA_orig_sequences'][ctx.sgRNA_ind]
    sgRNA_name = ref['sgRNA_names'][ctx.sgRNA_ind]
    label = "sgRNA_" + sgRNA
    if sgRNA_name != "":
        label = sgRNA_name
    return CRISPRessoShared.slugify(label)


def _sgRNA_legend(ctx: CorePlotContext) -> str:
    """Compute the legend string for the current sgRNA."""
    ref = _ref(ctx)
    sgRNA = ref['sgRNA_orig_sequences'][ctx.sgRNA_ind]
    sgRNA_name = ref['sgRNA_names'][ctx.sgRNA_ind]
    legend = "sgRNA " + sgRNA
    if sgRNA_name != "":
        legend = sgRNA_name + " (" + sgRNA + ")"
    return legend


_CLIPPING_NOTE = (
    " Note that histograms are clipped to show 99% of the data."
    " To show all data, run using the parameter '--plot_histogram_outliers'. "
)


def _build_indel_clipped_string(xmin, xmax, raw_xmin, raw_xmax):
    """Build the clipping annotation for plot_3a (indel size distribution).

    Compares the clipped xmin/xmax to the raw data range and returns a
    string describing which extremes were hidden, or ``""`` if nothing
    was clipped.
    """
    parts = []
    if xmax < raw_xmax:
        parts.append(f" (Maximum {int(raw_xmax)} not shown)")
    if xmin > raw_xmin:
        parts.append(f" (Minimum {int(raw_xmin)} not shown)")
    if parts:
        return _CLIPPING_NOTE + "".join(parts)
    return ""


def _build_freq_del_ins_clipped_string(xmax_ins, xmax_del, xmax_mut,
                                       raw_xmax_ins, raw_xmax_del,
                                       raw_xmax_mut):
    """Build the clipping annotation for plot_3b (ins/del/sub histograms).

    Compares the clipped xmax values to the raw maxima and returns a
    string describing which categories were truncated, or ``""`` if
    nothing was clipped.
    """
    parts = []
    if xmax_ins < raw_xmax_ins:
        parts.append(f" (Insertion maximum {int(raw_xmax_ins)} not shown)")
    if xmax_del < raw_xmax_del:
        parts.append(f" (Deletion minimum -{int(raw_xmax_del)} not shown)")
    if xmax_mut < raw_xmax_mut:
        parts.append(f" (Mutation maximum {int(raw_xmax_mut)} not shown)")
    if parts:
        return _CLIPPING_NOTE + "".join(parts)
    return ""


def _build_nuc_freq_df(ctx: CorePlotContext, ref_name: str | None = None):
    """Build nucleotide frequency and percentage DataFrames for a reference.

    Returns ``(df_nuc_freq, df_nuc_pct)`` where rows are nucleotides
    ``['A', 'C', 'G', 'T', 'N', '-']`` and columns are the reference
    sequence positions.
    """
    if ref_name is None:
        ref_name = ctx.ref_name
    ref_seq = ctx.refs[ref_name]['sequence']
    tot = float(ctx.counts_total[ref_name])

    df_nuc_freq = pd.DataFrame([
        ctx.all_base_count_vectors[ref_name + '_A'],
        ctx.all_base_count_vectors[ref_name + '_C'],
        ctx.all_base_count_vectors[ref_name + '_G'],
        ctx.all_base_count_vectors[ref_name + '_T'],
        ctx.all_base_count_vectors[ref_name + '_N'],
        ctx.all_base_count_vectors[ref_name + '_-'],
    ])
    df_nuc_freq.index = ['A', 'C', 'G', 'T', 'N', '-']
    df_nuc_freq.columns = list(ref_seq)
    df_nuc_pct = df_nuc_freq.divide(tot) if tot > 0 else df_nuc_freq * 0
    return df_nuc_freq, df_nuc_pct


def _build_mod_pct_rows(
    ref_name: str,
    total: float,
    seq_len: int,
    insertion_vectors: dict,
    insertion_left_vectors: dict,
    deletion_vectors: dict,
    substitution_vectors: dict,
    indelsub_vectors: dict,
    counts_total: dict,
    include_batch_label: bool = True,
):
    """Build the 6 modification-percentage rows for one reference.

    Returns a list of 6 ``np.ndarray`` rows (one per modification type).
    When *include_batch_label* is True (default), each row is prefixed
    with ``[ref_name, label, ...]``; when False, each row is prefixed
    with ``[label, ...]`` only.

    The caller can concatenate rows from multiple references and wrap
    them in a DataFrame.

    Used by :func:`prep_nucleotide_quilt` (single-ref, ``all_*`` vectors)
    and :func:`prep_hdr_nucleotide_quilt` (multi-ref, ``ref1_all_*``
    vectors) to avoid duplicating the row-construction logic.
    """
    tot = float(total)
    prefix = [ref_name] if include_batch_label else []
    return [
        np.concatenate((prefix + ['Insertions'], np.array(insertion_vectors[ref_name]).astype(float) / tot)),
        np.concatenate((prefix + ['Insertions_Left'], np.array(insertion_left_vectors[ref_name]).astype(float) / tot)),
        np.concatenate((prefix + ['Deletions'], np.array(deletion_vectors[ref_name]).astype(float) / tot)),
        np.concatenate((prefix + ['Substitutions'], np.array(substitution_vectors[ref_name]).astype(float) / tot)),
        np.concatenate((prefix + ['All_modifications'], np.array(indelsub_vectors[ref_name]).astype(float) / tot)),
        np.concatenate((prefix + ['Total'], [counts_total[ref_name]] * seq_len)),
    ]


def _compute_half_windows(cut_point, plot_window_size, ref_len):
    """Compute asymmetric left/right window sizes, clamping to amplicon bounds.

    Returns ``(left, right, window_truncated)``.
    """
    window_truncated = False
    if cut_point - plot_window_size + 1 >= 0:
        left = plot_window_size
    else:
        left = cut_point + 1
        window_truncated = True
    if cut_point + plot_window_size < ref_len:
        right = plot_window_size
    else:
        right = ref_len - cut_point - 1
        window_truncated = True
    return left, right, window_truncated


# =============================================================================
# Serialization-boundary helpers (not CorePlotContext-based)
#
# These run in the main thread to shrink DataFrames before sending to
# plot worker threads.  They are called by CorePlotContext-based prep
# functions and by CRISPRessoPlot's combined plot functions.
# =============================================================================

# Shared DNA-to-number mapping used by prep_alleles_table and friends.
_DNA_TO_NUMBERS = {'-': 0, 'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 5}

# Shared regex for finding insertions (runs of dashes in the reference).
_RE_FIND_INDELS = re.compile(r"(-*-)")


def _seq_to_numbers(seq):
    """Convert a DNA sequence string to a list of integers."""
    return [_DNA_TO_NUMBERS[x] for x in seq]


# Amino acid alphabet used by prep_amino_acid_table_for_plot.
_AMINO_ACIDS = [
    '*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '', '-',
]
_AA_TO_NUMBERS = {aa: i for i, aa in enumerate(_AMINO_ACIDS)}


def amino_acids_to_numbers(seq):
    """Convert an amino acid sequence to a list of integers.

    Used by :func:`prep_amino_acid_table_for_plot` and by
    ``CRISPRessoPlot.plot_amino_acid_heatmap``.
    """
    return [_AA_TO_NUMBERS[aa] for aa in seq]


def prep_alleles_table(df_alleles, reference_seq, MAX_N_ROWS, MIN_FREQUENCY):
    """Prepare a DataFrame of alleles for heatmap plotting.

    Converts aligned sequences to numeric arrays, detects insertions and
    substitutions, and builds annotation metadata.  Runs in the main thread
    so that only the compact output (not the full DataFrame) is serialised
    to plot worker threads.

    Parameters
    ----------
    df_alleles : pd.DataFrame
        Allele table indexed by ``Aligned_Sequence``, with columns
        ``Reference_Sequence``, ``#Reads``, ``%Reads``.
    reference_seq : str
        The unmodified reference sequence for this window.
    MAX_N_ROWS : int
        Maximum rows to include.
    MIN_FREQUENCY : float
        Minimum ``%Reads`` for a row to be included.

    Returns
    -------
    tuple
        ``(X, annot, y_labels, insertion_dict, per_element_annot_kws,
        is_reference)``

    """
    X = []
    annot = []
    y_labels = []
    insertion_dict = defaultdict(list)
    per_element_annot_kws = []
    is_reference = []

    idx_row = 0
    for idx, row in df_alleles[df_alleles['%Reads'] >= MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(_seq_to_numbers(idx.upper()))
        annot.append(list(idx))

        has_indels = False
        for p in _RE_FIND_INDELS.finditer(row['Reference_Sequence']):
            has_indels = True
            insertion_dict[idx_row].append((p.start(), p.end()))

        y_labels.append('%.2f%% (%d reads)' % (row['%Reads'], row['#Reads']))
        if idx == reference_seq and not has_indels:
            is_reference.append(True)
        else:
            is_reference.append(False)

        idx_row += 1

        idxs_sub = [i_sub for i_sub in range(len(idx)) if
                   (row['Reference_Sequence'][i_sub] != idx[i_sub]) and
                   (row['Reference_Sequence'][i_sub] != '-') and
                   (idx[i_sub] != '-')]
        to_append = np.array([{}] * len(idx), dtype=object)
        to_append[idxs_sub] = {'weight': 'bold', 'color': 'black', 'size': 16}
        per_element_annot_kws.append(to_append)

    return X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference


def prep_alleles_table_compare(df_alleles, sample_name_1, sample_name_2,
                               MAX_N_ROWS, MIN_FREQUENCY):
    """Prepare a merged allele table for comparison heatmap plotting.

    Same shape as :func:`prep_alleles_table` but labels rows with counts
    from both samples.  Used by CRISPRessoCompare.

    Parameters
    ----------
    df_alleles : pd.DataFrame
        Merged allele table with columns ``#Reads_<sample>``,
        ``%Reads_<sample>`` for each sample.
    sample_name_1, sample_name_2 : str
        Sample names (column suffixes).
    MAX_N_ROWS : int
        Maximum rows to include.
    MIN_FREQUENCY : float
        Minimum combined ``%Reads`` for a row to be included.

    Returns
    -------
    tuple
        ``(X, annot, y_labels, insertion_dict, per_element_annot_kws)``

    """
    X = []
    annot = []
    y_labels = []
    insertion_dict = defaultdict(list)
    per_element_annot_kws = []

    idx_row = 0
    for idx, row in df_alleles[df_alleles['%Reads_' + sample_name_1] + df_alleles['%Reads_' + sample_name_2] >= MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(_seq_to_numbers(idx.upper()))
        annot.append(list(idx))
        y_labels.append('%.2f%% (%d reads) %.2f%% (%d reads) ' % (
            row['%Reads_' + sample_name_1], row['#Reads_' + sample_name_1],
            row['%Reads_' + sample_name_2], row['#Reads_' + sample_name_2],
        ))

        for p in _RE_FIND_INDELS.finditer(row['Reference_Sequence']):
            insertion_dict[idx_row].append((p.start(), p.end()))

        idx_row += 1

        idxs_sub = [i_sub for i_sub in range(len(idx)) if
                   (row['Reference_Sequence'][i_sub] != idx[i_sub]) and
                   (row['Reference_Sequence'][i_sub] != '-') and
                   (idx[i_sub] != '-')]
        to_append = np.array([{}] * len(idx), dtype=object)
        to_append[idxs_sub] = {'weight': 'bold', 'color': 'black', 'size': 16}
        per_element_annot_kws.append(to_append)

    return X, annot, y_labels, insertion_dict, per_element_annot_kws


def prep_amino_acid_table_for_plot(df_alleles, reference_seq, MAX_N_ROWS,
                                   MIN_FREQUENCY):
    """Prepare an amino acid allele table for heatmap plotting.

    Same role as :func:`prep_alleles_table` but for amino acid sequences.
    Detects silent edits in addition to insertions and substitutions.

    Parameters
    ----------
    df_alleles : pd.DataFrame
        Amino acid allele table indexed by amino acid sequence, with
        columns ``Reference_Sequence``, ``#Reads``, ``%Reads``,
        ``silent_edit_inds``.
    reference_seq : str
        The unmodified reference amino acid sequence.
    MAX_N_ROWS : int
        Maximum rows to include.
    MIN_FREQUENCY : float
        Minimum ``%Reads`` for a row to be included.

    Returns
    -------
    tuple
        ``(X, annot, y_labels, insertion_dict, silent_edit_dict,
        per_element_annot_kws, is_reference, reference_seq)``

    """
    X = []
    annot = []
    y_labels = []
    insertion_dict = defaultdict(list)
    silent_edit_dict = defaultdict(list)
    per_element_annot_kws = []
    is_reference = []

    idx_row = 0

    for seq, row in df_alleles[df_alleles['%Reads'] >= MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(amino_acids_to_numbers(seq))
        annot.append(list(seq))

        silent_edit_dict[idx_row] = row['silent_edit_inds']

        has_indels = False
        for p in _RE_FIND_INDELS.finditer(row['Reference_Sequence']):
            has_indels = True
            insertion_dict[idx_row].append((p.start(), p.end()))

        y_labels.append('%.2f%% (%d reads)' % (row['%Reads'], row['#Reads']))
        if seq == reference_seq and not has_indels:
            is_reference.append(True)
        else:
            is_reference.append(False)

        idx_row += 1

        idxs_sub = [i_sub for i_sub in range(len(seq)) if
                   (row['Reference_Sequence'][i_sub] != seq[i_sub].upper()) and
                   (row['Reference_Sequence'][i_sub] != '-') and
                   (seq[i_sub] != '-')]
        to_append = np.array([{}] * len(seq), dtype=object)
        to_append[idxs_sub] = {'weight': 'bold', 'color': 'black', 'size': 16}
        per_element_annot_kws.append(to_append)

    for i, (x, a) in enumerate(zip(X, annot)):
        X[i] = x + amino_acids_to_numbers([''] * (len(reference_seq) - len(a)))
        annot[i] = a + [''] * (len(reference_seq) - len(a))

    return X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, reference_seq


# =============================================================================
# Public prep functions — each takes a single CorePlotContext
# =============================================================================


def prep_indel_size_distribution(ctx: CorePlotContext):
    """Prepare kwargs for plot_indel_size_distribution (plot_3a).

    Requires ``ctx.ref_name``.

    Computes xmin/xmax by clipping to the 99th percentile of the
    density distribution (unless ``args.plot_histogram_outliers`` is True),
    then clamping to at least ±15.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)

    hdensity = ref['hdensity']
    hlengths = ref['hlengths']
    center_index = ref['center_index']
    n_this_category = ctx.counts_total[ref_name]
    plot_histogram_outliers = ctx.args.plot_histogram_outliers

    raw_xmin = min(hlengths)
    raw_xmax = max(hlengths)
    xmin = raw_xmin
    xmax = raw_xmax

    if not plot_histogram_outliers:
        xmax = _clip_to_percentile(hlengths, hdensity, 0.99)
        xmin = _clip_to_percentile(hlengths[::-1], hdensity[::-1], 0.99)

    xmin = min(xmin, -15)
    xmax = max(xmax, 15)

    clipped_string = _build_indel_clipped_string(xmin, xmax, raw_xmin, raw_xmax)

    num_refs = len(ctx.ref_names)
    fig_filename_root = _make_fig_filename_root(
        ctx, '3a.' + _ref_plot_name(ctx) + 'Indel_size_distribution',
    )

    indel_hist_filename = _get_ref_info(ctx, ref_name, 'indel_histogram_filename', '')
    return {
        'hdensity': hdensity,
        'hlengths': hlengths,
        'center_index': center_index,
        'n_this_category': n_this_category,
        'xmin': xmin,
        'xmax': xmax,
        'title': plot_title_with_ref_name(
            'Indel size distribution', ref_name, num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'ref_name': ref_name,
        'clipped_string': clipped_string,
        'caption': (
            "Figure 3a: Frequency distribution of alleles with indels (blue) "
            "and without indels (red)." + clipped_string
        ),
        'data_files': [('Indel histogram', os.path.basename(indel_hist_filename))] if indel_hist_filename else [],
    }


def prep_frequency_deletions_insertions(ctx: CorePlotContext):
    """Prepare kwargs for plot_frequency_deletions_insertions (plot_3b).

    Requires ``ctx.ref_name``.

    Clips each histogram's xmax to the 99th percentile of the overall
    density (unless ``args.plot_histogram_outliers`` is True), clamped to 15.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)

    x_bins_ins = ref['x_bins_ins']
    y_values_ins = ref['y_values_ins']
    x_bins_del = ref['x_bins_del']
    y_values_del = ref['y_values_del']
    x_bins_mut = ref['x_bins_mut']
    y_values_mut = ref['y_values_mut']
    hdensity = ref['hdensity']
    counts_total = ctx.counts_total[ref_name]
    plot_histogram_outliers = ctx.args.plot_histogram_outliers

    raw_xmax_ins = max(x_bins_ins)
    raw_xmax_del = max(x_bins_del)
    raw_xmax_mut = max(x_bins_mut)
    xmax_ins = raw_xmax_ins
    xmax_del = raw_xmax_del
    xmax_mut = raw_xmax_mut

    if not plot_histogram_outliers:
        density_total = hdensity.sum()
        xmax_ins = _clip_to_percentile(x_bins_ins, y_values_ins, 0.99, total=density_total)
        xmax_del = _clip_to_percentile(x_bins_del, y_values_del, 0.99, total=density_total)
        xmax_mut = _clip_to_percentile(x_bins_mut, y_values_mut, 0.99, total=density_total)

    xmax_ins = max(15, xmax_ins)
    xmax_del = max(15, xmax_del)
    xmax_mut = max(15, xmax_mut)

    clipped_string = _build_freq_del_ins_clipped_string(
        xmax_ins, xmax_del, xmax_mut,
        raw_xmax_ins, raw_xmax_del, raw_xmax_mut,
    )

    num_refs = len(ctx.ref_names)
    fig_filename_root = _make_fig_filename_root(
        ctx,
        '3b.' + _ref_plot_name(ctx) + 'Insertion_deletion_substitutions_size_hist',
    )

    ins_hist_fn = _get_ref_info(ctx, ref_name, 'insertion_histogram_filename', '')
    del_hist_fn = _get_ref_info(ctx, ref_name, 'deletion_histogram_filename', '')
    sub_hist_fn = _get_ref_info(ctx, ref_name, 'substitution_histogram_filename', '')
    data_files = []
    if ins_hist_fn:
        data_files.append(('Insertions frequency', os.path.basename(ins_hist_fn)))
    if del_hist_fn:
        data_files.append(('Deletions Frequency', os.path.basename(del_hist_fn)))
    if sub_hist_fn:
        data_files.append(('Substitutions Frequency', os.path.basename(sub_hist_fn)))

    return {
        'ref': ref,
        'counts_total': counts_total,
        'fig_filename_root': fig_filename_root,
        'plot_titles': {
            'ins': plot_title_with_ref_name('Insertions', ref_name, num_refs),
            'del': plot_title_with_ref_name('Deletions', ref_name, num_refs),
            'mut': plot_title_with_ref_name('Substitutions', ref_name, num_refs),
        },
        'xmax_ins': xmax_ins,
        'xmax_del': xmax_del,
        'xmax_mut': xmax_mut,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'ref_name': ref_name,
        'clipped_string': clipped_string,
        'caption': (
            "Figure 3b: Left panel, frequency distribution of sequence modifications that "
            "increase read length with respect to the reference amplicon, classified as "
            "insertions (positive indel size). Middle panel, frequency distribution of "
            "sequence modifications that reduce read length with respect to the reference "
            "amplicon, classified as deletions (negative indel size). Right panel, frequency "
            "distribution of sequence modifications that do not alter read length with respect "
            "to the reference amplicon, which are classified as substitutions (number of "
            "substituted positions shown)." + clipped_string
        ),
        'data_files': data_files,
    }


def prep_amplicon_modifications(ctx: CorePlotContext):
    """Prepare kwargs for plot_amplicon_modifications (plot_4a).

    Requires ``ctx.ref_name``.

    Pattern A — packaging with y_max computation and title generation.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    all_indelsub_count_vector = ctx.all_indelsub_count_vectors[ref_name]
    include_idxs_list = ref['include_idxs']
    cut_points = ref['sgRNA_cut_points']
    plot_cut_points = ref['sgRNA_plot_cut_points']
    sgRNA_intervals = ref['sgRNA_intervals']
    ref_len = ref['sequence_length']

    y_max = max(all_indelsub_count_vector) * 1.1

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '4a.' + _ref_plot_name(ctx) + 'Combined_insertion_deletion_substitution_locations',
    )

    return {
        'all_indelsub_count_vectors': all_indelsub_count_vector,
        'include_idxs_list': include_idxs_list,
        'cut_points': cut_points,
        'plot_cut_points': plot_cut_points,
        'sgRNA_intervals': sgRNA_intervals,
        'n_total': ctx.N_TOTAL,
        'n_this_category': ctx.counts_total[ref_name],
        'ref_name': ref_name,
        'num_refs': num_refs,
        'ref_len': ref_len,
        'y_max': y_max,
        'plot_titles': {
            'combined': plot_title_with_ref_name(
                'Combined Insertions/Deletions/Substitutions', ref_name, num_refs,
            ),
            'main': plot_title_with_ref_name(
                'Mutation position distribution', ref_name, num_refs,
            ),
        },
        'fig_filename_root': fig_filename_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
        'caption': "Figure 4a: Combined frequency of any modification across the amplicon. Modifications outside of the quantification window are also shown.",
        'data_files': [('Combined modification frequency', os.path.basename(_get_ref_info(ctx, ref_name, 'combined_pct_vector_filename', '')))] if _get_ref_info(ctx, ref_name, 'combined_pct_vector_filename', '') else [],
    }


def prep_modification_frequency(ctx: CorePlotContext):
    """Prepare kwargs for plot_modification_frequency (plot_4b).

    Requires ``ctx.ref_name``.

    Pattern A — packaging with title generation. ``y_max`` is computed
    from the combined indelsub count vector (same as plot_4a).
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    y_max = max(ctx.all_indelsub_count_vectors[ref_name]) * 1.1

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '4b.' + _ref_plot_name(ctx) + 'Insertion_deletion_substitution_locations',
    )

    return {
        'include_idxs_list': ref['include_idxs'],
        'all_insertion_count_vectors': ctx.all_insertion_count_vectors[ref_name],
        'all_deletion_count_vectors': ctx.all_deletion_count_vectors[ref_name],
        'all_substitution_count_vectors': ctx.all_substitution_count_vectors[ref_name],
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'ref_len': ref['sequence_length'],
        'ref_name': ref_name,
        'num_refs': num_refs,
        'n_total': ctx.N_TOTAL,
        'n_this_category': ctx.counts_total[ref_name],
        'cut_points': ref['sgRNA_cut_points'],
        'plot_cut_points': ref['sgRNA_plot_cut_points'],
        'y_max': y_max,
        'plot_title': plot_title_with_ref_name(
            'Mutation position distribution', ref_name, num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
        'caption': "Figure 4b: Frequency of insertions, deletions, and substitutions across the entire amplicon, including modifications outside of the quantification window.",
        'data_files': [('Modification frequency', os.path.basename(_get_ref_info(ctx, ref_name, 'mod_count_filename', '')))] if _get_ref_info(ctx, ref_name, 'mod_count_filename', '') else [],
    }


def prep_dsODN_piechart(ctx: CorePlotContext):
    """Prepare kwargs for plot_class_dsODN_piechart (plot_1d).

    No scope fields required.

    Computes labels and sizes from df_alleles 'contains dsODN' column.
    """
    df_alleles = ctx.df_alleles
    N_TOTAL = ctx.N_TOTAL
    fig_filename_root = _make_fig_filename_root(ctx, '1d.Detection_of_dsODN')

    n_contain = df_alleles[df_alleles['contains dsODN'] == True]['#Reads'].sum()
    n_not_contain = df_alleles[df_alleles['contains dsODN'] == False]['#Reads'].sum()
    labels = [
        'Contains dsODN\n(' + str(n_contain) + ' reads)',
        'No dsODN\n(' + str(n_not_contain) + ' reads)',
    ]
    sizes = [
        100 * n_contain / float(N_TOTAL),
        100 * n_not_contain / float(N_TOTAL),
    ]
    allele_frequency_table_filename = _get_run_info(ctx, 'allele_frequency_table_filename', '')
    return {
        'sizes': sizes,
        'labels': labels,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'caption': "Figure 1d: Frequency of detection of dsODN " + str(getattr(ctx.args, 'dsODN', '')),
        'data_files': [('Allele table', os.path.basename(allele_frequency_table_filename))] if allele_frequency_table_filename else [],
    }


def prep_nucleotide_quilt(ctx: CorePlotContext):
    """Prepare kwargs for plot_nucleotide_quilt (plot_2a).

    Requires ``ctx.ref_name``.

    Builds ``nuc_pct_df`` and ``mod_pct_df`` from per-position count vectors.

    The returned dict is passed directly to ``plot_nucleotide_quilt``.
    CORE should also read ``nuc_pct_df`` and ``mod_pct_df`` back from the
    returned dict to pass as inputs to ``prep_nucleotide_quilt_around_sgRNA``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    ref_seq = ref['sequence']

    _df_nuc_freq_all, df_nuc_pct_all = _build_nuc_freq_df(ctx)

    mod_pcts = _build_mod_pct_rows(
        ref_name=ref_name,
        total=ctx.counts_total[ref_name],
        seq_len=len(ref_seq),
        insertion_vectors=ctx.all_insertion_count_vectors,
        insertion_left_vectors=ctx.all_insertion_left_count_vectors,
        deletion_vectors=ctx.all_deletion_count_vectors,
        substitution_vectors=ctx.all_substitution_count_vectors,
        indelsub_vectors=ctx.all_indelsub_count_vectors,
        counts_total=ctx.counts_total,
        include_batch_label=False,
    )
    colnames = ['Modification'] + list(ref_seq)
    modification_percentage_summary_df = _to_numeric_ignore_columns(
        pd.DataFrame(mod_pcts, columns=colnames), {'Modification'},
    )

    nuc_df_for_plot = df_nuc_pct_all.reset_index().rename(columns={'index': 'Nucleotide'})
    nuc_df_for_plot.insert(0, 'Batch', ref_name)
    mod_df_for_plot = modification_percentage_summary_df.copy()
    mod_df_for_plot.insert(0, 'Batch', ref_name)

    fig_filename_root = _make_fig_filename_root(
        ctx, '2a.' + _ref_plot_name(ctx) + 'Nucleotide_percentage_quilt',
    )

    nuc_freq_filename = _get_ref_info(ctx, ref_name, 'nuc_freq_filename', '')
    return {
        'nuc_pct_df': nuc_df_for_plot,
        'mod_pct_df': mod_df_for_plot,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'sgRNA_names': ref['sgRNA_names'],
        'sgRNA_mismatches': ref['sgRNA_mismatches'],
        'sgRNA_sequences': ref['sgRNA_sequences'],
        'quantification_window_idxs': ref['include_idxs'],
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': (
            "Figure 2a: Nucleotide distribution across amplicon. At each base in the reference "
            "amplicon, the percentage of each base as observed in sequencing reads is shown "
            "(A = green; C = orange; G = yellow; T = purple). Black bars show the percentage of "
            "reads for which that base was deleted. Brown bars between bases show the percentage "
            "of reads having an insertion at that position."
        ),
        'data_files': [('Nucleotide frequency table', os.path.basename(nuc_freq_filename))] if nuc_freq_filename else [],
    }


def prep_nucleotide_quilt_around_sgRNA(ctx: CorePlotContext):
    """Prepare kwargs for plot_nucleotide_quilt around one sgRNA (plot_2b).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Internally calls :func:`prep_nucleotide_quilt` to obtain the full-amplicon
    DataFrames, then slices them to the window around the cut point and
    adjusts sgRNA interval coordinates to the local frame.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind
    cut_point = ref['sgRNA_cut_points'][sgRNA_ind]
    plot_half_window = max(1, ctx.args.plot_window_size)
    ref_len = ref['sequence_length']
    sgRNA_intervals = ref['sgRNA_intervals']
    include_idxs_list = ref['include_idxs']

    # Build full-amplicon DataFrames
    quilt_data = prep_nucleotide_quilt(ctx)
    nuc_df_for_plot = quilt_data['nuc_pct_df']
    mod_df_for_plot = quilt_data['mod_pct_df']

    # Columns 0 and 1 are 'Batch' and 'Nucleotide'/'Modification' labels;
    # sequence position p maps to column index p + 2.
    new_sel_cols_start = max(2, cut_point - plot_half_window + 1)
    new_sel_cols_end = min(ref_len, cut_point + plot_half_window + 1)
    sel_cols = [0, 1] + list(range(new_sel_cols_start + 2, new_sel_cols_end + 2))

    new_sgRNA_intervals = [
        (int_start - new_sel_cols_start, int_end - new_sel_cols_start)
        for (int_start, int_end) in sgRNA_intervals
    ]
    new_include_idx = [x - new_sel_cols_start for x in include_idxs_list]

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '2b.' + _ref_plot_name(ctx) + 'Nucleotide_percentage_quilt_around_' + _sgRNA_label(ctx),
    )

    sgRNA_legend = _sgRNA_legend(ctx)
    quant_window_nuc_freq_filename = _get_ref_info(ctx, ref_name, 'quant_window_nuc_freq_filename', '')
    return {
        'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
        'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': new_sgRNA_intervals,
        'sgRNA_names': ref['sgRNA_names'],
        'sgRNA_mismatches': ref['sgRNA_mismatches'],
        'sgRNA_sequences': ref['sgRNA_sequences'],
        'quantification_window_idxs': new_include_idx,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': 'Figure 2b: Nucleotide distribution around the ' + sgRNA_legend + '.',
        'data_files': [('Nucleotide frequency in quantification window', os.path.basename(quant_window_nuc_freq_filename))] if quant_window_nuc_freq_filename else [],
    }


def prep_hdr_nucleotide_quilt(ctx: CorePlotContext):
    """Prepare kwargs for plot_nucleotide_quilt for HDR comparison (plot_4g).

    No scope fields required (iterates over all references with reads).

    Builds multi-batch ``nuc_pct_df`` and ``mod_pct_df`` across all references
    with non-zero read counts, aligned to the first reference's coordinate
    system.
    """
    ref_names_for_hdr = [r for r in ctx.ref_names if ctx.counts_total[r] > 0]
    ref0 = ref_names_for_hdr[0]
    ref_seq = ctx.refs[ref0]['sequence']
    seq_len = ctx.refs[ref0]['sequence_length']

    nuc_pcts = []
    for rn in ref_names_for_hdr:
        tot = float(ctx.counts_total[rn])
        for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
            nuc_pcts.append(np.concatenate(
                ([rn, nuc],
                 np.array(ctx.ref1_all_base_count_vectors[rn + '_' + nuc]).astype(float) / tot),
            ))
    colnames = ['Batch', 'Nucleotide'] + list(ref_seq)
    nuc_pct_df = _to_numeric_ignore_columns(
        pd.DataFrame(nuc_pcts, columns=colnames), {'Batch', 'Nucleotide'},
    )

    mod_pcts = []
    for rn in ref_names_for_hdr:
        mod_pcts.extend(_build_mod_pct_rows(
            ref_name=rn,
            total=ctx.counts_total[rn],
            seq_len=seq_len,
            insertion_vectors=ctx.ref1_all_insertion_count_vectors,
            insertion_left_vectors=ctx.ref1_all_insertion_left_count_vectors,
            deletion_vectors=ctx.ref1_all_deletion_count_vectors,
            substitution_vectors=ctx.ref1_all_substitution_count_vectors,
            indelsub_vectors=ctx.ref1_all_indelsub_count_vectors,
            counts_total=ctx.counts_total,
        ))
    colnames = ['Batch', 'Modification'] + list(ref_seq)
    mod_pct_df = _to_numeric_ignore_columns(
        pd.DataFrame(mod_pcts, columns=colnames), {'Batch', 'Modification'},
    )

    fig_filename_root = _make_fig_filename_root(ctx, '4g.HDR_nucleotide_percentage_quilt')

    return {
        'nuc_pct_df': nuc_pct_df,
        'mod_pct_df': mod_pct_df,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': ctx.refs[ref0]['sgRNA_intervals'],
        'sgRNA_names': ctx.refs[ref0]['sgRNA_names'],
        'sgRNA_mismatches': ctx.refs[ref0]['sgRNA_mismatches'],
        'sgRNA_sequences': ctx.refs[ref0]['sgRNA_sequences'],
        'quantification_window_idxs': [],  # windows may differ between amplicons
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': (
            "Figure 4g: Nucleotide distribution across all amplicons. At each base in the "
            "reference amplicon, the percentage of each base as observed in sequencing reads "
            "is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the "
            "percentage of reads for which that base was deleted. Brown bars between bases "
            "show the percentage of reads having an insertion at that position."
        ),
        'data_files': [],
    }


def prep_pe_nucleotide_quilt(ctx: CorePlotContext):
    """Prepare kwargs for plot_nucleotide_quilt for PE comparison (plot_11a).

    No scope fields required.

    Same computation as ``prep_hdr_nucleotide_quilt`` but uses the first
    reference's ``include_idxs`` for ``quantification_window_idxs`` (PE
    amplicons share the same quantification window, unlike HDR).
    """
    result = prep_hdr_nucleotide_quilt(ctx)
    result['quantification_window_idxs'] = ctx.refs[ctx.ref_names[0]]['include_idxs']
    result['fig_filename_root'] = _make_fig_filename_root(
        ctx, '11a.Prime_editing_nucleotide_percentage_quilt',
    )
    result['caption'] = (
        "Figure 11a: Nucleotide distribution across all amplicons. At each base in the "
        "reference amplicon, the percentage of each base as observed in sequencing reads "
        "is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the "
        "percentage of reads for which that base was deleted. Brown bars between bases "
        "show the percentage of reads having an insertion at that position."
    )
    result['data_files'] = []
    return result


def prep_pe_nucleotide_quilt_around_sgRNA(ctx: CorePlotContext):
    """Prepare kwargs for plot_nucleotide_quilt around one sgRNA for PE (plot_11b).

    Requires ``ctx.sgRNA_ind``. Uses the first reference for sgRNA metadata.

    Internally calls :func:`prep_pe_nucleotide_quilt` to obtain the multi-batch
    DataFrames, then applies the same window-slicing as
    :func:`prep_nucleotide_quilt_around_sgRNA`.
    """
    ref0_name = ctx.ref_names[0]
    ref0 = ctx.refs[ref0_name]
    sgRNA_ind = ctx.sgRNA_ind
    cut_point = ref0['sgRNA_cut_points'][sgRNA_ind]
    plot_half_window = max(1, ctx.args.plot_window_size)
    ref_len = ref0['sequence_length']
    sgRNA_intervals = ref0['sgRNA_intervals']
    include_idxs_list = ref0['include_idxs']

    # Build full-amplicon PE DataFrames
    pe_data = prep_pe_nucleotide_quilt(ctx)
    nuc_df_for_plot = pe_data['nuc_pct_df']
    mod_df_for_plot = pe_data['mod_pct_df']

    # Window slicing (same logic as prep_nucleotide_quilt_around_sgRNA)
    new_sel_cols_start = max(2, cut_point - plot_half_window + 1)
    new_sel_cols_end = min(ref_len, cut_point + plot_half_window + 1)
    sel_cols = [0, 1] + list(range(new_sel_cols_start + 2, new_sel_cols_end + 2))

    new_sgRNA_intervals = [
        (int_start - new_sel_cols_start, int_end - new_sel_cols_start)
        for (int_start, int_end) in sgRNA_intervals
    ]
    new_include_idx = [x - new_sel_cols_start for x in include_idxs_list]

    # Build sgRNA label from first reference
    sgRNA = ref0['sgRNA_orig_sequences'][sgRNA_ind]
    sgRNA_name = ref0['sgRNA_names'][sgRNA_ind]
    label = "sgRNA_" + sgRNA
    if sgRNA_name != "":
        label = sgRNA_name
    label = CRISPRessoShared.slugify(label)

    # Human-readable legend for caption (with spaces, not slugified)
    legend = "sgRNA " + sgRNA
    if sgRNA_name != "":
        legend = sgRNA_name + " (" + sgRNA + ")"

    fig_filename_root = _make_fig_filename_root(
        ctx, '11b.Nucleotide_percentage_quilt_around_' + label,
    )

    # Build data_files from quant_window_nuc_freq_filename for each PE ref
    pe_ref_names = [r for r in ctx.ref_names if ctx.counts_total[r] > 0]
    pe_data_files = []
    for rn in pe_ref_names:
        qw_nuc_fn = _get_ref_info(ctx, rn, 'quant_window_nuc_freq_filename', '')
        if qw_nuc_fn:
            pe_data_files.append(('Nucleotide frequency in quantification window for ' + rn, os.path.basename(qw_nuc_fn)))

    return {
        'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
        'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': new_sgRNA_intervals,
        'sgRNA_names': pe_data['sgRNA_names'],
        'sgRNA_mismatches': pe_data['sgRNA_mismatches'],
        'sgRNA_sequences': pe_data['sgRNA_sequences'],
        'quantification_window_idxs': new_include_idx,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': 'Figure 11b: Nucleotide distribution around the ' + legend + '.',
        'data_files': pe_data_files,
    }


def prep_global_frameshift_data(ctx: CorePlotContext):
    """Aggregate frameshift/splice/inframe counts across all coding-seq refs.

    No scope fields required.

    Used by plots 5a (global frameshift pie), 6a (global frameshift profiles),
    and 8a (global splice sites). For HDR refs, unmodified reads are added to
    the non-modified-non-frameshift bucket.

    Returns a dict with all aggregated values needed by the three plot
    functions.
    """
    global_MODIFIED_FRAMESHIFT = 0
    global_MODIFIED_NON_FRAMESHIFT = 0
    global_NON_MODIFIED_NON_FRAMESHIFT = 0
    global_SPLICING_SITES_MODIFIED = 0

    global_hists_frameshift = Counter()
    global_hists_frameshift[0] = 0
    global_hists_inframe = Counter()
    global_hists_inframe[0] = 0

    global_count_total = 0
    global_count_modified = 0
    global_count_unmodified = 0

    for ref_name in ctx.ref_names:
        if ctx.refs[ref_name]['contains_coding_seq']:
            global_MODIFIED_FRAMESHIFT += ctx.counts_modified_frameshift[ref_name]
            global_MODIFIED_NON_FRAMESHIFT += ctx.counts_modified_non_frameshift[ref_name]
            global_NON_MODIFIED_NON_FRAMESHIFT += ctx.counts_non_modified_non_frameshift[ref_name]
            global_SPLICING_SITES_MODIFIED += ctx.counts_splicing_sites_modified[ref_name]

            if ref_name == "HDR":
                # for HDR, add all unmodified reads to those that have
                # modifications not in exons
                global_NON_MODIFIED_NON_FRAMESHIFT += ctx.counts_unmodified[ref_name]

            for exon_len, count in ctx.hists_frameshift[ref_name].items():
                global_hists_frameshift[exon_len] += count
            for exon_len, count in ctx.hists_inframe[ref_name].items():
                global_hists_inframe[exon_len] += count

            global_count_total += ctx.counts_total[ref_name]
            global_count_modified += ctx.counts_modified[ref_name]
            global_count_unmodified += ctx.counts_unmodified[ref_name]

    return {
        'global_modified_frameshift': global_MODIFIED_FRAMESHIFT,
        'global_modified_non_frameshift': global_MODIFIED_NON_FRAMESHIFT,
        'global_non_modified_non_frameshift': global_NON_MODIFIED_NON_FRAMESHIFT,
        'global_splicing_sites_modified': global_SPLICING_SITES_MODIFIED,
        'global_hists_frameshift': global_hists_frameshift,
        'global_hists_inframe': global_hists_inframe,
        'global_count_total': global_count_total,
        'global_count_modified': global_count_modified,
        'global_count_unmodified': global_count_unmodified,
    }


def prep_global_modifications_reference(ctx: CorePlotContext):
    """Prepare kwargs for plot_global_modifications_reference (plot_4e/4f).

    Requires ``ctx.ref_name`` (either the primary amplicon or ``"HDR"``).

    The plot root and title are determined automatically based on whether
    ``ctx.ref_name`` is the first reference (4e) or ``"HDR"`` (4f).
    """
    ref_name = ctx.ref_name
    ref0 = ctx.ref_names[0]

    if ref_name == ref0:
        fig_filename_root = _make_fig_filename_root(ctx, '4e.' + ref0 + '.Global_mutations_in_all_reads')
        plot_title = 'Mutation position distribution in all reads with reference to %s' % ref0
    else:
        fig_filename_root = _make_fig_filename_root(
            ctx,
            '4f.' + ref0 + '.Global_mutations_in_HDR_reads_with_reference_to_' + ref0,
        )
        plot_title = 'Mutation position distribution in %s reads with reference to %s' % (ref_name, ref0)

    return {
        'ref1_all_insertion_count_vectors': ctx.ref1_all_insertion_count_vectors[ref_name],
        'ref1_all_deletion_count_vectors': ctx.ref1_all_deletion_count_vectors[ref_name],
        'ref1_all_substitution_count_vectors': ctx.ref1_all_substitution_count_vectors[ref_name],
        'ref1': ctx.refs[ref0],
        'include_idxs_list': _ref(ctx)['include_idxs'],
        'n_total': ctx.N_TOTAL,
        'ref_len': _ref(ctx)['sequence_length'],
        'ref_name': ref0,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
        'plot_title': plot_title,
        'fig_filename_root': fig_filename_root,
        'caption': (
            ("Figure 4e: Positions of modifications in all reads when aligned to the reference "
             "sequence (" + ref0 + "). Insertions: red, deletions: purple, substitutions: green. "
             "All modifications (including those outside the quantification window) are shown.")
            if ref_name == ref0 else
            ("Figure 4f: Positions of modifications in HDR reads with respect to the reference "
             "sequence (" + ref0 + "). All modifications (including those outside the "
             "quantification window) are shown.")
        ),
        'data_files': [],
    }


def prep_log_nuc_freqs(ctx: CorePlotContext):
    """Prepare kwargs for plot_log_nuc_freqs (plot_10d).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Computes ``plot_quant_window_idxs`` — the indices in the plotting window
    that overlap with the quantification window.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind

    include_idxs_list = ref['include_idxs']
    plot_idxs = ref['sgRNA_plot_idxs'][sgRNA_ind]
    ref_len = ref['sequence_length']
    tot_aln_reads = ctx.counts_total[ref_name]

    # Build nucleotide frequency DataFrame and slice to plot window
    df_nuc_freq_all, _df_nuc_pct_all = _build_nuc_freq_df(ctx)
    df_nuc_freq = df_nuc_freq_all.iloc[:, plot_idxs]

    is_window = np.zeros(ref_len)
    for include_idx in include_idxs_list:
        is_window[include_idx] = 1

    plot_quant_window_idxs = []
    for plot_ind, loc in enumerate(plot_idxs):
        if is_window[loc]:
            plot_quant_window_idxs.append(plot_ind - 2)

    num_refs = len(ctx.ref_names)
    sgRNA_leg = _sgRNA_legend(ctx)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '10d.' + _ref_plot_name(ctx) + 'Log2_nucleotide_frequency_around_' + _sgRNA_label(ctx),
    )

    return {
        'df_nuc_freq': df_nuc_freq,
        'tot_aln_reads': tot_aln_reads,
        'plot_title': plot_title_with_ref_name(
            'Log2 Nucleotide Frequencies Around the ' + sgRNA_leg,
            ref_name,
            num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'quantification_window_idxs': plot_quant_window_idxs,
        'caption': (
            "Figure 10d: Log2 nucleotide frequencies for each position in the plotting window "
            "around the " + sgRNA_leg + ". The quantification window is outlined by the dotted box."
        ),
        'data_files': [('Nucleotide frequencies in quantification window', os.path.basename(_get_ref_info(ctx, ref_name, 'quant_window_nuc_freq_filename', '')))] if _get_ref_info(ctx, ref_name, 'quant_window_nuc_freq_filename', '') else [],
    }


def find_closest_sgRNA_cut_point(exon_start, exon_end, sgRNA_cut_points, sgRNA_plot_cut_points):
    """Find the sgRNA with cut_point closest to the given exon interval.

    If no sgRNA cut points are provided, defaults to the exon midpoint.

    Parameters
    ----------
    exon_start : int
        Start position of the exon interval.
    exon_end : int
        End position of the exon interval.
    sgRNA_cut_points : list of int
        List of sgRNA cut point positions.
    sgRNA_plot_cut_points : list
        List of sgRNA plot cut point values (corresponding to sgRNA_cut_points).

    Returns
    -------
    tuple
        (best_cut_point, best_plot_cut_point) where best_cut_point is the
        closest sgRNA cut point (or exon midpoint if no guides), and
        best_plot_cut_point is the corresponding plot cut point (or False if
        no guides).

    """
    best_cut_point = (exon_start + exon_end) // 2  # default to exon midpoint if no guides
    best_plot_cut_point = False
    if len(sgRNA_cut_points) > 0:
        best_sgRNA_idx = 0
        best_distance = float('inf')
        for sgRNA_idx, cp in enumerate(sgRNA_cut_points):
            if exon_start <= cp <= exon_end:
                distance = 0
            else:
                distance = min(abs(cp - exon_start), abs(cp - exon_end))
            if distance < best_distance:
                best_distance = distance
                best_sgRNA_idx = sgRNA_idx
        best_cut_point = sgRNA_cut_points[best_sgRNA_idx]
        best_plot_cut_point = sgRNA_plot_cut_points[best_sgRNA_idx]
    return best_cut_point, best_plot_cut_point


def prep_amino_acid_table(ctx: CorePlotContext):
    """Prepare amino acid table data for plot_9a and CSV export.

    Requires ``ctx.ref_name`` and ``ctx.coding_seq_ind``.

    Converts a nucleotide coding sequence to amino acids, finds the closest
    sgRNA cut point to the coding sequence's exon interval, computes the
    amino acid cut point, builds the amino acid allele DataFrame via
    ``CRISPRessoShared.get_amino_acid_dataframe``, and calls
    :func:`prep_amino_acid_table_for_plot` to produce a
    serialization-friendly representation for the plot worker thread.

    Returns a dict with:

    - ``coding_seq_amino_acids``: amino acid string for the coding sequence
    - ``amino_acid_cut_point``: cut point in amino acid coordinates
    - ``df_to_plot``: DataFrame of amino acid alleles (for CSV export)
    - ``plot_input``: dict of kwargs for ``plot_amino_acid_heatmap``, or
      ``None`` if no rows pass the frequency threshold
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    coding_seq_ind = ctx.coding_seq_ind

    coding_seqs = ctx.run_data['running_info']['coding_seqs']
    coding_seq = coding_seqs[coding_seq_ind]
    coding_seq_names = ctx.run_data['running_info'].get('coding_seq_names', [str(i) for i in range(len(coding_seqs))])
    coding_seq_label = coding_seq_names[coding_seq_ind]
    sgRNA_cut_points = ref['sgRNA_cut_points']
    sgRNA_plot_cut_points = ref['sgRNA_plot_cut_points']
    df_alleles_for_ref = ctx.df_alleles.loc[ctx.df_alleles['Reference_Name'] == ref_name]
    exon_interval_start = ref['exon_intervals'][coding_seq_ind][0]

    # Find the sgRNA with cut_point closest to this coding sequence's exon interval
    exon_start, exon_end = ref['exon_intervals'][coding_seq_ind]
    best_cut_point, best_plot_cut_point = find_closest_sgRNA_cut_point(
        exon_start, exon_end, sgRNA_cut_points, sgRNA_plot_cut_points,
    )

    _ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    blosum_path = os.path.join(_ROOT, "BLOSUM62")

    coding_seq_amino_acids = CRISPRessoShared.get_amino_acids_from_nucs(coding_seq)
    amino_acid_cut_point = (best_cut_point - exon_start + 1) // 3
    amino_acid_cut_point = max(0, min(amino_acid_cut_point, len(coding_seq_amino_acids) - 1))
    df_to_plot = CRISPRessoShared.get_amino_acid_dataframe(
        df_alleles_for_ref,
        exon_interval_start,
        len(coding_seq_amino_acids),
        blosum_path,
        amino_acid_cut_point,
    )

    # Build serialization-friendly plot input
    n_good = df_to_plot[
        df_to_plot['%Reads'] >= ctx.args.min_frequency_alleles_around_cut_to_plot
    ].shape[0]

    plot_input = None
    if n_good > 0:
        X, annot, y_labels, insertion_dict, silent_edit_dict, per_element_annot_kws, is_reference, ref_seq_aa = \
            prep_amino_acid_table_for_plot(
                df_to_plot,
                coding_seq_amino_acids,
                ctx.args.max_rows_alleles_around_cut_to_plot,
                ctx.args.min_frequency_alleles_around_cut_to_plot,
            )

        # Apply wildtype allele annotation (same logic as plot_amino_acid_table bundle)
        annotate_wt = ctx.args.annotate_wildtype_allele
        if annotate_wt != '':
            for ix, is_ref in enumerate(is_reference):
                if is_ref:
                    y_labels[ix] += annotate_wt

        fig_filename_root = _make_fig_filename_root(
            ctx,
            '9a.' + _ref_plot_name(ctx) + 'amino_acid_table_around_' + coding_seq_label,
        )

        plot_input = {
            'reference_seq_amino_acids': ref_seq_aa,
            'fig_filename_root': fig_filename_root,
            'X': X,
            'annot': annot,
            'y_labels': y_labels,
            'insertion_dict': insertion_dict,
            'silent_edit_dict': silent_edit_dict,
            'per_element_annot_kws': per_element_annot_kws,
            'custom_colors': ctx.custom_config.get('colors', {}),
            'SAVE_ALSO_PNG': ctx.save_png,
            'plot_cut_point': best_plot_cut_point,
            'annotate_wildtype_allele': ctx.args.annotate_wildtype_allele,
            'amino_acid_cut_point': amino_acid_cut_point,
        }

    return {
        'coding_seq_amino_acids': coding_seq_amino_acids,
        'amino_acid_cut_point': amino_acid_cut_point,
        'df_to_plot': df_to_plot,
        'plot_input': plot_input,
        'caption': (
            "Figure 9a: Visualization of the distribution of identified amino acids "
            "based on the coding sequence " + coding_seq_label + " (" + coding_seq + "). "
            "The vertical dashed line indicates the predicted cleavage site."
        ),
        'data_files': [],
    }


def _prep_windowed_alleles(
    df_alleles_around_cut,
    cut_point,
    window_left,
    window_right,
    ref_sequence,
    sgRNA_intervals,
    count_total,
    allele_plot_pcts_only_for_assigned_reference,
    expand_allele_plots_by_quantification,
):
    """Shared logic for ``prep_alleles_around_cut`` and ``prep_base_edit_quilt``.

    .. warning::

        **Mutates** *df_alleles_around_cut* in place when
        *allele_plot_pcts_only_for_assigned_reference* is True (adds a
        ``%AllReads`` column and overwrites ``%Reads``).  CORE relies on
        this mutation — it reads the modified DataFrame back from the
        returned dict for CSV writing.

    Applies:

    1. Optional percentage adjustment when
       *allele_plot_pcts_only_for_assigned_reference* is True
    2. Reference sequence slicing for the window
    3. Optional groupby collapse (when *expand_allele_plots_by_quantification*
       is False)
    4. sgRNA interval coordinate adjustment to the local window frame

    Returns ``(df_alleles_around_cut, df_to_plot, ref_seq_around_cut,
    new_sgRNA_intervals, new_sel_cols_start)``.
    """
    if allele_plot_pcts_only_for_assigned_reference:
        df_alleles_around_cut['%AllReads'] = df_alleles_around_cut['%Reads']
        df_alleles_around_cut['%Reads'] = df_alleles_around_cut['#Reads'] / count_total * 100

    ref_seq_around_cut = ref_sequence[
        cut_point - window_left + 1:cut_point + window_right + 1
    ]

    df_to_plot = df_alleles_around_cut
    if not expand_allele_plots_by_quantification:
        df_to_plot = df_alleles_around_cut.groupby(
            ['Aligned_Sequence', 'Reference_Sequence'],
        ).sum().reset_index().set_index('Aligned_Sequence')
        df_to_plot.sort_values(
            by=['#Reads', 'Aligned_Sequence', 'Reference_Sequence'],
            inplace=True,
            ascending=[False, True, True],
        )

    new_sgRNA_intervals = []
    new_sel_cols_start = cut_point - window_left
    for (int_start, int_end) in sgRNA_intervals:
        new_sgRNA_intervals.append(
            (int_start - new_sel_cols_start - 1, int_end - new_sel_cols_start - 1),
        )

    return (
        df_alleles_around_cut,
        df_to_plot,
        ref_seq_around_cut,
        new_sgRNA_intervals,
        new_sel_cols_start,
    )


def prep_alleles_around_cut(ctx: CorePlotContext):
    """Prepare alleles-around-cut data for plot_9 and CSV export.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Computes asymmetric window sizes, calls
    ``CRISPRessoShared.get_dataframe_around_cut_asymmetrical``, then applies
    percentage adjustment, optional groupby collapse, and sgRNA interval
    recomputation via ``_prep_windowed_alleles``.

    Internally calls :func:`prep_alleles_table` to produce a
    serialization-friendly representation for the plot worker thread.
    If no rows pass the frequency threshold, ``plot_input`` is ``None``.

    .. warning::

        **Mutates** the intermediate ``df_alleles_around_cut`` in place when
        ``args.allele_plot_pcts_only_for_assigned_reference`` is True.
        See ``_prep_windowed_alleles`` for details.

    Returns a dict with:

    - ``df_alleles_around_cut``: the (possibly mutated) alleles DataFrame
    - ``ref_seq_around_cut``: reference sequence in the window
    - ``new_sgRNA_intervals``: sgRNA intervals in local coordinates
    - ``new_cut_point``: cut point in local coordinates, or None if the cut
      point is not inside any sgRNA interval
    - ``window_truncated``: True if the window was truncated due to proximity
      to amplicon edges
    - ``plot_input``: dict of kwargs for ``plot_alleles_table_prepped``, or
      ``None`` if no rows pass the frequency threshold
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind

    cut_point = ref['sgRNA_cut_points'][sgRNA_ind]
    plot_cut_point = ref['sgRNA_plot_cut_points'][sgRNA_ind]
    ref_len = ref['sequence_length']
    plot_window_size = ctx.args.plot_window_size

    plot_half_window_left, plot_half_window_right, window_truncated = \
        _compute_half_windows(cut_point, plot_window_size, ref_len)

    df_alleles_around_cut = CRISPRessoShared.get_dataframe_around_cut_asymmetrical(
        ctx.df_alleles.loc[ctx.df_alleles['Reference_Name'] == ref_name],
        cut_point,
        plot_half_window_left,
        plot_half_window_right,
    )

    (
        df_alleles_around_cut,
        df_to_plot,
        ref_seq_around_cut,
        new_sgRNA_intervals,
        new_sel_cols_start,
    ) = _prep_windowed_alleles(
        df_alleles_around_cut=df_alleles_around_cut,
        cut_point=cut_point,
        window_left=plot_half_window_left,
        window_right=plot_half_window_right,
        ref_sequence=ref['sequence'],
        sgRNA_intervals=ref['sgRNA_intervals'],
        count_total=ctx.counts_total[ref_name],
        allele_plot_pcts_only_for_assigned_reference=ctx.args.allele_plot_pcts_only_for_assigned_reference,
        expand_allele_plots_by_quantification=ctx.args.expand_allele_plots_by_quantification,
    )

    new_cut_point = None
    for (int_start, int_end) in ref['sgRNA_intervals']:
        if int_start <= cut_point <= int_end:
            new_cut_point = cut_point - new_sel_cols_start - 1

    # rebase onto the window using the SAME convention as new_sgRNA_intervals
    # and new_cut_point above (note the extra -1); plain subtraction, as the
    # nucleotide quilt uses, would land this one column to the right here
    new_include_idx = [x - new_sel_cols_start - 1 for x in ref['include_idxs']]

    # Build serialization-friendly plot input via prep_alleles_table
    n_good = df_alleles_around_cut[
        df_alleles_around_cut['%Reads'] >= ctx.args.min_frequency_alleles_around_cut_to_plot
    ].shape[0]

    plot_input = None
    if n_good > 0:
        prepped_alleles, annotations, y_labels, insertion_dict, per_element_annot_kws, is_reference = prep_alleles_table(
            df_to_plot,
            ref_seq_around_cut,
            ctx.args.max_rows_alleles_around_cut_to_plot,
            ctx.args.min_frequency_alleles_around_cut_to_plot,
        )

        fig_filename_root = _make_fig_filename_root(
            ctx,
            '9.' + _ref_plot_name(ctx) + 'Alleles_frequency_table_around_' + _sgRNA_label(ctx),
        )

        plot_input = {
            'reference_seq': ref_seq_around_cut,
            'prepped_df_alleles': prepped_alleles,
            'annotations': annotations,
            'y_labels': y_labels,
            'insertion_dict': insertion_dict,
            'per_element_annot_kws': per_element_annot_kws,
            'is_reference': is_reference,
            'fig_filename_root': fig_filename_root,
            'custom_colors': ctx.custom_config.get('colors', {}),
            'SAVE_ALSO_PNG': ctx.save_png,
            'plot_cut_point': plot_cut_point,
            'cut_point_ind': new_cut_point if window_truncated else None,
            'sgRNA_intervals': new_sgRNA_intervals,
            'sgRNA_names': ref['sgRNA_names'],
            'sgRNA_mismatches': ref['sgRNA_mismatches'],
            'annotate_wildtype_allele': ctx.args.annotate_wildtype_allele,
        }

    sgRNA_legend = _sgRNA_legend(ctx)
    return {
        'df_alleles_around_cut': df_alleles_around_cut,
        'ref_seq_around_cut': ref_seq_around_cut,
        'new_sgRNA_intervals': new_sgRNA_intervals,
        'quantification_window_idxs': new_include_idx,
        'new_cut_point': new_cut_point,
        'window_truncated': window_truncated,
        'plot_input': plot_input,
        # exposed so CRISPRessoPro can compose its own caption for the
        # interactive table without re-deriving the legend format
        'sgRNA_legend': sgRNA_legend,
        # NOTE: this caption describes the STATIC matplotlib figure, which really
        # does draw red insertion rectangles, bold substitutions and a dashed
        # cleavage line. The Pro interactive table draws none of those and
        # supplies its own caption -- see CRISPRessoPro.plots.data_funcs.
        'caption': (
            "Figure 9: Visualization of the distribution of identified alleles around the "
            "cleavage site for the " + sgRNA_legend + ". Nucleotides are indicated by unique "
            "colors (A = green; C = red; G = yellow; T = purple). Substitutions are shown in "
            "bold font. Red rectangles highlight inserted sequences. Horizontal dashed lines "
            "indicate deleted sequences. The vertical dashed line indicates the predicted "
            "cleavage site."
        ),
        'data_files': [('Allele frequency table', _data_file_basename(ctx, _ref_plot_name(ctx) + 'Alleles_frequency_table_around_' + _sgRNA_label(ctx) + '.txt'))],
    }


def prep_base_edit_quilt(ctx: CorePlotContext):
    """Prepare base edit quilt data for plot_10h and CSV export.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Internally calls ``CRISPRessoShared.get_base_edit_dataframe_around_cut``
    to produce the sliced DataFrame, then applies the same windowing logic
    as ``prep_alleles_around_cut`` (via ``_prep_windowed_alleles``).

    Internally calls :func:`prep_alleles_table` to produce a
    serialization-friendly representation for the plot worker thread.
    If no rows pass the frequency threshold, ``plot_input`` is ``None``.

    Uses a symmetric window (``args.plot_window_size``) and computes
    ``x_labels`` — the 1-indexed positions of the conversion nucleotide in
    the full reference sequence.

    .. warning::

        **Mutates** the intermediate ``df_alleles_around_cut`` in place when
        ``args.allele_plot_pcts_only_for_assigned_reference`` is True.
        See ``_prep_windowed_alleles`` for details.

    Returns a dict with:

    - ``df_alleles_around_cut``: the (possibly mutated) alleles DataFrame
    - ``ref_seq_around_cut``: reference sequence in the window
    - ``plot_input``: dict of kwargs for ``plot_alleles_table_prepped``, or
      ``None`` if no rows pass the frequency threshold
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)

    cut_point = ref['sgRNA_cut_points'][ctx.sgRNA_ind]
    plot_half_window = max(1, ctx.args.plot_window_size)
    conversion_nuc_from = ctx.args.conversion_nuc_from

    df_alleles_around_cut = CRISPRessoShared.get_base_edit_dataframe_around_cut(
        ctx.df_alleles.loc[ctx.df_alleles['Reference_Name'] == ref_name],
        conversion_nuc_from,
    )

    (
        df_alleles_around_cut,
        df_to_plot,
        ref_seq_around_cut,
        new_sgRNA_intervals,
        _new_sel_cols_start,
    ) = _prep_windowed_alleles(
        df_alleles_around_cut=df_alleles_around_cut,
        cut_point=cut_point,
        window_left=plot_half_window,
        window_right=plot_half_window,
        ref_sequence=ref['sequence'],
        sgRNA_intervals=ref['sgRNA_intervals'],
        count_total=ctx.counts_total[ref_name],
        allele_plot_pcts_only_for_assigned_reference=ctx.args.allele_plot_pcts_only_for_assigned_reference,
        expand_allele_plots_by_quantification=ctx.args.expand_allele_plots_by_quantification,
    )

    x_labels = [
        ind for ind, a in enumerate(ref['sequence'], start=1)
        if a == conversion_nuc_from
    ]

    # Build serialization-friendly plot input via prep_alleles_table
    n_good = df_alleles_around_cut[
        df_alleles_around_cut['%Reads'] >= ctx.args.min_frequency_alleles_around_cut_to_plot
    ].shape[0]

    plot_input = None
    if n_good > 0:
        prepped_alleles, annotations, y_labels, insertion_dict, per_element_annot_kws, is_reference = prep_alleles_table(
            df_to_plot,
            ref_seq_around_cut,
            ctx.args.max_rows_alleles_around_cut_to_plot,
            ctx.args.min_frequency_alleles_around_cut_to_plot,
        )

        fig_filename_root = _make_fig_filename_root(
            ctx,
            '10h.' + _ref_plot_name(ctx) + 'base_edit_' + conversion_nuc_from + 's_quilt',
        )

        plot_input = {
            'reference_seq': ref_seq_around_cut,
            'prepped_df_alleles': prepped_alleles,
            'annotations': annotations,
            'y_labels': y_labels,
            'insertion_dict': insertion_dict,
            'per_element_annot_kws': per_element_annot_kws,
            'is_reference': is_reference,
            'fig_filename_root': fig_filename_root,
            'custom_colors': ctx.custom_config.get('colors', {}),
            'SAVE_ALSO_PNG': ctx.save_png,
            'plot_cut_point': None,
            'sgRNA_intervals': None,
            'sgRNA_names': None,
            'sgRNA_mismatches': None,
            'annotate_wildtype_allele': '',
            'plot_reference_sequence_above': False,
            'x_labels': x_labels,
        }

    conversion_nuc_from = ctx.args.conversion_nuc_from
    return {
        'df_alleles_around_cut': df_alleles_around_cut,
        'ref_seq_around_cut': ref_seq_around_cut,
        'plot_input': plot_input,
        'caption': (
            "Figure 10h: Quilt of target nucleotide: " + conversion_nuc_from
            + " across entire amplicon. The x-axis shows the corresponding position of the "
            "nucleotide in the reference amplicon (1-indexed). Nucleotides are indicated by "
            "unique colors (A = green; C = red; G = yellow; T = purple). Substitutions are "
            "shown in bold font. Red rectangles highlight inserted sequences. Horizontal "
            "dashed lines indicate deleted sequences."
        ),
        'data_files': [],
    }


def prep_alternate_allele_counts(sub_base_vectors, ref_name, ref_sequence):
    """Aggregate per-position substitution counts by reference nucleotide.

    For each reference base, sums the substitution counts across all
    positions where that base appears. Returns a nested dict
    ``{ref_nuc: {obs_nuc: count}}``.

    This is a low-level utility (not CorePlotContext-based) used by
    ``count_alternate_alleles`` in CORE.
    """
    alph = ['A', 'C', 'G', 'T', 'N']
    alt_nuc_counts = {a: {b: 0 for b in alph} for a in alph}
    for pos, ref_base in enumerate(ref_sequence):
        ref_base = ref_base.upper()
        if ref_base not in alt_nuc_counts:
            continue
        for nuc in alph:
            alt_nuc_counts[ref_base][nuc] += int(sub_base_vectors[ref_name + '_' + nuc][pos])
    return alt_nuc_counts


def prep_class_piechart_and_barplot(class_counts_order, class_counts,
                                    ref_names, expected_hdr_amplicon_seq,
                                    N_TOTAL):
    """Compute labels and sizes for class piechart/barplot (plot_1b/1c).

    This is a standalone utility (not CorePlotContext-based) that transforms
    raw class counts into display-ready labels and percentage sizes.

    Returns dict with 'labels', 'sizes', 'N_TOTAL'.
    """
    labels = []
    sizes = []
    for class_name in class_counts_order:
        if expected_hdr_amplicon_seq != "" and class_name == ref_names[0] + "_MODIFIED":
            labels.append("NHEJ\n(" + str(class_counts[class_name]) + " reads)")
        elif expected_hdr_amplicon_seq != "" and class_name == "HDR_MODIFIED":
            labels.append("Imperfect HDR\n(" + str(class_counts[class_name]) + " reads)")
        elif expected_hdr_amplicon_seq != "" and class_name == "HDR_UNMODIFIED":
            labels.append("HDR\n(" + str(class_counts[class_name]) + " reads)")
        else:
            display_class_name = class_name
            if len(ref_names) == 1:
                display_class_name = display_class_name.replace('Reference_', '')
            labels.append(display_class_name + "\n(" + str(class_counts[class_name]) + " reads)")
        sizes.append(100 * class_counts[class_name] / float(N_TOTAL))
    return {'labels': labels, 'sizes': sizes, 'N_TOTAL': N_TOTAL}


# =============================================================================
# Base editing utility functions
#
# These are used by both CRISPRessoCORE (which imports them from here)
# and by prep_base_edit_upset below.
# =============================================================================


def get_refpos_values(ref_aln_seq, read_aln_seq):
    """Given a reference alignment return a dict mapping reference positions to read bases.

    ``refpos_dict[ind]`` is the value of the read at the position
    corresponding to the *ind*-th base in the reference.  Any additional
    bases in the read (gaps in the ref) are assigned to the first position
    of the ref (i.e. ``refpos_dict[0]``).  For gaps in the read, the value
    is appended to the last non-gap reference position to the left.

    Example::

        ref_seq  = '--A-TGC-'
        read_seq = 'GGAGTCGA'
        get_refpos_values(ref_seq, read_seq)
        # {0: 'GGAG', 1: 'T', 2: 'C', 3: 'GA'}
    """
    refpos_dict = defaultdict(str)

    # First, if there are insertions in read, add those to the first position in ref
    if ref_aln_seq[0] == '-':
        aln_index = 0
        read_start_bases = ""
        while aln_index < len(ref_aln_seq) and ref_aln_seq[aln_index] == '-':
            read_start_bases += read_aln_seq[aln_index]
            aln_index += 1
        refpos_dict[0] = read_start_bases
        ref_aln_seq = ref_aln_seq[aln_index:]
        read_aln_seq = read_aln_seq[aln_index:]

    ref_pos = 0
    last_nongap_ref_pos = 0
    for ind in range(len(ref_aln_seq)):
        ref_base = ref_aln_seq[ind]
        read_base = read_aln_seq[ind]
        if ref_base == '-':
            refpos_dict[last_nongap_ref_pos] += read_base
        else:
            refpos_dict[ref_pos] += read_base
            last_nongap_ref_pos = ref_pos
            ref_pos += 1
    return refpos_dict


def get_bp_substitutions(ref_changes_dict, ref_seq, ref_positions_to_include):
    """Discover positions and bases that differ between reference and target (substitutions)."""
    bp_substitutions_arr = []
    for idx in ref_positions_to_include:
        ref_base = ref_seq[idx]
        if ref_changes_dict[idx] != ref_base:
            bp_substitutions_arr.append((idx, ref_base, ref_changes_dict[idx]))
    return bp_substitutions_arr


def get_upset_plot_counts(df_alleles, bp_substitutions_arr, wt_ref_name):
    """Count allele combinations for upset-style base editing plots."""
    # set up counters
    binary_allele_counts = defaultdict(int)  # e.g. T,T,X,T > 100 where each item is a string of the base at each position in bp_substitutions_arr, and 'X' is nontarget
    category_allele_counts = defaultdict(int)  # e.g. T,T,R,T > 100 where each item is a string of the base at each position in bp_substitutions_arr, and 'T' is Target, 'R' is Reference, 'D' is Deletion, 'I' is insertion, and 'N' is anything else
    precise_allele_counts = defaultdict(int)  # e.g. A,A,C,AA > 100 where each item is a string of the base at each position in bp_substitutions_arr

    total_alleles = 0
    total_alleles_reads = 0
    total_alleles_on_ref = 0
    total_alleles_reads_on_ref = 0

    total_target_noindel_reads = 0
    total_target_indel_reads = 0
    total_reference_noindel_reads = 0
    total_reference_indel_reads = 0
    total_other_noindel_reads = 0
    total_other_indel_reads = 0

    target_base_counts = [0] * len(bp_substitutions_arr)
    reference_base_counts = [0] * len(bp_substitutions_arr)
    deletion_base_counts = [0] * len(bp_substitutions_arr)
    insertion_base_counts = [0] * len(bp_substitutions_arr)
    other_base_counts = [0] * len(bp_substitutions_arr)

    # iterate all alleles in input allele table
    for idx, allele in df_alleles.iterrows():
        total_alleles += 1
        total_alleles_reads += allele['#Reads']

        if allele.Reference_Name != wt_ref_name:
            continue
        total_alleles_on_ref += 1
        total_alleles_reads_on_ref += allele['#Reads']

        has_indel_guide = False
        if allele.n_deleted > 0:
            has_indel_guide = True
        if allele.n_inserted > 0:
            has_indel_guide = True

        has_indel = has_indel_guide

        ref_aln = allele.Reference_Sequence
        read_aln = allele.Aligned_Sequence
        ref_base_position_lookup = get_refpos_values(ref_aln, read_aln)

        binary_arr = []
        cat_arr = []
        val_arr = []
        for ind, (ref_ind, ref_base, mod_base) in enumerate(bp_substitutions_arr):
            base_at_pos = ref_base_position_lookup[ref_ind]
            this_binary = 'X'
            this_category = 'N'
            if base_at_pos == ref_base:
                this_category = 'R'
                reference_base_counts[ind] += allele['#Reads']
            elif base_at_pos == mod_base:
                this_category = 'T'
                this_binary = 'T'
                target_base_counts[ind] += allele['#Reads']
            elif base_at_pos == '-':
                this_category = 'D'
                deletion_base_counts[ind] += allele['#Reads']
            elif len(base_at_pos) != 1:
                this_category = 'I'
                insertion_base_counts[ind] += allele['#Reads']
            else:
                this_category = 'N'
                other_base_counts[ind] += allele['#Reads']
            binary_arr.append(this_binary)
            cat_arr.append(this_category)
            val_arr.append(base_at_pos)

        if cat_arr.count('R') == len(cat_arr):
            if not has_indel:
                total_reference_noindel_reads += allele['#Reads']
            else:
                total_reference_indel_reads += allele['#Reads']
        elif cat_arr.count('T') == len(cat_arr):
            if not has_indel:
                total_target_noindel_reads += allele['#Reads']
            else:
                total_target_indel_reads += allele['#Reads']
        elif not has_indel:
            total_other_noindel_reads += allele['#Reads']
        else:
            total_other_indel_reads += allele['#Reads']

        binary_arr_str = "\t".join(binary_arr) + "\t" + str(has_indel)
        cat_arr_str = "\t".join(cat_arr) + "\t" + str(has_indel)
        val_arr_str = "\t".join(val_arr) + "\t" + str(has_indel)

        binary_allele_counts[binary_arr_str] += allele['#Reads']
        category_allele_counts[cat_arr_str] += allele['#Reads']
        precise_allele_counts[val_arr_str] += allele['#Reads']

    total_counts = [total_alleles_reads] * len(bp_substitutions_arr)

    return {
        "binary_allele_counts": binary_allele_counts,
        "category_allele_counts": category_allele_counts,
        "precise_allele_counts": precise_allele_counts,
        "total_alleles": total_alleles,
        "total_alleles_reads": total_alleles_reads,
        "total_alleles_on_ref": total_alleles_on_ref,
        "total_alleles_reads_on_ref": total_alleles_reads_on_ref,
        "total_target_noindel_reads": total_target_noindel_reads,
        "total_target_indel_reads": total_target_indel_reads,
        "total_reference_noindel_reads": total_reference_noindel_reads,
        "total_reference_indel_reads": total_reference_indel_reads,
        "total_other_noindel_reads": total_other_noindel_reads,
        "total_other_indel_reads": total_other_indel_reads,
        "target_base_counts": target_base_counts,
        "reference_base_counts": reference_base_counts,
        "deletion_base_counts": deletion_base_counts,
        "insertion_base_counts": insertion_base_counts,
        "other_base_counts": other_base_counts,
        "total_counts": total_counts,
    }


def get_base_edit_target_sequence(ref_seq, df_alleles, base_editor_target_ref_skip_allele_count):
    """Find the target (edited) sequence from the allele table.

    Scans *df_alleles* for the first modified allele whose unaligned
    sequence differs from *ref_seq*, skipping the first
    *base_editor_target_ref_skip_allele_count* such alleles.

    Returns the target sequence string, or ``""`` if none is found.
    """
    target_seq = ""
    seen_nonref_allele_count = 0
    for idx, allele in df_alleles.iterrows():
        if allele.Aligned_Sequence.replace("-", "") != ref_seq and allele.Read_Status == 'MODIFIED':
            if seen_nonref_allele_count >= base_editor_target_ref_skip_allele_count:
                target_seq = allele.Aligned_Sequence.replace("-", "")
                break
            else:
                logger.debug('Skipping allele ' + str(idx) + ' with sequence ' + allele.Aligned_Sequence)
            seen_nonref_allele_count += 1
    if target_seq == "":
        logger.warning('Target reference sequence not found in allele table (all reads were equal to the reference sequence)')

    return target_seq


def write_base_edit_counts(ref_name, counts_dict, bp_substitutions_arr, _jp):
    """Write base editing count tables to disk.

    Writes binary, category, precise allele counts and summary arrays
    to TSV files under the path constructed by *_jp*.
    """
    prefix = '10i.' + ref_name

    with open(_jp(prefix + '.binary_allele_counts.txt'), 'w') as fout:
        sorted_binary_allele_counts = sorted(counts_dict['binary_allele_counts'].keys(), key=lambda x: counts_dict['binary_allele_counts'][x], reverse=True)
        fout.write("\t".join([str(x) for x in bp_substitutions_arr]) + '\thas_indel\tcount\n')
        for allele_str in sorted_binary_allele_counts:
            fout.write(allele_str + '\t' + str(counts_dict['binary_allele_counts'][allele_str]) + '\n')

    with open(_jp(prefix + '.category_allele_counts.txt'), 'w') as fout:
        sorted_category_allele_counts = sorted(counts_dict['category_allele_counts'].keys(), key=lambda x: counts_dict['category_allele_counts'][x], reverse=True)
        fout.write("\t".join([str(x) for x in bp_substitutions_arr]) + '\thas_indel\tcount\n')
        for allele_str in sorted_category_allele_counts:
            fout.write(allele_str + '\t' + str(counts_dict['category_allele_counts'][allele_str]) + '\n')

    with open(_jp(prefix + '.precise_allele_counts.txt'), 'w') as fout:
        sorted_precise_allele_counts = sorted(counts_dict['precise_allele_counts'].keys(), key=lambda x: counts_dict['precise_allele_counts'][x], reverse=True)
        fout.write("\t".join([str(x) for x in bp_substitutions_arr]) + '\thas_indel\tcount\n')
        for allele_str in sorted_precise_allele_counts:
            fout.write(allele_str + '\t' + str(counts_dict['precise_allele_counts'][allele_str]) + '\n')

    with open(_jp(prefix + '.arrays.txt'), 'w') as fout:
        fout.write('Class\t' + "\t".join([str(x) for x in bp_substitutions_arr]) + '\n')
        fout.write('total_counts\t' + "\t".join([str(x) for x in counts_dict['total_counts']]) + '\n')
        fout.write('reference_counts\t' + "\t".join([str(x) for x in counts_dict['reference_base_counts']]) + '\n')
        fout.write('target_counts\t' + "\t".join([str(x) for x in counts_dict['target_base_counts']]) + '\n')
        fout.write('deletion_counts\t' + "\t".join([str(x) for x in counts_dict['deletion_base_counts']]) + '\n')
        fout.write('insertion_counts\t' + "\t".join([str(x) for x in counts_dict['insertion_base_counts']]) + '\n')
        fout.write('other_counts\t' + "\t".join([str(x) for x in counts_dict['other_base_counts']]) + '\n')

    with open(_jp(prefix + '.counts.txt'), 'w') as fout:
        target_name = 'Target'
        fout.write("\t".join([ref_name, ref_name + "_indels", target_name, target_name + "_indels", "other", "other_indels"]) + '\n')
        fout.write("\t".join([str(x) for x in [counts_dict['total_reference_noindel_reads'], counts_dict['total_reference_indel_reads'], counts_dict['total_target_noindel_reads'], counts_dict['total_target_indel_reads'], counts_dict['total_other_noindel_reads'], counts_dict['total_other_indel_reads']]]) + '\n')


def prep_base_edit_upset(ref_seq, df_alleles, ref_name, sgRNA_interval,
                         args, gap_incentive=None, sgRNA_legend=''):
    """Prepare data for base edit upset plot (plot_10i).

    Finds the target sequence, aligns it to reference, identifies
    base changes, and computes combination counts.

    This is a standalone utility (not CorePlotContext-based) that encapsulates
    the inline logic from CRISPRessoCORE for plot_10i.

    Parameters
    ----------
    ref_seq : str
        Reference amplicon sequence.
    df_alleles : pd.DataFrame
        Allele table (full, not filtered by ref_name).
    ref_name : str
        Name of the reference (wild-type) amplicon.
    sgRNA_interval : tuple
        (start, end) of the sgRNA interval for this sgRNA.
    args : argparse.Namespace
        Must have: base_editor_target_ref_skip_allele_count,
        base_editor_consider_changes_outside_qw,
        needleman_wunsch_gap_open, needleman_wunsch_gap_extend,
        needleman_wunsch_aln_matrix_loc, quantification_window_coordinates.
    gap_incentive : np.ndarray or None
        Per-position gap incentive array for alignment.  If ``None``,
        zeros are used.

    Returns
    -------
    dict or None
        Dict with ``'bp_substitutions_arr'`` and ``'counts_dict'``
        (the full output of :func:`get_upset_plot_counts`),
        or ``None`` if no target sequence found or
        ``quantification_window_coordinates`` is set.

    """
    target_seq = get_base_edit_target_sequence(
        ref_seq, df_alleles,
        args.base_editor_target_ref_skip_allele_count,
    )

    if not target_seq:
        return None
    if args.quantification_window_coordinates is not None:
        return None

    # Lazy import — CRISPResso2Align is a compiled C extension
    from CRISPResso2 import CRISPResso2Align

    # Align target to reference
    aln_matrix_loc = args.needleman_wunsch_aln_matrix_loc
    if aln_matrix_loc == 'EDNAFULL':
        aln_matrix = CRISPResso2Align.make_matrix()
    else:
        if not os.path.exists(aln_matrix_loc):
            return None
        aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    if gap_incentive is None:
        gap_incentive = np.zeros(len(ref_seq) + 1, dtype=int)

    aln_target_seq, aln_ref_seq, _aln_score = CRISPResso2Align.global_align(
        target_seq,
        ref_seq,
        matrix=aln_matrix,
        gap_incentive=gap_incentive,
        gap_open=args.needleman_wunsch_gap_open,
        gap_extend=args.needleman_wunsch_gap_extend,
    )

    # Determine which positions to include
    if args.base_editor_consider_changes_outside_qw:
        ref_positions_to_include = list(range(len(ref_seq)))
    else:
        this_start, this_stop = sgRNA_interval
        ref_positions_to_include = list(range(this_start, this_stop + 1))

    ref_changes_dict = get_refpos_values(aln_ref_seq, aln_target_seq)
    bp_substitutions_arr = get_bp_substitutions(
        ref_changes_dict, ref_seq, ref_positions_to_include,
    )

    counts_dict = get_upset_plot_counts(
        df_alleles, bp_substitutions_arr, ref_name,
    )

    conversion_nuc_from = getattr(args, 'conversion_nuc_from', 'C')
    caption = (
        f"Figure 10i: Upset plot of Base Edits for {conversion_nuc_from} "
        f"around cut site for {sgRNA_legend}. Each dot matrix at the bottom "
        f"represents a specific combination of base edits (colored by target "
        f"position), and the bar plot at the top shows the number of reads "
        f"with each combination."
    ) if sgRNA_legend else ''

    return {
        'bp_substitutions_arr': bp_substitutions_arr,
        'counts_dict': counts_dict,
        'caption': caption,
        'data_files': [],
    }


def prep_conversion_at_sel_nucs(ctx: CorePlotContext):
    """Compute from_nuc_indices and selected-nucleotide DataFrames for plots 10e/10f/10g.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Finds the column positions matching ``conversion_nuc_from`` and returns
    sliced, re-labeled DataFrames used for both CSV writes and plot inputs.

    Returns a dict with:
    - ``from_nuc_indices``: list of int column positions matching the conversion nucleotide
    - ``just_sel_nuc_pcts``: DataFrame sliced to those columns, columns renamed
    - ``just_sel_nuc_freqs``: same slicing applied to frequency DataFrame
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    sgRNA_ind = ctx.sgRNA_ind
    conversion_nuc_from = ctx.args.conversion_nuc_from

    plot_idxs = ref['sgRNA_plot_idxs'][sgRNA_ind]

    # Build nucleotide freq/pct DataFrames and slice to plot window
    df_nuc_freq_all, df_nuc_pct_all = _build_nuc_freq_df(ctx)
    plot_nuc_pcts = df_nuc_pct_all.iloc[:, plot_idxs]
    plot_nuc_freqs = df_nuc_freq_all.iloc[:, plot_idxs]

    from_nuc_indices = [
        pos for pos, char in enumerate(list(plot_nuc_pcts.columns.values))
        if char == conversion_nuc_from
    ]

    just_sel_nuc_pcts = plot_nuc_pcts.iloc[:, from_nuc_indices].copy()
    just_sel_nuc_pcts.columns = [
        conversion_nuc_from + str(pos + 1) for pos in from_nuc_indices
    ]

    just_sel_nuc_freqs = plot_nuc_freqs.iloc[:, from_nuc_indices].copy()
    just_sel_nuc_freqs.columns = [
        conversion_nuc_from + str(pos + 1) for pos in from_nuc_indices
    ]

    return {
        'from_nuc_indices': from_nuc_indices,
        'just_sel_nuc_pcts': just_sel_nuc_pcts,
        'just_sel_nuc_freqs': just_sel_nuc_freqs,
    }


# =============================================================================
# Remaining prep functions — packaging for uniform prep → plot → metadata shape
# =============================================================================


def prep_read_barplot(ctx: CorePlotContext):
    """Prepare kwargs for plot_read_barplot (plot_1a)."""
    aln_stats = ctx.run_data['running_info']['alignment_stats']
    fig_filename_root = _make_fig_filename_root(ctx, '1a.Read_barplot')
    mapping_stats_filename = _get_run_info(ctx, 'mapping_stats_filename', '')
    return {
        'N_READS_INPUT': aln_stats['N_READS_INPUT'],
        'N_READS_AFTER_PREPROCESSING': aln_stats['N_READS_AFTER_PREPROCESSING'],
        'N_TOTAL': ctx.N_TOTAL,
        'fig_filename_root': fig_filename_root,
        'save_png': ctx.save_png,
        'caption': "Figure 1a: The number of reads in input fastqs, after preprocessing, and after alignment to amplicons.",
        'data_files': [('Mapping statistics', os.path.basename(mapping_stats_filename))] if mapping_stats_filename else [],
    }


def prep_class_piechart_and_barplot_plot(ctx: CorePlotContext):
    """Prepare kwargs for plot_class_piechart_and_barplot (plot_1b/1c).

    Returns captions for both plots: ``piechart_caption`` and ``barplot_caption``.
    The piechart caption has an HDR-specific variant when ``expected_hdr_amplicon_seq`` is set.
    """
    piechart_plot_root = _make_fig_filename_root(ctx, '1b.Alignment_pie_chart')
    barplot_plot_root = _make_fig_filename_root(ctx, '1c.Alignment_barplot')
    quant_of_editing_freq_filename = _get_run_info(ctx, 'quant_of_editing_freq_filename', '')

    piechart_caption = (
        "Figure 1b: Alignment and editing frequency of reads as determined by the "
        "percentage and number of sequence reads showing unmodified and modified alleles."
    )
    if ctx.args.expected_hdr_amplicon_seq != "":
        piechart_caption = (
            "Figure 1b: Alignment and editing frequency of reads as determined by the "
            "percentage and number of sequence reads showing unmodified and modified alleles. "
            "NHEJ reads align more closely to the unmodified reference sequence, but have mutations "
            "present in the specified quantification window. HDR reads align to the HDR reference "
            "sequence and have no mutations in the specified quantification window. Imperfect HDR "
            "reads have mutations in the specified window. AMBIGUOUS reads align equally well to "
            "the unmodified and HDR reference sequences."
        )

    barplot_caption = (
        "Figure 1c: Alignment and editing frequency of reads as determined by the "
        "percentage and number of sequence reads showing unmodified and modified alleles."
    )

    data_files = [('Quantification of editing', os.path.basename(quant_of_editing_freq_filename))] if quant_of_editing_freq_filename else []

    return {
        'class_counts_order': ctx.class_counts_order,
        'class_counts': ctx.class_counts,
        'ref_names': ctx.ref_names,
        'expected_hdr_amplicon_seq': ctx.args.expected_hdr_amplicon_seq,
        'N_TOTAL': ctx.N_TOTAL,
        'piechart_plot_root': piechart_plot_root,
        'barplot_plot_root': barplot_plot_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_png': ctx.save_png,
        'piechart_caption': piechart_caption,
        'barplot_caption': barplot_caption,
        'piechart_data_files': data_files,
        'barplot_data_files': data_files,
    }


def prep_alleles_homology_histogram(ctx: CorePlotContext):
    """Prepare kwargs for plot_alleles_homology_histogram (plot_1e)."""
    fig_filename_root = _make_fig_filename_root(ctx, '1e.Allele_homology_histogram')
    return {
        'fig_root': fig_filename_root,
        'homology_scores': ctx.homology_scores,
        'counts': ctx.homology_counts,
        'min_homology': ctx.args.default_min_aln_score,
        'save_also_png': ctx.save_png,
        'caption': (
            "Figure 1e: Distribution of read alignment homology scores, showing the "
            "best-scoring alignment of each sequencing read to the provided amplicons. "
            "The dashed line indicates the minimum alignment score threshold used to "
            "discard low-quality alignments."
        ),
        'data_files': [('Alleles Homology Scores', _data_file_basename(ctx, 'Alleles_homology_scores.txt'))],
    }


def prep_quantification_window_locations(ctx: CorePlotContext):
    """Prepare kwargs for plot_quantification_window_locations (plot_4c).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '4c.' + _ref_plot_name(ctx) + 'Quantification_window_insertion_deletion_substitution_locations',
    )

    return {
        'insertion_count_vectors': ctx.insertion_count_vectors[ref_name],
        'deletion_count_vectors': ctx.deletion_count_vectors[ref_name],
        'substitution_count_vectors': ctx.substitution_count_vectors[ref_name],
        'include_idxs_list': ref['include_idxs'],
        'cut_points': ref['sgRNA_cut_points'],
        'plot_cut_points': ref['sgRNA_plot_cut_points'],
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'ref_len': ref['sequence_length'],
        'num_refs': num_refs,
        'n_total': ctx.N_TOTAL,
        'n_this_category': ctx.counts_total[ref_name],
        'plot_title': plot_title_with_ref_name(
            'Mutation position distribution', ref_name, num_refs,
        ),
        'ref_name': ref_name,
        'fig_filename_root': fig_filename_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
        'caption': "Figure 4c: Frequency of insertions, deletions, and substitutions across the entire amplicon, considering only modifications that overlap with the quantification window.",
        'data_files': [('Modification frequency in quantification window', os.path.basename(_get_ref_info(ctx, ref_name, 'quant_window_mod_count_filename', '')))] if _get_ref_info(ctx, ref_name, 'quant_window_mod_count_filename', '') else [],
    }


def prep_position_dependent_indels(ctx: CorePlotContext):
    """Prepare kwargs for plot_position_dependent_indels (plot_4d).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '4d.' + _ref_plot_name(ctx) + 'Position_dependent_average_indel_size',
    )

    return {
        'insertion_length_vectors': ctx.insertion_length_vectors[ref_name],
        'deletion_length_vectors': ctx.deletion_length_vectors[ref_name],
        'cut_points': ref['sgRNA_cut_points'],
        'plot_cut_points': ref['sgRNA_plot_cut_points'],
        'ref_len': ref['sequence_length'],
        'plot_titles': {
            'ins': plot_title_with_ref_name(
                'Position dependent insertion size', ref_name, num_refs,
            ),
            'del': plot_title_with_ref_name(
                'Position dependent deletion size', ref_name, num_refs,
            ),
        },
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'ref_name': ref_name,
        'caption': "Figure 4d: Position dependent insertion size(left) and deletion size (right), including only modifications that overlap with the quantification window.",
        'data_files': [],
    }


def prep_frameshift_analysis(ctx: CorePlotContext):
    """Prepare kwargs for plot_frameshift_analysis (plot_5).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '5.' + _ref_plot_name(ctx) + 'Frameshift_in-frame_mutations_pie_chart',
    )

    return {
        'modified_frameshift': ctx.counts_modified_frameshift[ref_name],
        'modified_non_frameshift': ctx.counts_modified_non_frameshift[ref_name],
        'non_modified_non_frameshift': ctx.counts_non_modified_non_frameshift[ref_name],
        'cut_points': ref['sgRNA_cut_points'],
        'plot_cut_points': ref['sgRNA_plot_cut_points'],
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'exon_intervals': ref['exon_intervals'],
        'ref_len': ref['sequence_length'],
        'ref_name': ref_name,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': "Figure 5: Frameshift analysis of coding sequence reads affected by modifications (unmodified reads are excluded from this analysis).",
        'data_files': [('Frameshift analysis', os.path.basename(_get_ref_info(ctx, ref_name, 'frameshift_analysis_filename', '')))] if _get_ref_info(ctx, ref_name, 'frameshift_analysis_filename', '') else [],
    }


def prep_frameshift_frequency(ctx: CorePlotContext):
    """Prepare kwargs for plot_frameshift_frequency (plot_6).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '6.' + _ref_plot_name(ctx) + 'Frameshift_in-frame_mutation_profiles',
    )

    hists_inframe_count = ctx.hists_inframe[ref_name][0]
    return {
        'hists_frameshift': ctx.hists_frameshift[ref_name],
        'hists_inframe': ctx.hists_inframe[ref_name],
        'plot_titles': {
            'fs': plot_title_with_ref_name('Frameshift profile', ref_name, num_refs),
            'if': plot_title_with_ref_name('In-frame profile', ref_name, num_refs),
        },
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'ref_name': ref_name,
        'caption': (
            "Figure 6: Frameshift and in-frame mutagenesis profiles indicating position "
            "affected by modification. The y axis shows the number of reads and percentage "
            "of all reads in that category (frameshifted (top) or in-frame (bottom)). "
            "%d reads with no length modifications are not shown." % hists_inframe_count
        ),
        'data_files': [],
    }


def prep_non_coding_mutations(ctx: CorePlotContext):
    """Prepare kwargs for plot_non_coding_mutations (plot_7).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '7.' + _ref_plot_name(ctx) + 'Insertion_deletion_substitution_locations_noncoding',
    )

    return {
        'insertion_count_vectors_noncoding': ctx.insertion_count_vectors_noncoding[ref_name],
        'deletion_count_vectors_noncoding': ctx.deletion_count_vectors_noncoding[ref_name],
        'substitution_count_vectors_noncoding': ctx.substitution_count_vectors_noncoding[ref_name],
        'include_idxs_list': ref['include_idxs'],
        'cut_points': ref['sgRNA_cut_points'],
        'plot_cut_points': ref['sgRNA_plot_cut_points'],
        'ref_len': ref['sequence_length'],
        'sgRNA_intervals': ref['sgRNA_intervals'],
        'plot_title': plot_title_with_ref_name(
            'Noncoding mutation position distribution', ref_name, num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'save_also_png': ctx.save_png,
        'ref_name': ref_name,
        'caption': "Figure 7: Reads with insertions, deletions, and substitutions mapped to reference amplicon position exclusively in noncoding region/s (that is, without mutations affecting coding sequences). The predicted cleavage site is indicated by a vertical dashed line. Only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are plotted.",
        'data_files': [],
    }


def prep_potential_splice_sites(ctx: CorePlotContext):
    """Prepare kwargs for plot_potential_splice_sites (plot_8).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '8.' + _ref_plot_name(ctx) + 'Potential_splice_sites_pie_chart',
    )

    return {
        'splicing_sites_modified': ctx.counts_splicing_sites_modified[ref_name],
        'count_total': ctx.counts_total[ref_name],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'ref_name': ref_name,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': "Figure 8: Predicted impact on splice sites. Potential splice sites modified refers to reads in which the either of the two intronic positions adjacent to exon junctions are disrupted.",
        'data_files': [],
    }


def prep_subs_across_ref(ctx: CorePlotContext):
    """Prepare kwargs for plot_subs_across_ref (plot_10a).

    Requires ``ctx.ref_name``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '10a.' + _ref_plot_name(ctx) + 'Substitution_frequencies_at_each_bp',
    )

    return {
        'ref_len': ref['sequence_length'],
        'ref_seq': ref['sequence'],
        'ref_name': ref_name,
        'ref_count': ctx.counts_total[ref_name],
        'all_substitution_base_vectors': ctx.all_substitution_base_vectors,
        'plot_title': plot_title_with_ref_name(
            'Substitution frequency', ref_name, num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'quantification_window_idxs': ref['include_idxs'],
        'custom_colors': ctx.custom_config.get('colors', {}),
        'ref_name': ref_name,
        'caption': "Figure 10a: Substitution frequencies across the amplicon.",
        'data_files': [],
    }


def prep_sub_freq_barplot(ctx: CorePlotContext):
    """Prepare kwargs for plot_sub_freqs — full amplicon (plot_10b).

    Requires ``ctx.ref_name``. Uses ``ctx.alt_nuc_counts_all``.
    """
    ref_name = ctx.ref_name
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '10b.' + _ref_plot_name(ctx) + 'Substitution_frequency_barplot',
    )

    return {
        'alt_nuc_counts': ctx.alt_nuc_counts_all[ref_name],
        'plot_title': plot_title_with_ref_name(
            'Substitution frequency\nin entire amplicon', ref_name, num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': "Figure 10b: Substitution frequencies across the amplicon.",
        'data_files': [('Nucleotide frequencies', os.path.basename(_get_ref_info(ctx, ref_name, 'nuc_freq_filename', '')))] if _get_ref_info(ctx, ref_name, 'nuc_freq_filename', '') else [],
    }


def prep_sub_freq_barplot_quant_window(ctx: CorePlotContext):
    """Prepare kwargs for plot_sub_freqs — quantification window (plot_10c).

    Requires ``ctx.ref_name``. Uses ``ctx.alt_nuc_counts``.
    """
    ref_name = ctx.ref_name
    num_refs = len(ctx.ref_names)

    fig_filename_root = _make_fig_filename_root(
        ctx,
        '10c.' + _ref_plot_name(ctx) + 'Substitution_frequency_barplot_in_quantification_window',
    )

    return {
        'alt_nuc_counts': ctx.alt_nuc_counts[ref_name],
        'plot_title': plot_title_with_ref_name(
            'Substitution frequency\nin quantification window', ref_name, num_refs,
        ),
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': "Figure 10c: Substitution frequencies in the quantification window",
        'data_files': [('Nucleotide frequencies in quantification window', os.path.basename(_get_ref_info(ctx, ref_name, 'quant_window_sub_freq_filename', '')))] if _get_ref_info(ctx, ref_name, 'quant_window_sub_freq_filename', '') else [],
    }


def _prep_conversion_at_sel_nucs_common(ctx: CorePlotContext, plot_number: str, variant: str = ''):
    """Shared logic for plots 10e, 10f, and 10g.

    Builds the nucleotide percentage DataFrame sliced to the sgRNA plot
    window and returns the common kwargs dict.  The three public
    ``prep_conversion_at_sel_nucs_*`` functions differ only in the
    *plot_number* (e.g. ``'10e'``) and optional *variant* infix
    (e.g. ``'no_ref_'``) they pass here.

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.
    """
    ref_name = ctx.ref_name
    ref = _ref(ctx)
    num_refs = len(ctx.ref_names)
    plot_idxs = ref['sgRNA_plot_idxs'][ctx.sgRNA_ind]

    _df_nuc_freq_all, df_nuc_pct_all = _build_nuc_freq_df(ctx)
    plot_nuc_pcts = df_nuc_pct_all.iloc[:, plot_idxs]
    ref_seq_slice = ''.join([ref['sequence'][i] for i in plot_idxs])

    fig_filename_root = _make_fig_filename_root(
        ctx,
        plot_number + '.' + _ref_plot_name(ctx) + 'Selected_conversion_' + variant + 'at_'
        + ctx.args.conversion_nuc_from + 's_around_' + _sgRNA_label(ctx),
    )

    return {
        'df_subs': plot_nuc_pcts,
        'ref_name': ref_name,
        'ref_sequence': ref_seq_slice,
        'plot_title': plot_title_with_ref_name(
            'Substitution Frequencies at ' + ctx.args.conversion_nuc_from
            + 's around the ' + _sgRNA_legend(ctx),
            ref_name, num_refs,
        ),
        'conversion_nuc_from': ctx.args.conversion_nuc_from,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
    }


def _sel_nuc_freq_data_files(ctx: CorePlotContext):
    """Build data_files for plots 10e/10f/10g (selected nucleotide frequency table)."""
    sel_nuc_freq_basename = _data_file_basename(ctx, _ref_plot_name(ctx) + 'Selected_nucleotide_frequency_table_around_' + _sgRNA_label(ctx) + '.txt')
    return [('Nucleotide frequencies at ' + ctx.args.conversion_nuc_from + 's', sel_nuc_freq_basename)]


def prep_conversion_at_sel_nucs_plot(ctx: CorePlotContext):
    """Prepare kwargs for plot_conversion_at_sel_nucs (plot_10e).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.

    Builds nucleotide percentage DataFrame sliced to the sgRNA plot window.
    """
    result = _prep_conversion_at_sel_nucs_common(ctx, '10e')
    sgRNA_legend = _sgRNA_legend(ctx)
    result['caption'] = (
        "Figure 10e: Proportion of each base at each nucleotide targeted by base editors "
        "in the plotting window around the " + sgRNA_legend + ". The number of each target "
        "base is annotated on the reference sequence at the bottom of the plot."
    )
    result['data_files'] = _sel_nuc_freq_data_files(ctx)
    return result


def prep_conversion_at_sel_nucs_not_include_ref(ctx: CorePlotContext):
    """Prepare kwargs for plot_conversion_at_sel_nucs_not_include_ref (plot_10f).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.
    """
    result = _prep_conversion_at_sel_nucs_common(ctx, '10f', variant='no_ref_')
    sgRNA_legend = _sgRNA_legend(ctx)
    result['caption'] = (
        "Figure 10f: Non-reference base proportions. For target nucleotides in the plotting "
        "window, this plot shows the proportion of non-reference (non-"
        + ctx.args.conversion_nuc_from + ") bases as a percentage of all non-reference "
        "sequences. The number of each target base is annotated on the reference sequence "
        "at the bottom of the plot."
    )
    result['data_files'] = _sel_nuc_freq_data_files(ctx)
    return result


def prep_conversion_at_sel_nucs_not_include_ref_scaled(ctx: CorePlotContext):
    """Prepare kwargs for plot_conversion_at_sel_nucs_not_include_ref_scaled (plot_10g).

    Requires ``ctx.ref_name`` and ``ctx.sgRNA_ind``.
    """
    result = _prep_conversion_at_sel_nucs_common(ctx, '10g', variant='no_ref_scaled_')
    sgRNA_legend = _sgRNA_legend(ctx)
    result['caption'] = (
        "Figure 10g: Non-reference base counts. For target nucleotides in the plotting "
        "window, this plot shows the number of non-reference (non-"
        + ctx.args.conversion_nuc_from + ") bases. The number of each target base is "
        "annotated on the reference sequence at the bottom of the plot."
    )
    result['data_files'] = _sel_nuc_freq_data_files(ctx)
    return result


def prep_global_frameshift_analysis(ctx: CorePlotContext, global_data: dict | None = None):
    """Prepare kwargs for plot_global_frameshift_analysis (plot_5a).

    *global_data* may be passed from a prior
    :func:`prep_global_frameshift_data` call to avoid recomputing it.
    If ``None``, it is computed here automatically.
    """
    if global_data is None:
        global_data = prep_global_frameshift_data(ctx)

    fig_filename_root = _make_fig_filename_root(
        ctx, '5a.Global_frameshift_in-frame_mutations_pie_chart',
    )

    return {
        'global_modified_frameshift': global_data['global_modified_frameshift'],
        'global_modified_non_frameshift': global_data['global_modified_non_frameshift'],
        'global_non_modified_non_frameshift': global_data['global_non_modified_non_frameshift'],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': "Figure 5a: Frameshift analysis of coding sequence reads affected by modifications for all reads. Unmodified reference reads are excluded from this plot, and all HDR reads are included in this plot.",
        'data_files': [],
    }


def prep_global_frameshift_in_frame_mutations(ctx: CorePlotContext, global_data: dict | None = None):
    """Prepare kwargs for plot_global_frameshift_in_frame_mutations (plot_6a).

    *global_data* may be passed from a prior
    :func:`prep_global_frameshift_data` call to avoid recomputing it.
    If ``None``, it is computed here automatically.
    """
    if global_data is None:
        global_data = prep_global_frameshift_data(ctx)

    fig_filename_root = _make_fig_filename_root(
        ctx, '6a.Global_frameshift_in-frame_mutation_profiles',
    )

    global_hists_inframe_count = global_data['global_hists_inframe'][0]
    return {
        'global_hists_frameshift': global_data['global_hists_frameshift'],
        'global_hists_inframe': global_data['global_hists_inframe'],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'caption': (
            "Figure 6a: Frameshift and in-frame mutagenesis profiles for all reads indicating "
            "position affected by modification. The y axis shows the number of reads and "
            "percentage of all reads in that category (frameshifted (top) or in-frame (bottom)). "
            "%d reads with no length modifications are not shown." % global_hists_inframe_count
        ),
        'data_files': [],
    }


def prep_impact_on_splice_sites(ctx: CorePlotContext, global_data: dict | None = None):
    """Prepare kwargs for plot_impact_on_splice_sites (plot_8a).

    *global_data* may be passed from a prior
    :func:`prep_global_frameshift_data` call to avoid recomputing it.
    If ``None``, it is computed here automatically.
    """
    if global_data is None:
        global_data = prep_global_frameshift_data(ctx)

    fig_filename_root = _make_fig_filename_root(
        ctx, '8a.Global_potential_splice_sites_pie_chart',
    )

    return {
        'global_splicing_sites_modified': global_data['global_splicing_sites_modified'],
        'global_count_total': global_data['global_count_total'],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'caption': "Figure 8a: Predicted impact on splice sites for all reads. Potential splice sites modified refers to reads in which the either of the two intronic positions adjacent to exon junctions are disrupted.",
        'data_files': [],
    }


def prep_scaffold_indel_lengths(ctx: CorePlotContext):
    """Prepare kwargs for plot_scaffold_indel_lengths (plot_11c)."""
    fig_filename_root = _make_fig_filename_root(
        ctx, '11c.Prime_editing_scaffold_insertion_sizes',
    )

    scaffold_insertion_sizes_filename = _get_run_info(ctx, 'scaffold_insertion_sizes_filename', '')
    return {
        'df_scaffold_insertion_sizes': ctx.df_scaffold_insertion_sizes,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'caption': (
            "Figure 11c: Scaffold insertion lengths and deletion lengths in reads that contain "
            "a scaffold insertion. 'Length matching scaffold' shows the number of basepairs "
            "immediately after the pegRNA extension sequence that exactly match the scaffold RNA "
            "sequence. 'Insertion length' shows the length of the insertion immediately after "
            "the pegRNA extension sequence (including bases that do not match the scaffold sequence)."
        ),
        'data_files': [('Scaffold insertion alleles with insertion sizes', os.path.basename(scaffold_insertion_sizes_filename))] if scaffold_insertion_sizes_filename else [],
    }


# =============================================================================
# Batch / Aggregate prep functions
# =============================================================================


def _compute_sub_sgRNA_intervals(sgRNA_plot_idxs, consensus_sgRNA_intervals, include_idxs):
    """Map full-amplicon sgRNA intervals to sub-indexed coordinates.

    When plotting around a single sgRNA, the DataFrame columns are a subset
    of the full amplicon. This function computes where the sgRNA intervals
    fall in the sub-indexed coordinate system.

    Returns (sub_sgRNA_intervals, sub_include_idxs).
    """
    sub_sgRNA_intervals = []
    for sgRNA_interval in consensus_sgRNA_intervals:
        newstart = None
        newend = None
        for idx, i in enumerate(sgRNA_plot_idxs):
            if i <= sgRNA_interval[0]:
                newstart = idx
            if newend is None and i >= sgRNA_interval[1]:
                newend = idx
        # if guide doesn't overlap with plot idxs
        if newend == 0 or newstart == len(sgRNA_plot_idxs):
            continue
        # correct partial overlaps
        elif newstart is None and newend is None:
            newstart = 0
            newend = len(include_idxs) - 1
        elif newstart is None:
            newstart = 0
        elif newend is None:
            newend = len(include_idxs) - 1
        sub_sgRNA_intervals.append((newstart, newend))
    # Scale include_idxs to sub-indexed coordinates
    sub_include_idxs = include_idxs - sgRNA_plot_idxs[0]
    return sub_sgRNA_intervals, sub_include_idxs


def prep_batch_nuc_quilt_around_sgRNA(ctx):
    """Prepare kwargs for plot_nucleotide_quilt around one sgRNA (Batch/Aggregate).

    Parameters
    ----------
    ctx : BatchPlotContext | AggregatePlotContext
        Must have ``amplicon_name`` and ``sgRNA_ind`` scope fields set.

    Returns
    -------
    dict
        kwargs for ``CRISPRessoPlot.plot_nucleotide_quilt``, plus
        ``sub_consensus_guides`` for metadata use.

    """
    amp = ctx.amplicon_name
    sgRNA_ind = ctx.sgRNA_ind
    sgRNA_plot_idxs = ctx.consensus_sgRNA_plot_idxs[amp][sgRNA_ind]
    include_idxs = ctx.consensus_include_idxs[amp]

    # Slice DataFrames to sgRNA columns
    plot_idxs_flat = [0, 1]  # Batch, Nucleotide/Modification columns
    plot_idxs_flat.extend([plot_idx + 2 for plot_idx in sgRNA_plot_idxs])

    nuc_pct_df = ctx.nucleotide_percentage_summary_dfs[amp].iloc[:, plot_idxs_flat]
    nuc_pct_df = pd.concat([nuc_pct_df.iloc[:, 0:2], nuc_pct_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)

    mod_pct_df = ctx.modification_percentage_summary_dfs[amp].iloc[:, plot_idxs_flat]
    mod_pct_df = pd.concat([mod_pct_df.iloc[:, 0:2], mod_pct_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)

    # Compute sub-indexed sgRNA intervals
    sub_sgRNA_intervals, sub_include_idxs = _compute_sub_sgRNA_intervals(
        sgRNA_plot_idxs, ctx.consensus_sgRNA_intervals[amp], include_idxs,
    )

    # Collect sub_consensus_guides for metadata
    sub_consensus_guides = []
    for sgRNA_index, sgRNA_interval in enumerate(ctx.consensus_sgRNA_intervals[amp]):
        for idx_val, i in enumerate(sgRNA_plot_idxs):
            if i >= sgRNA_interval[0]:
                sub_consensus_guides.append(ctx.consensus_guides[amp][sgRNA_index])
                break

    sgRNA = ctx.consensus_guides[amp][sgRNA_ind]
    amplicon_plot_name = amp + "."
    if len(ctx.amplicon_names) == 1 and amp == "Reference":
        amplicon_plot_name = ""
    fig_filename_root = ctx._jp(amplicon_plot_name.replace('.', '') + 'Nucleotide_percentage_quilt_around_sgRNA_' + sgRNA)

    # Detect group column name from DataFrame (Batch='Batch', Aggregate='Folder')
    group_column = ctx.nucleotide_percentage_summary_dfs[amp].columns[0]

    return {
        'nuc_pct_df': nuc_pct_df,
        'mod_pct_df': mod_pct_df,
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': sub_sgRNA_intervals,
        'sgRNA_sequences': sub_consensus_guides,
        'quantification_window_idxs': sub_include_idxs,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'group_column': group_column,
    }


def prep_batch_nuc_quilt(ctx):
    """Prepare kwargs for plot_nucleotide_quilt for full region (Batch/Aggregate).

    Parameters
    ----------
    ctx : BatchPlotContext | AggregatePlotContext
        Must have ``amplicon_name`` scope field set.

    """
    amp = ctx.amplicon_name
    amplicon_plot_name = amp + "."
    if len(ctx.amplicon_names) == 1 and amp == "Reference":
        amplicon_plot_name = ""
    fig_filename_root = ctx._jp(amplicon_plot_name.replace('.', '') + 'Nucleotide_percentage_quilt')

    # Detect group column name from DataFrame (Batch='Batch', Aggregate='Folder')
    group_column = ctx.nucleotide_percentage_summary_dfs[amp].columns[0]

    result = {
        'nuc_pct_df': ctx.nucleotide_percentage_summary_dfs[amp],
        'mod_pct_df': ctx.modification_percentage_summary_dfs[amp],
        'fig_filename_root': fig_filename_root,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
        'group_column': group_column,
    }
    if ctx.guides_all_same[amp] and ctx.consensus_guides[amp]:
        result['sgRNA_intervals'] = ctx.consensus_sgRNA_intervals[amp]
        result['sgRNA_sequences'] = ctx.consensus_guides[amp]
        result['quantification_window_idxs'] = ctx.consensus_include_idxs[amp]
    return result


def prep_batch_conversion_map_around_sgRNA(ctx):
    """Prepare kwargs for plot_conversion_map around one sgRNA (Batch only).

    Parameters
    ----------
    ctx : BatchPlotContext
        Must have ``amplicon_name`` and ``sgRNA_ind`` scope fields set.

    """
    amp = ctx.amplicon_name
    sgRNA_ind = ctx.sgRNA_ind
    sgRNA_plot_idxs = ctx.consensus_sgRNA_plot_idxs[amp][sgRNA_ind]
    include_idxs = ctx.consensus_include_idxs[amp]

    plot_idxs_flat = [0, 1]
    plot_idxs_flat.extend([plot_idx + 2 for plot_idx in sgRNA_plot_idxs])

    nuc_pct_df = ctx.nucleotide_percentage_summary_dfs[amp].iloc[:, plot_idxs_flat]
    nuc_pct_df = pd.concat([nuc_pct_df.iloc[:, 0:2], nuc_pct_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)

    sub_sgRNA_intervals, sub_include_idxs = _compute_sub_sgRNA_intervals(
        sgRNA_plot_idxs, ctx.consensus_sgRNA_intervals[amp], include_idxs,
    )

    sgRNA = ctx.consensus_guides[amp][sgRNA_ind]
    amplicon_plot_name = amp + "."
    if len(ctx.amplicon_names) == 1 and amp == "Reference":
        amplicon_plot_name = ""
    fig_filename_root = ctx._jp(amplicon_plot_name + 'Nucleotide_conversion_map_around_sgRNA_' + sgRNA)

    return {
        'nuc_pct_df': nuc_pct_df,
        'fig_filename_root': fig_filename_root,
        'conversion_nuc_from': ctx.args.conversion_nuc_from,
        'conversion_nuc_to': ctx.args.conversion_nuc_to,
        'save_also_png': ctx.save_png,
        'sgRNA_intervals': sub_sgRNA_intervals,
        'quantification_window_idxs': sub_include_idxs,
        'custom_colors': ctx.custom_config.get('colors', {}),
    }


def prep_batch_conversion_map(ctx):
    """Prepare kwargs for plot_conversion_map for full region (Batch only).

    Parameters
    ----------
    ctx : BatchPlotContext
        Must have ``amplicon_name`` scope field set.

    """
    amp = ctx.amplicon_name
    amplicon_plot_name = amp + "."
    if len(ctx.amplicon_names) == 1 and amp == "Reference":
        amplicon_plot_name = ""
    fig_filename_root = ctx._jp(amplicon_plot_name + 'Nucleotide_conversion_map')

    result = {
        'nuc_pct_df': ctx.nucleotide_percentage_summary_dfs[amp],
        'fig_filename_root': fig_filename_root,
        'conversion_nuc_from': ctx.args.conversion_nuc_from,
        'conversion_nuc_to': ctx.args.conversion_nuc_to,
        'save_also_png': ctx.save_png,
        'custom_colors': ctx.custom_config.get('colors', {}),
    }
    if ctx.guides_all_same[amp] and ctx.consensus_guides[amp]:
        result['sgRNA_intervals'] = ctx.consensus_sgRNA_intervals[amp]
        result['quantification_window_idxs'] = ctx.consensus_include_idxs[amp]
    return result


def prep_batch_allele_modification_heatmap(ctx):
    """Prepare kwargs for plot_allele_modification_heatmap (Batch/Aggregate).

    Parameters
    ----------
    ctx : BatchPlotContext | AggregatePlotContext
        Must have ``amplicon_name`` and ``mod_type`` scope fields set.

    """
    amp = ctx.amplicon_name
    mod_type = ctx.mod_type
    mod_freq_df = ctx.modification_frequency_summary_dfs[amp]

    modification_df = mod_freq_df[mod_freq_df['Modification'] == mod_type].copy()
    modification_df.index = [
        '{0} ({1})'.format(batch, batch_index)
        for batch_index, batch in enumerate(modification_df['Batch'], 1)
    ]
    modification_df = modification_df.drop(['Modification', 'Batch'], axis=1)
    modification_df.columns = [
        '{0} ({1})'.format(column, position)
        for position, column in enumerate(modification_df.columns, 1)
    ]

    if ctx.guides_all_same[amp]:
        sgRNA_intervals = [ctx.consensus_sgRNA_intervals[amp]] * modification_df.shape[0]
    else:
        sgRNA_intervals = [ctx.consensus_sgRNA_intervals[amp]]

    amplicon_plot_name = amp + "."
    if len(ctx.amplicon_names) == 1 and amp == "Reference":
        amplicon_plot_name = ""

    plot_name = 'CRISPRessoBatch_percentage_of_{0}_across_alleles_{1}_heatmap'.format(
        mod_type.lower(), amp,
    )
    plot_path = '{0}.html'.format(ctx._jp(plot_name))
    div_id = '{0}-allele-modification-heatmap-{1}'.format(amp.lower(), mod_type.lower())

    return {
        'sample_values': modification_df,
        'sample_sgRNA_intervals': sgRNA_intervals,
        'plot_path': plot_path,
        'title': mod_type,
        'div_id': div_id,
        'amplicon_name': amp,
        'plot_name': plot_name,
    }


def prep_batch_allele_modification_line(ctx):
    """Prepare kwargs for plot_allele_modification_line (Batch/Aggregate).

    Parameters
    ----------
    ctx : BatchPlotContext | AggregatePlotContext
        Must have ``amplicon_name`` and ``mod_type`` scope fields set.

    """
    amp = ctx.amplicon_name
    mod_type = ctx.mod_type
    mod_freq_df = ctx.modification_frequency_summary_dfs[amp]

    modification_df = mod_freq_df[mod_freq_df['Modification'] == mod_type].copy()
    modification_df.index = [
        '{0} ({1})'.format(batch, batch_index)
        for batch_index, batch in enumerate(modification_df['Batch'], 1)
    ]
    modification_df = modification_df.drop(['Modification', 'Batch'], axis=1)
    modification_df.columns = [
        '{0} ({1})'.format(column, position)
        for position, column in enumerate(modification_df.columns, 1)
    ]

    if ctx.guides_all_same[amp]:
        sgRNA_intervals = [ctx.consensus_sgRNA_intervals[amp]] * modification_df.shape[0]
    else:
        sgRNA_intervals = [ctx.consensus_sgRNA_intervals[amp]]

    amplicon_plot_name = amp + "."
    if len(ctx.amplicon_names) == 1 and amp == "Reference":
        amplicon_plot_name = ""

    plot_name = 'CRISPRessoBatch_percentage_of_{0}_across_alleles_{1}_line'.format(
        mod_type.lower(), amp,
    )
    plot_path = '{0}.html'.format(ctx._jp(plot_name))
    div_id = '{0}-allele-modification-line-{1}'.format(amp.lower(), mod_type.lower())

    return {
        'sample_values': modification_df,
        'sample_sgRNA_intervals': sgRNA_intervals,
        'plot_path': plot_path,
        'title': mod_type,
        'div_id': div_id,
        'amplicon_name': amp,
        'plot_name': plot_name,
    }


# =============================================================================
# Multi-mode summary plot prep functions (Pooled, WGS, Aggregate)
# =============================================================================


def _resolve_min_reads_cutoff(ctx):
    """Return the minimum-reads cutoff for the given mode's args namespace.

    Pooled and WGS use ``min_reads_to_use_region``; Aggregate uses
    ``min_reads_for_inclusion``. Returns whichever attribute exists.
    """
    args = ctx.args
    if hasattr(args, 'min_reads_to_use_region'):
        return args.min_reads_to_use_region
    if hasattr(args, 'min_reads_for_inclusion'):
        return args.min_reads_for_inclusion
    raise AttributeError(
        "args namespace lacks both 'min_reads_to_use_region' and "
        "'min_reads_for_inclusion'"
    )


def _sanitize_summary_names(df):
    r"""Return a copy of ``df`` with the ``Name`` column passed through
    ``CRISPRessoShared.slugify``.

    This is a display-label normalization: it strips emoji and other
    non-ASCII characters from user-supplied sample names before they
    are embedded in plot labels.  Without it, plotly's JSON encoder
    produces environment-dependent output (``\uXXXX`` escapes with
    stdlib json vs raw UTF-8 with orjson), which causes integration
    tests to flap between local and CI.

    The source DataFrame on the context is NOT mutated.  CSVs, log
    messages, and ``crispresso2_info`` keys continue to hold the raw
    user-supplied name; only the plot kwargs see the sanitized form.

    Non-string entries in the ``Name`` column are coerced via ``str()``
    before slugify.  If the column is absent, the copy is returned
    unchanged.
    """
    out = df.copy()
    if 'Name' in out.columns:
        out['Name'] = out['Name'].apply(lambda n: CRISPRessoShared.slugify(str(n)))
    return out


def prep_reads_total(ctx, prefix: str, name: str | None = None):
    """Prepare kwargs for plot_reads_total (reads summary bar chart).

    Used by Pooled, WGS, and Aggregate modes.

    Parameters
    ----------
    ctx : PooledPlotContext | WGSPlotContext | AggregatePlotContext
        Must have ``df_summary_quantification``, ``save_png``, ``_jp``,
        and the appropriate min-reads arg for its mode.
    prefix : str
        Filename prefix, e.g. ``'CRISPRessoPooled'``, ``'CRISPRessoWGS'``,
        ``'CRISPRessoAggregate'``.
    name : str, optional
        Override for the bare filename root (no directory). When provided,
        the generated path is ``ctx._jp(name)``; otherwise the path defaults
        to ``ctx._jp(f'{prefix}_reads_summary')``.

    """
    bare = name if name is not None else f'{prefix}_reads_summary'
    return {
        'df_summary_quantification': _sanitize_summary_names(ctx.df_summary_quantification),
        'fig_filename_root': ctx._jp(bare),
        'save_png': ctx.save_png,
        'cutoff': _resolve_min_reads_cutoff(ctx),
    }


def prep_unmod_mod_pcts(ctx, prefix: str, name: str | None = None):
    """Prepare kwargs for plot_unmod_mod_pcts (modification summary bar chart).

    Used by Pooled, WGS, and Aggregate modes.

    Parameters
    ----------
    ctx : PooledPlotContext | WGSPlotContext | AggregatePlotContext
        Must have ``df_summary_quantification``, ``save_png``, ``_jp``,
        and the appropriate min-reads arg for its mode.
    prefix : str
        Filename prefix, e.g. ``'CRISPRessoPooled'``, ``'CRISPRessoWGS'``,
        ``'CRISPRessoAggregate'``.
    name : str, optional
        Override for the bare filename root (no directory). When provided,
        the generated path is ``ctx._jp(name)``; otherwise the path defaults
        to ``ctx._jp(f'{prefix}_modification_summary')``.

    """
    bare = name if name is not None else f'{prefix}_modification_summary'
    return {
        'df_summary_quantification': _sanitize_summary_names(ctx.df_summary_quantification),
        'fig_filename_root': ctx._jp(bare),
        'save_png': ctx.save_png,
        'cutoff': _resolve_min_reads_cutoff(ctx),
    }


# =============================================================================
# Compare prep functions
# =============================================================================


def prep_compare_editing_barchart(ctx):
    """Prepare kwargs for plot_quantification_comparison_barchart.

    Parameters
    ----------
    ctx : ComparePlotContext
        Must have ``amplicon_name``, ``amplicon_info_1/2``,
        ``sample_1_name``, ``sample_2_name`` set.

    Returns
    -------
    dict
        kwargs for ``CRISPRessoPlot.plot_quantification_comparison_barchart``.

    """
    amp = ctx.amplicon_name
    n_refs = len(ctx.amplicon_names)

    n_total_1 = float(ctx.amplicon_info_1[amp]['Reads_aligned'])
    n_unmodified_1 = float(ctx.amplicon_info_1[amp]['Unmodified'])
    n_modified_1 = float(ctx.amplicon_info_1[amp]['Modified'])

    n_total_2 = float(ctx.amplicon_info_2[amp]['Reads_aligned'])
    n_unmodified_2 = float(ctx.amplicon_info_2[amp]['Unmodified'])
    n_modified_2 = float(ctx.amplicon_info_2[amp]['Modified'])

    amplicon_plot_name = amp + "."
    if n_refs == 1 and amp == "Reference":
        amplicon_plot_name = ""

    plot_titles = {
        'vs': plot_title_with_ref_name(
            '%s VS %s' % (ctx.sample_1_name, ctx.sample_2_name), amp, n_refs,
        ),
        'diff': plot_title_with_ref_name(
            '%s - %s' % (ctx.sample_1_name, ctx.sample_2_name), amp, n_refs,
        ),
    }

    return {
        'n_total_1': n_total_1,
        'n_unmodified_1': n_unmodified_1,
        'n_modified_1': n_modified_1,
        'n_total_2': n_total_2,
        'n_unmodified_2': n_unmodified_2,
        'n_modified_2': n_modified_2,
        'sample_1_name': ctx.sample_1_name,
        'sample_2_name': ctx.sample_2_name,
        'plot_titles': plot_titles,
        'plot_path': ctx._jp('1.' + amplicon_plot_name + 'Editing_comparison'),
        'save_also_png': ctx.save_png,
    }


def prep_compare_modification_positions(ctx):
    """Compute Fisher tests and prepare kwargs for plot_quantification_positions.

    Parameters
    ----------
    ctx : ComparePlotContext
        Must have ``amplicon_name``, ``mod_freqs_1/2``, ``mod_type``,
        ``consensus_sequence``, ``cut_points``, ``sgRNA_intervals``,
        ``quant_windows_1/2``, ``sample_1_name/2_name`` set.

    Returns
    -------
    dict with keys:
        ``plot_kwargs`` : dict for ``CRISPRessoPlot.plot_quantification_positions``
        ``mod_df`` : pd.DataFrame for CSV export
        ``sig_count`` : int — significant positions (Bonferroni)
        ``sig_count_quant_window`` : int or np.nan

    """
    from scipy import stats

    amp = ctx.amplicon_name
    mod = ctx.mod_type
    n_refs = len(ctx.amplicon_names)

    mod_counts_1 = np.array(ctx.mod_freqs_1[mod], dtype=float)
    tot_counts_1 = np.array(ctx.mod_freqs_1['Total'], dtype=float)
    unmod_counts_1 = tot_counts_1 - mod_counts_1

    mod_counts_2 = np.array(ctx.mod_freqs_2[mod], dtype=float)
    tot_counts_2 = np.array(ctx.mod_freqs_2['Total'], dtype=float)
    unmod_counts_2 = tot_counts_2 - mod_counts_2

    fisher_results = [
        stats.fisher_exact([[z[0], z[1]], [z[2], z[3]]]) if max(z) > 0 else [np.nan, 1.0]
        for z in zip(mod_counts_1, unmod_counts_1, mod_counts_2, unmod_counts_2)
    ]
    oddsratios, pvalues = [a for a, b in fisher_results], [b for a, b in fisher_results]

    # Build mod_df for CSV export
    rows = [
        [ctx.sample_1_name + '_' + mod] + list(mod_counts_1),
        [ctx.sample_1_name + '_total'] + list(tot_counts_1),
        [ctx.sample_2_name + '_' + mod] + list(mod_counts_2),
        [ctx.sample_2_name + '_total'] + list(tot_counts_2),
        ['odds_ratios'] + list(oddsratios),
        ['pvalues'] + list(pvalues),
    ]

    m, pvals = len(pvalues), np.asarray(pvalues)
    qval_bonferroni = pvals * float(m)
    qval_bonferroni[np.where(qval_bonferroni > 1)] = 1
    rows.append(['qval_bonferroni'] + list(qval_bonferroni))

    colnames = ['Reference'] + list(ctx.consensus_sequence)
    mod_df = pd.DataFrame(rows, columns=colnames)

    # Significance counts
    sig_count = len(np.where(qval_bonferroni <= ctx.args.reported_qvalue_cutoff)[0])

    sig_count_quant_window = np.nan
    quant_windows_are_equal = np.array_equal(ctx.quant_windows_1, ctx.quant_windows_2)
    if quant_windows_are_equal:
        qvals_in_quant = np.take(qval_bonferroni, ctx.quant_windows_1)
        sig_count_quant_window = len(np.where(qvals_in_quant <= ctx.args.reported_qvalue_cutoff)[0])

    amplicon_plot_name = amp + "."
    if n_refs == 1 and amp == "Reference":
        amplicon_plot_name = ""

    mod_name = mod
    if mod == "All_modifications":
        mod_name = "Combined modifications (insertions, deletions and substitutions)"

    plot_title = plot_title_with_ref_name(
        '%s: %s - %s' % (mod, ctx.sample_1_name, ctx.sample_2_name), amp, n_refs,
    )

    plot_kwargs = {
        'mod_counts_1': mod_counts_1,
        'tot_counts_1': tot_counts_1,
        'mod_counts_2': mod_counts_2,
        'tot_counts_2': tot_counts_2,
        'len_amplicon': len(ctx.consensus_sequence),
        'pvalues': pvalues,
        'consensus_sequence_len': len(ctx.consensus_sequence),
        'cut_points': ctx.cut_points,
        'sgRNA_intervals': ctx.sgRNA_intervals,
        'plot_title': plot_title,
        'plot_path': ctx._jp('2.' + amplicon_plot_name + mod + '_quantification'),
        'save_also_png': ctx.save_png,
    }

    return {
        'plot_kwargs': plot_kwargs,
        'mod_df': mod_df,
        'mod_name': mod_name,
        'sig_count': sig_count,
        'sig_count_quant_window': sig_count_quant_window,
        'amplicon_plot_name': amplicon_plot_name,
    }


def prep_compare_allele_table(ctx, pair_index):
    """Merge allele tables and prepare kwargs for plot_alleles_table_compare.

    Parameters
    ----------
    ctx : ComparePlotContext
        Must have ``allele_pairs``, ``sample_1_name/2_name``,
        ``consensus_sequence``, ``sgRNA_intervals``, ``cut_points``,
        and ``args`` set.
    pair_index : int
        Index into ``ctx.allele_pairs``.

    Returns
    -------
    dict with keys:
        ``merged_df`` : pd.DataFrame for CSV export
        ``ref_seq_around_cut`` : str
        ``plot_top_kwargs`` : dict for plot_alleles_table_compare (LFC desc)
        ``plot_bottom_kwargs`` : dict for plot_alleles_table_compare (LFC asc)
        ``is_base_edit`` : bool
        ``file_root`` : str — allele file root for filenames

    """
    allele_file_1, allele_file_2, df1, df2 = ctx.allele_pairs[pair_index]

    # Find unmodified reference for comparison
    ref_seq_around_cut = ""
    if len(df1.loc[df1['Reference_Sequence'].str.contains('-') == False]) > 0:
        ref_seq_around_cut = df1.loc[df1['Reference_Sequence'].str.contains('-') == False]['Reference_Sequence'].iloc[0]
    elif len(df2.loc[df2['Reference_Sequence'].str.contains('-') == False]) > 0:
        ref_seq_around_cut = df2.loc[df2['Reference_Sequence'].str.contains('-') == False]['Reference_Sequence'].iloc[0]
    else:
        seq_len = df2[df2['Unedited'] == True]['Reference_Sequence'].iloc[0]
        for sgRNA_interval, cut_point in zip(ctx.sgRNA_intervals, ctx.cut_points):
            sgRNA_seq = ctx.consensus_sequence[sgRNA_interval[0]:sgRNA_interval[1]]
            if sgRNA_seq in allele_file_1:
                ref_seq_around_cut = ctx.consensus_sequence[
                    max(0, cut_point - ctx.args.offset_around_cut_to_plot + 1):
                    min(seq_len, cut_point + ctx.args.offset_around_cut_to_plot + 1)
                ]
                break

    # Merge
    merged = pd.merge(
        df1, df2,
        on=['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated'],
        suffixes=('_' + ctx.sample_1_name, '_' + ctx.sample_2_name),
        how='outer',
    )
    quant_cols = [
        '#Reads_' + ctx.sample_1_name, '%Reads_' + ctx.sample_1_name,
        '#Reads_' + ctx.sample_2_name, '%Reads_' + ctx.sample_2_name,
    ]
    merged[quant_cols] = merged[quant_cols].fillna(0)
    lfc_error = 0.1
    merged['each_LFC'] = np.log2(
        ((merged['%Reads_' + ctx.sample_1_name] + lfc_error) /
         (merged['%Reads_' + ctx.sample_2_name] + lfc_error)).astype(float)
    ).replace([np.inf, np.nan], 0)
    merged = merged.sort_values(
        ['%Reads_' + ctx.sample_1_name, 'Aligned_Sequence', 'Reference_Sequence'],
        ascending=[False, True, True],
    )
    merged = merged.reset_index(drop=True).set_index('Aligned_Sequence')

    file_root = os.path.split(allele_file_1)[1].replace(".txt", "")
    is_base_edit = 'base_edit' in allele_file_1

    common_kwargs = {
        'reference_seq': ref_seq_around_cut,
        'sample_name_1': ctx.sample_1_name,
        'sample_name_2': ctx.sample_2_name,
        'MIN_FREQUENCY': ctx.args.min_frequency_alleles_around_cut_to_plot,
        'MAX_N_ROWS': ctx.args.max_rows_alleles_around_cut_to_plot,
        'SAVE_ALSO_PNG': ctx.save_png,
    }

    # Tiebreaker sort: Aligned_Sequence/Reference_Sequence ascending.
    # Keep in sync with the original CRISPRessoCompareCORE implementation
    # so the alleles selected (MAX_N_ROWS) and their order match the
    # expected PDF text streams in the integration tests.
    plot_top_kwargs = dict(common_kwargs)
    plot_top_kwargs['df_alleles'] = merged.sort_values(
        ['each_LFC', 'Aligned_Sequence', 'Reference_Sequence'],
        ascending=[False, True, True],
    )
    plot_top_kwargs['fig_filename_root'] = ctx._jp('3.' + file_root + '_top')

    plot_bottom_kwargs = dict(common_kwargs)
    plot_bottom_kwargs['df_alleles'] = merged.sort_values(
        ['each_LFC', 'Aligned_Sequence', 'Reference_Sequence'],
        ascending=[True, True, True],
    )
    plot_bottom_kwargs['fig_filename_root'] = ctx._jp('3.' + file_root + '_bottom')

    return {
        'merged_df': merged,
        'ref_seq_around_cut': ref_seq_around_cut,
        'plot_top_kwargs': plot_top_kwargs,
        'plot_bottom_kwargs': plot_bottom_kwargs,
        'is_base_edit': is_base_edit,
        'file_root': file_root,
    }


# =============================================================================
# Shared data-file exports (used by both CORE and Pro hooks)
# =============================================================================


def write_core_amplicon_data_files(ctx: CorePlotContext, crispresso2_info: dict):
    """Write CORE's per-amplicon CSV data files.

    Replicates the CSV exports that ``CRISPRessoCORE.main``'s ``elif not
    args.suppress_plots:`` branch performs at the top of its per-reference
    plotting loop.  Both CORE (non-Pro path) and Pro's ``on_plots_complete``
    hook call this so that the same data files are produced regardless of
    which plotting backend handles the figures.

    Writes, for each reference with ``counts_total[ref_name] > 0``:

    - ``{ref_plot_name}Nucleotide_frequency_table.txt``
    - ``{ref_plot_name}Nucleotide_percentage_table.txt``
    - ``{ref_plot_name}Quantification_window_nucleotide_frequency_table.txt``
    - ``{ref_plot_name}Quantification_window_nucleotide_percentage_table.txt``

    And, when ``args.base_editor_output`` is set:

    - ``{ref_plot_name}Substitution_frequency_table.txt``
    - ``{ref_plot_name}Quantification_window_substitution_frequency_table.txt``

    Each file's basename is recorded under the standard key in
    ``crispresso2_info['results']['refs'][ref_name]``.

    Skipped entirely when ``args.crispresso1_mode`` is True (mirrors CORE).

    Parameters
    ----------
    ctx : CorePlotContext
        Context exposing ``ref_names``, ``refs``, ``counts_total``,
        ``all_base_count_vectors``, ``all_substitution_base_vectors``,
        ``substitution_base_vectors``, and ``_jp``.
    crispresso2_info : dict
        Run info dict whose ``results.refs[ref_name]`` entries are updated
        with the written filenames (matching CORE's key names).

    """
    if getattr(ctx.args, 'crispresso1_mode', False):
        return

    for ref_name in ctx.ref_names:
        if ctx.counts_total[ref_name] < 1:
            continue

        ref = ctx.refs[ref_name]
        # Gracefully skip refs that lack the fields we need.  This happens
        # in minimal test fixtures; production refs always have these.
        if not all(k in ref for k in ('ref_plot_name', 'sequence', 'include_idxs')):
            continue
        # Skip if the per-amplicon count vectors aren't available.
        if (ref_name + '_A') not in ctx.all_base_count_vectors:
            continue

        ref_plot_name = ref['ref_plot_name']
        ref_seq = ref['sequence']
        include_idxs_list = ref['include_idxs']
        quantification_window_ref_seq = [ref_seq[x] for x in include_idxs_list]
        tot_aln_reads = ctx.counts_total[ref_name]
        info_ref = crispresso2_info['results']['refs'][ref_name]

        # Derive quantification-window counts from the full-amplicon vectors.
        window_rows = [
            [
                ctx.all_base_count_vectors[ref_name + '_' + nuc][x]
                for x in include_idxs_list
            ]
            for nuc in ['A', 'C', 'G', 'T', 'N', '-']
        ]
        df_nuc_freq = pd.DataFrame(
            window_rows,
            index=['A', 'C', 'G', 'T', 'N', '-'],
            columns=quantification_window_ref_seq,
        )
        quant_window_nuc_freq_filename = ctx._jp(
            ref_plot_name + 'Quantification_window_nucleotide_frequency_table.txt'
        )
        df_nuc_freq.to_csv(quant_window_nuc_freq_filename, sep='\t', header=True, index=True)
        info_ref['quant_window_nuc_freq_filename'] = os.path.basename(quant_window_nuc_freq_filename)

        df_nuc_pct = df_nuc_freq.divide(tot_aln_reads)
        quant_window_nuc_pct_filename = ctx._jp(
            ref_plot_name + 'Quantification_window_nucleotide_percentage_table.txt'
        )
        df_nuc_pct.to_csv(quant_window_nuc_pct_filename, sep='\t', header=True, index=True)
        info_ref['quant_window_nuc_pct_filename'] = os.path.basename(quant_window_nuc_pct_filename)

        df_nuc_freq_all = pd.DataFrame(
            [
                ctx.all_base_count_vectors[ref_name + '_' + nuc]
                for nuc in ['A', 'C', 'G', 'T', 'N', '-']
            ],
            index=['A', 'C', 'G', 'T', 'N', '-'],
            columns=list(ref_seq),
        )
        nuc_freq_filename = ctx._jp(ref_plot_name + 'Nucleotide_frequency_table.txt')
        df_nuc_freq_all.to_csv(nuc_freq_filename, sep='\t', header=True, index=True)
        info_ref['nuc_freq_filename'] = os.path.basename(nuc_freq_filename)

        df_nuc_pct_all = df_nuc_freq_all.divide(tot_aln_reads)
        nuc_pct_filename = ctx._jp(ref_plot_name + 'Nucleotide_percentage_table.txt')
        df_nuc_pct_all.to_csv(nuc_pct_filename, sep='\t', header=True, index=True)
        info_ref['nuc_pct_filename'] = os.path.basename(nuc_pct_filename)

        if getattr(ctx.args, 'base_editor_output', False):
            alph = ['A', 'C', 'G', 'T', 'N']
            df_sub_freq = pd.DataFrame(
                [list(ctx.substitution_base_vectors[ref_name + '_' + a]) for a in alph],
                index=alph,
                columns=quantification_window_ref_seq,
            )
            quant_window_sub_freq_filename = ctx._jp(
                ref_plot_name + 'Quantification_window_substitution_frequency_table.txt'
            )
            df_sub_freq.to_csv(quant_window_sub_freq_filename, sep='\t', header=True, index=True)
            info_ref['quant_window_sub_freq_filename'] = os.path.basename(quant_window_sub_freq_filename)

            df_sub_freq_all = pd.DataFrame(
                [list(ctx.all_substitution_base_vectors[ref_name + '_' + a]) for a in alph],
                index=alph,
                columns=list(ref_seq),
            )
            sub_freq_table_filename = ctx._jp(
                ref_plot_name + 'Substitution_frequency_table.txt'
            )
            df_sub_freq_all.to_csv(sub_freq_table_filename, sep='\t', header=True, index=True)
            info_ref['sub_freq_table_filename'] = os.path.basename(sub_freq_table_filename)


def write_core_per_sgRNA_data_files(ctx: CorePlotContext, crispresso2_info: dict):
    """Write CORE's per-sgRNA CSV data files.

    Mirrors the per-sgRNA CSV exports from ``CRISPRessoCORE.main``'s ``elif
    not args.suppress_plots:`` branch.  Called by Pro's ``on_plots_complete``
    so Pro mode produces the same on-disk data files as the non-Pro path.

    For each reference amplicon with reads, iterates ``sgRNA_sequences`` and
    writes:

    - ``{ref_plot_name}Alleles_frequency_table_around_{sgRNA_label}.txt`` (plot 9)

    When base-editor output is requested and not crispresso1 mode:

    - ``{ref_plot_name}Selected_nucleotide_percentage_table_around_{sgRNA_label}.txt`` (plot 10e/f/g)
    - ``{ref_plot_name}Selected_nucleotide_frequency_table_around_{sgRNA_label}.txt`` (plot 10e/f/g)
    - ``10i.{ref_name}.{sgRNA_label}.{binary_allele_counts,category_allele_counts,precise_allele_counts,arrays,counts}.txt`` (plot 10i)

    Writes per-amplicon base-edit alleles-around-cut table (plot 10h) when
    base-editor output is requested.

    Requires ``ctx.df_alleles`` to be populated.  Mutates
    ``ctx.ref_name`` and ``ctx.sgRNA_ind`` during iteration but restores
    them to ``None`` on return.
    """
    saved_ref_name = ctx.ref_name
    saved_sgRNA_ind = ctx.sgRNA_ind
    is_base_editor = getattr(ctx.args, 'base_editor_output', False)
    is_crispresso1 = getattr(ctx.args, 'crispresso1_mode', False)

    try:
        for ref_name in ctx.ref_names:
            if ctx.counts_total[ref_name] < 1:
                continue

            ref = ctx.refs[ref_name]
            # Gracefully skip refs that lack the fields we need.
            if 'ref_plot_name' not in ref or 'sgRNA_sequences' not in ref:
                continue
            if ctx.df_alleles is None:
                continue

            ref_plot_name = ref['ref_plot_name']
            info_ref = crispresso2_info['results']['refs'][ref_name]
            info_ref.setdefault('allele_frequency_files', [])
            ctx.ref_name = ref_name

            # Plot 10h (base-edit quilt, per amplicon, uses last sgRNA)
            if is_base_editor and not is_crispresso1 and len(ref.get('sgRNA_cut_points', [])) > 0:
                ctx.sgRNA_ind = len(ref['sgRNA_cut_points']) - 1
                base_edit_data = prep_base_edit_quilt(ctx)
                base_edit_allele_filename = ctx._jp(
                    ref_plot_name + 'base_edit_'
                    + ctx.args.conversion_nuc_from + 's_quilt.txt'
                )
                base_edit_data['df_alleles_around_cut'].to_csv(
                    base_edit_allele_filename,
                    sep='\t',
                    header=True,
                    index=True,
                    index_label=base_edit_data['df_alleles_around_cut'].index.name
                )
                info_ref['allele_frequency_files'].append(
                    os.path.basename(base_edit_allele_filename)
                )

            # Per sgRNA writes
            sgRNA_sequences = ref['sgRNA_sequences']
            for sgRNA_ind, _sgRNA_seq in enumerate(sgRNA_sequences):
                sgRNA = ref['sgRNA_orig_sequences'][sgRNA_ind]
                sgRNA_name = ref['sgRNA_names'][sgRNA_ind]
                sgRNA_label = "sgRNA_" + sgRNA
                if sgRNA_name:
                    sgRNA_label = sgRNA_name
                sgRNA_label = CRISPRessoShared.slugify(sgRNA_label)
                ctx.sgRNA_ind = sgRNA_ind

                # Plot 9: alleles around cut
                alleles_data = prep_alleles_around_cut(ctx)
                allele_filename = ctx._jp(
                    ref_plot_name + 'Alleles_frequency_table_around_'
                    + sgRNA_label + '.txt'
                )
                alleles_data['df_alleles_around_cut'].to_csv(
                    allele_filename, sep='\t', header=True, index=True, index_label=alleles_data['df_alleles_around_cut'].index.name
                )
                info_ref['allele_frequency_files'].append(
                    os.path.basename(allele_filename)
                )

                # Plot 10e/f/g: selected nucleotide frequency/percentage
                if is_base_editor and not is_crispresso1:
                    sel_nuc_data = prep_conversion_at_sel_nucs(ctx)
                    just_sel_nuc_pcts = sel_nuc_data['just_sel_nuc_pcts']
                    just_sel_nuc_freqs = sel_nuc_data['just_sel_nuc_freqs']
                    quant_window_sel_nuc_pct_filename = ctx._jp(
                        ref_plot_name + 'Selected_nucleotide_percentage_table_around_'
                        + sgRNA_label + '.txt'
                    )
                    just_sel_nuc_pcts.to_csv(
                        quant_window_sel_nuc_pct_filename, sep='\t', header=True, index=True,
                    )
                    quant_window_sel_nuc_freq_filename = ctx._jp(
                        ref_plot_name + 'Selected_nucleotide_frequency_table_around_'
                        + sgRNA_label + '.txt'
                    )
                    just_sel_nuc_freqs.to_csv(
                        quant_window_sel_nuc_freq_filename, sep='\t', header=True, index=True,
                    )

                    # Plot 10i: base-edit upset counts
                    upset_data = prep_base_edit_upset(
                        ref['sequence'],
                        ctx.df_alleles,
                        ref_name,
                        ref['sgRNA_intervals'][sgRNA_ind],
                        ctx.args,
                        gap_incentive=ref['gap_incentive'],
                        sgRNA_legend=sgRNA_label,
                    )
                    if upset_data is not None:
                        write_base_edit_counts(
                            ref_name + '.' + sgRNA_label,
                            upset_data['counts_dict'],
                            upset_data['bp_substitutions_arr'],
                            ctx._jp,
                        )
    finally:
        ctx.ref_name = saved_ref_name
        ctx.sgRNA_ind = saved_sgRNA_ind


def write_core_coding_data_files(ctx: CorePlotContext, crispresso2_info: dict):
    """Write CORE's per-coding-sequence CSV data files (amino acid tables).

    Mirrors the per-coding-sequence CSV exports from CORE's ``elif`` branch.
    Called by Pro's ``on_plots_complete`` so Pro mode produces the same
    on-disk data files as the non-Pro path.

    For each reference amplicon that ``contains_coding_seq``, iterates
    ``running_info['coding_seqs']`` and writes
    ``{ref_plot_name}amino_acid_table_for_{coding_seq_label}.txt``.

    Mutates ``ctx.ref_name`` and ``ctx.coding_seq_ind`` during iteration
    but restores them on return.
    """
    saved_ref_name = ctx.ref_name
    saved_coding_seq_ind = ctx.coding_seq_ind
    coding_seqs = crispresso2_info['running_info'].get('coding_seqs', [])
    coding_seq_names = crispresso2_info['running_info'].get('coding_seq_names', [])
    if not coding_seqs:
        return

    try:
        for ref_name in ctx.ref_names:
            if ctx.counts_total[ref_name] < 1:
                continue
            ref = ctx.refs[ref_name]
            if not ref.get('contains_coding_seq', False):
                continue
            if 'ref_plot_name' not in ref:
                continue

            ref_plot_name = ref['ref_plot_name']
            ctx.ref_name = ref_name

            for i, _coding_seq in enumerate(coding_seqs):
                coding_seq_label = (
                    coding_seq_names[i]
                    if i < len(coding_seq_names)
                    else str(i)
                )
                ctx.coding_seq_ind = i
                amino_acid_data = prep_amino_acid_table(ctx)
                df_to_plot = amino_acid_data['df_to_plot']
                amino_acid_filename = ctx._jp(
                    ref_plot_name + 'amino_acid_table_for_'
                    + coding_seq_label + '.txt'
                )
                df_to_plot.to_csv(
                    amino_acid_filename, sep='\t', header=True, index=True,
                )
    finally:
        ctx.ref_name = saved_ref_name
        ctx.coding_seq_ind = saved_coding_seq_ind


def write_core_hdr_data_files(ctx: CorePlotContext, crispresso2_info: dict):
    """Write HDR cross-amplicon nucleotide/modification summary tables.

    Mirrors plot_4g's CSV exports from CORE's elif branch.  Only writes when
    ``args.expected_hdr_amplicon_seq`` is non-empty.  Writes:

    - ``{ref_plot_name}Reads_from_all_amplicons_nucleotide_percent_table.txt``
    - ``{ref_plot_name}Reads_from_all_amplicons_modification_percent_table.txt``
    """
    if not getattr(ctx.args, 'expected_hdr_amplicon_seq', ''):
        return

    ref_names_for_hdr = [r for r in ctx.ref_names if ctx.counts_total[r] > 0]
    if not ref_names_for_hdr:
        return
    if 'ref_plot_name' not in ctx.refs[ref_names_for_hdr[0]]:
        return

    saved_ref_name = ctx.ref_name
    ctx.ref_name = ctx.ref_names[0]
    try:
        plot_4g_input = prep_hdr_nucleotide_quilt(ctx)
    finally:
        ctx.ref_name = saved_ref_name

    hdr_nuc_pct_df = plot_4g_input['nuc_pct_df']
    hdr_mod_pct_df = plot_4g_input['mod_pct_df']
    ref_plot_name = ctx.refs[ref_names_for_hdr[0]]['ref_plot_name']

    mod_freq_filename = ctx._jp(
        ref_plot_name + 'Reads_from_all_amplicons_modification_percent_table.txt'
    )
    hdr_mod_pct_df.rename(columns={'Batch': 'Amplicon'}).to_csv(
        mod_freq_filename, sep='\t', header=True, index=False,
    )
    nuc_freq_filename = ctx._jp(
        ref_plot_name + 'Reads_from_all_amplicons_nucleotide_percent_table.txt'
    )
    hdr_nuc_pct_df.rename(columns={'Batch': 'Amplicon'}).to_csv(
        nuc_freq_filename, sep='\t', header=True, index=False,
    )


def get_pe_scaffold_search(prime_edited_ref_sequence, prime_editing_pegRNA_extension_seq, prime_editing_pegRNA_scaffold_seq, prime_editing_pegRNA_scaffold_min_match_length):
    """For prime editing, determines the scaffold string to search for (the shortest substring of args.prime_editing_pegRNA_scaffold_seq not in the prime-edited reference sequence)
    params:
     prime_edited_ref_sequence: reference sequence of the prime-edited sequence
     prime_editing_extension_seq: RNA sequence of extension sequence
     prime_editing_pegRNA_scaffold_seq: sequence of the scaffold sequence
     prime_editing_pegRNA_scaffold_min_match_length: minimum number of bases required to match between scaffold and read to count as scaffold-incorporated
    returns:
     tuple of(
         index of location in ref to find scaffold seq if it exists
         shortest dna sequence to identify scaffold sequence
     )
    """
    logger.info('Processing pegRNA scaffold sequence...')
    # first, define the sequence we are looking for (extension plus the first base(s) of the scaffold)
    scaffold_dna = CRISPRessoShared.reverse_complement(prime_editing_pegRNA_scaffold_seq.upper().replace('U', 'T'))

    extension_seq_dna_top_strand = prime_editing_pegRNA_extension_seq.upper().replace('U', 'T')
    prime_editing_extension_seq_dna = CRISPRessoShared.reverse_complement(extension_seq_dna_top_strand)

    scaffold_start_loc = prime_edited_ref_sequence.index(prime_editing_extension_seq_dna) + len(prime_editing_extension_seq_dna)

    # next find min length scaffold that when combined with min extension sequence is not in edited sequence
    len_scaffold_to_use = prime_editing_pegRNA_scaffold_min_match_length  # length of the scaffold sequence to be added to the extension sequence, start at 1bp
    scaffold_dna_search = prime_editing_extension_seq_dna + scaffold_dna[0:len_scaffold_to_use]
    while scaffold_dna_search in prime_edited_ref_sequence:
        if len_scaffold_to_use > len(scaffold_dna):
            raise CRISPRessoShared.BadParameterException('The DNA scaffold provided is found in the unedited reference sequence. Please provide a longer scaffold sequence.')
        len_scaffold_to_use += 1
        scaffold_dna_search = prime_editing_extension_seq_dna + scaffold_dna[0:len_scaffold_to_use]

    logger.info('Searching for scaffold-templated reads with the sequence: \'' + str(scaffold_dna[0:len_scaffold_to_use]) + '\' starting at position ' + str(scaffold_start_loc) + ' in reads that align to the prime-edited sequence')
    return (scaffold_start_loc, scaffold_dna[0:len_scaffold_to_use])


def write_core_prime_editing_data_files(ctx: CorePlotContext, crispresso2_info: dict):
    """Write prime editing scaffold insertion sizes table.

    Mirrors the plot_11c CSV export from CORE's elif branch.  Only writes
    when both ``args.prime_editing_pegRNA_extension_seq`` and
    ``args.prime_editing_pegRNA_scaffold_seq`` are set.

    Computes the DataFrame on the fly by iterating
    ``ctx.df_alleles`` restricted to ``Reference_Name == 'Scaffold-incorporated'``
    and measuring the match length of the observed scaffold against the
    expected (reverse-complemented) sequence.  Also stores the computed
    DataFrame on ``ctx.df_scaffold_insertion_sizes`` so subsequent plot
    preparation can reuse it.
    """
    args = ctx.args
    if not getattr(args, 'prime_editing_pegRNA_extension_seq', ''):
        return
    if not getattr(args, 'prime_editing_pegRNA_scaffold_seq', ''):
        return

    refs = ctx.refs
    if 'Prime-edited' not in refs:
        return

    scaffold_dna_seq = CRISPRessoShared.reverse_complement(
        args.prime_editing_pegRNA_scaffold_seq.upper().replace('U', 'T')
    )
    pe_seq = refs['Prime-edited']['sequence']
    pe_scaffold_dna_info = get_pe_scaffold_search(
        pe_seq,
        args.prime_editing_pegRNA_extension_seq,
        args.prime_editing_pegRNA_scaffold_seq,
        args.prime_editing_pegRNA_scaffold_min_match_length,
    )

    df_alleles_scaffold = ctx.df_alleles.loc[
        ctx.df_alleles['Reference_Name'] == 'Scaffold-incorporated'
    ]

    def get_scaffold_len(row, scaffold_start_loc, scaffold_seq):
        pe_read_possible_scaffold_loc = row['ref_positions'].index(scaffold_start_loc - 1) + 1
        i = 0
        num_match_scaffold = 0
        num_gaps = 0
        matches_scaffold_to_this_point = True
        has_gaps_to_this_point = True
        aln_seq_to_test = row['Aligned_Sequence'][pe_read_possible_scaffold_loc:].replace("-", "")
        while (
            i < len(scaffold_seq)
            and i < len(aln_seq_to_test)
            and (matches_scaffold_to_this_point or has_gaps_to_this_point)
        ):
            if matches_scaffold_to_this_point and aln_seq_to_test[i] == scaffold_seq[i]:
                num_match_scaffold += 1
            else:
                matches_scaffold_to_this_point = False

            if (
                has_gaps_to_this_point
                and row['Reference_Sequence'][pe_read_possible_scaffold_loc + i] == '-'
            ):
                num_gaps += 1
            else:
                has_gaps_to_this_point = False
            i += 1
        return (
            row['Aligned_Sequence'],
            row['Reference_Sequence'],
            num_match_scaffold,
            num_gaps,
            row['#Reads'],
            row['%Reads'],
        )

    df_scaffold_insertion_sizes = pd.DataFrame(
        list(
            df_alleles_scaffold.apply(
                lambda row: get_scaffold_len(row, pe_scaffold_dna_info[0], scaffold_dna_seq),
                axis=1,
            ).values
        ),
        columns=[
            'Aligned_Sequence', 'Reference_Sequence',
            'Num_match_scaffold', 'Num_gaps', '#Reads', '%Reads',
        ],
    )

    filename = ctx._jp('Scaffold_insertion_sizes.txt')
    df_scaffold_insertion_sizes.to_csv(filename, sep='\t', header=True, index=False)
    crispresso2_info['Scaffold_insertion_sizes_filename'] = filename
    # Stash so downstream plot code can reuse.
    ctx.df_scaffold_insertion_sizes = df_scaffold_insertion_sizes


def write_all_core_data_files(ctx: CorePlotContext, crispresso2_info: dict):
    """Write all CORE CSV data files that the non-Pro elif branch writes.

    Convenience wrapper that runs :func:`write_core_amplicon_data_files`,
    :func:`write_core_hdr_data_files`, :func:`write_core_per_sgRNA_data_files`,
    :func:`write_core_coding_data_files`, and
    :func:`write_core_prime_editing_data_files` in the correct order
    (amplicon writes first so HDR/per-sgRNA writes can reference the
    metadata they populate).
    """
    write_core_amplicon_data_files(ctx, crispresso2_info)
    write_core_hdr_data_files(ctx, crispresso2_info)
    write_core_per_sgRNA_data_files(ctx, crispresso2_info)
    write_core_coding_data_files(ctx, crispresso2_info)
    write_core_prime_editing_data_files(ctx, crispresso2_info)
