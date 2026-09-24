# -*- coding: utf-8 -*-
"""CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
"""

import os
from concurrent.futures import ProcessPoolExecutor, wait
from copy import deepcopy
from functools import partial
import sys
import re
import traceback
from datetime import datetime
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2.CRISPRessoReports import CRISPRessoReport
from CRISPResso2.plots.data_prep import (
    prep_batch_conversion_map,
    prep_batch_conversion_map_around_sgRNA,
    prep_batch_nuc_quilt,
    prep_batch_nuc_quilt_around_sgRNA,
)

C2PRO_INSTALLED = CRISPRessoShared.is_C2Pro_installed()

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error = logger.critical
warn = logger.warning
debug = logger.debug
info = logger.info

_ROOT = os.path.abspath(os.path.dirname(__file__))


# Support functions###
def check_library(library_name):
    try:
        return __import__(library_name)
    except:
        error('You need to install %s module to use CRISPRessoBatch!' % library_name)
        sys.exit(1)


pd = check_library('pandas')
np = check_library('numpy')


def should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff=300):
    """Determine if large plots should be plotted.

    Parameters
    ----------
    num_rows : int
        Number of rows in the dataframe.
    c2pro_installed : bool
        Whether CRISPRessoPro is installed.
    use_matplotlib : bool
        Whether to use matplotlib when CRISPRessoPro is installed, i.e. value
        of `--use_matplotlib`.
    large_plot_cutoff : int, optional
        Number of samples at which to not plot large plots with matplotlib.
        Note that each sample has 6 rows in the datafame. Defaults to 300.

    Returns
    -------
    bool
        Whether to plot large plots.

    """
    return (
        (not use_matplotlib and c2pro_installed)
        or (num_rows / 6) < large_plot_cutoff
    )


def main():
    try:
        start_time = datetime.now()
        start_time_string = start_time.strftime('%Y-%m-%d %H:%M:%S')

        # if no args are given, print a simplified help message
        if len(sys.argv) == 1:
            raise CRISPRessoShared.BadParameterException(CRISPRessoShared.format_cl_text('usage: CRISPRessoBatch  [-bs BATCH_SETTINGS]  [-n NAME]\n' +
                'commonly-used arguments:\n' +
                '-h, --help            show the full list of arguments\n' +
                '-v, --version         show program\'s version number and exit\n' +
                '-bs BATCH_SETTINGS    Tab-separated file where rows are samples and columns specify settings for each sample.\n' +
                '-n NAME, --name NAME  Name for the analysis (default: name based on input file name)'
            ))

        parser = CRISPRessoShared.getCRISPRessoArgParser("Batch", parser_title='CRISPRessoBatch Parameters')

        args = parser.parse_args()

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

        description = ['~~~CRISPRessoBatch~~~', '-Analysis of CRISPR/Cas9 outcomes from batch deep sequencing data-']
        batch_string = r'''
 _________________
| __    ___ __    |
||__) /\ | /  |__||
||__)/--\| \__|  ||
|_________________|
        '''
        info(CRISPRessoShared.get_crispresso_header(description, batch_string))

        crispresso_options = CRISPRessoShared.get_core_crispresso_options()
        options_to_ignore = {'name', 'output_folder', 'zip_output'}
        crispresso_options_for_batch = list(crispresso_options - options_to_ignore)

        CRISPRessoShared.check_file(args.batch_settings)
        custom_config = CRISPRessoShared.check_custom_config(args)

        if args.zip_output and not args.place_report_in_output_folder:
            warn('Invalid arguement combination: If zip_output is True then place_report_in_output_folder must also be True. Setting place_report_in_output_folder to True.')
            args.place_report_in_output_folder = True

        batch_folder_name = os.path.splitext(os.path.basename(args.batch_settings))[0]
        if args.name and args.name != "":
            clean_name = CRISPRessoShared.slugify(args.name)
            if args.name != clean_name:
                warn(
                    'The specified name {0} contained invalid characters and was changed to: {1}'.format(
                        args.name, clean_name,
                    ),
                )
            batch_folder_name = clean_name

        output_folder_name = 'CRISPRessoBatch_on_%s' % batch_folder_name
        OUTPUT_DIRECTORY = os.path.abspath(output_folder_name)

        if args.batch_output_folder:
                 OUTPUT_DIRECTORY = os.path.join(os.path.abspath(args.batch_output_folder), output_folder_name)

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename)  # handy function to put a file in the output directory

        from CRISPResso2.plots import CRISPRessoPlot
        CRISPRessoPlot.setMatplotlibDefaults()

        try:
            info('Creating Folder %s' % OUTPUT_DIRECTORY, {'percent_complete': 0})
            os.makedirs(OUTPUT_DIRECTORY)
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename = _jp('CRISPRessoBatch_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        status_handler = CRISPRessoShared.StatusHandler(os.path.join(OUTPUT_DIRECTORY, 'CRISPRessoBatch_status.json'))
        logger.addHandler(status_handler)

        with open(log_filename, 'w+') as outfile:
            outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Batch_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2Batch_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}}  # keep track of all information for this run to be saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        n_processes_for_batch = 1
        if args.n_processes == "max":
            n_processes_for_batch = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes_for_batch = int(args.n_processes)

        crispresso_cmd_to_write = ' '.join(sys.argv)
        if args.write_cleaned_report:
            cmd_copy = sys.argv[:]
            cmd_copy[0] = 'CRISPRessoBatch'
            for i in range(len(cmd_copy)):
                if os.sep in cmd_copy[i]:
                    cmd_copy[i] = os.path.basename(cmd_copy[i])

            crispresso_cmd_to_write = ' '.join(cmd_copy)  # clean command doesn't show the absolute path to the executable or other files
        crispresso2_info['running_info']['command_used'] = crispresso_cmd_to_write

        # parse excel sheet
        batch_params = pd.read_csv(args.batch_settings, comment='#', sep='\t')
        # pandas either allows for auto-detect sep or for comment. not both
#        batch_params=pd.read_csv(args.batch_settings,sep=None,engine='python',error_bad_lines=False)
        batch_params.columns = batch_params.columns.str.strip(' -\xd0')

        # if there are more processes than batches, use more processes on each sub-crispresso run
        # if there are more batches than processes, just use 1 process on each sub-crispresso run
        args.n_processes = 1
        num_batches = batch_params.shape[0]
        if int(n_processes_for_batch / num_batches) > 1:
            args.n_processes = int(n_processes_for_batch / num_batches)

        int_columns = ['default_min_aln_score', 'min_average_read_quality', 'min_single_bp_quality',
                       'min_bp_quality_or_N',
                       'quantification_window_size', 'quantification_window_center', 'exclude_bp_from_left',
                       'exclude_bp_from_right',
                       'plot_window_size', 'max_rows_alleles_around_cut_to_plot']
        for int_col in int_columns:
            if int_col in batch_params.columns:
                batch_params.fillna(value={int_col: getattr(args, int_col)}, inplace=True)
                batch_params[int_col] = batch_params[int_col].astype(int)

        # rename column "a" to "amplicon_seq", etc
        batch_params.rename(index=str, columns=CRISPRessoShared.get_crispresso_options_lookup("Core"), inplace=True)
        batch_count = batch_params.shape[0]
        batch_params.index = range(batch_count)

        if 'fastq_r1' not in batch_params and 'bam_input' not in batch_params:
            raise CRISPRessoShared.BadParameterException("fastq_r1 must be specified in the batch settings file. Current headings are: "
                    + str(batch_params.columns.values))

        # add args from the command line to batch_params_df
        for arg in vars(args):
            if arg not in batch_params:
                batch_params[arg] = getattr(args, arg)
            elif (getattr(args, arg) is not None):
                batch_params.fillna(value={arg: getattr(args, arg)}, inplace=True)

        # assert that all names are unique
        # and clean names

        for i in range(batch_count):
            if batch_params.loc[i, 'name'] == '':
                batch_params.at[i, 'name'] = i
            batch_params.at[i, 'name'] = CRISPRessoShared.clean_filename(batch_params.loc[i, 'name'])

        if batch_params.drop_duplicates('name').shape[0] != batch_params.shape[0]:
            raise CRISPRessoShared.BadParameterException('Batch input names must be unique. The given names are not unique: ' + str(batch_params.loc[:, 'name']))

        # for each row, check to make sure that run's parameters are correct
        for idx, row in batch_params.iterrows():
            # check parameters for batch input for each batch
            # Check presence of input fastq/bam files
            has_input = False
            if 'fastq_r1' in row and row.fastq_r1 != '':
                CRISPRessoShared.check_file(row.fastq_r1)
                has_input = True

            if 'fastq_r2' in row and row.fastq_r2 != '':
                CRISPRessoShared.check_file(row.fastq_r2)
                has_input = True

            if 'input_bam' in row and row.input_bam != '':
                CRISPRessoShared.check_file(row.input_bam)
                has_input = True

            if not has_input:
                raise CRISPRessoShared.BadParameterException("At least one fastq file must be given as a command line parameter or be specified in the batch settings file with the heading 'fastq_r1' (fastq_r1 on row %s '%s' is invalid)" % (int(idx) + 1, row.fastq_r1))

            if args.auto:
                continue

            curr_amplicon_seq_str = row.amplicon_seq
            if curr_amplicon_seq_str is None:
                raise CRISPRessoShared.BadParameterException("Amplicon sequence must be given as a command line parameter or be specified in the batch settings file with the heading 'amplicon_seq' (Amplicon seq on row %s '%s' is invalid)" % (int(idx) + 1, curr_amplicon_seq_str))

            # set quantification windows for each amplicon (needed to tell whether to discard a guide below)
            curr_amplicon_seq_arr = str(curr_amplicon_seq_str).split(',')
            curr_amplicon_quant_window_coordinates_arr = [None] * len(curr_amplicon_seq_arr)
            if row.quantification_window_coordinates is not None:
                for idx, coords in enumerate(row.quantification_window_coordinates.strip("'").strip('"').split(",")):
                    if coords != "":
                        curr_amplicon_quant_window_coordinates_arr[idx] = coords

            # assert that guides are in the amplicon sequences and that quantification windows are within the amplicon sequences
            guides_are_in_amplicon = {}  # dict of whether a guide is in at least one amplicon sequence
            # iterate through amplicons for this run
            for amp_idx, curr_amplicon_seq in enumerate(curr_amplicon_seq_arr):
                this_include_idxs = []  # mask for bp to include for this amplicon seq, as specified by sgRNA cut points
                this_sgRNA_intervals = []
                curr_amplicon_quant_window_coordinates = curr_amplicon_quant_window_coordinates_arr[amp_idx]
                wrong_nt = CRISPRessoShared.find_wrong_nt(curr_amplicon_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The amplicon sequence in row %d (%s) contains incorrect characters:%s' % (idx + 1, curr_amplicon_seq_str, ' '.join(wrong_nt)))

                # iterate through guides
                curr_guide_seq_string = row.guide_seq
                if curr_guide_seq_string is not None and re.match(r'(?i)^$|^nan?$', str(curr_guide_seq_string)) is None:
                    guides = str(curr_guide_seq_string).strip().upper().split(',')
                    for curr_guide_seq in guides:
                        wrong_nt = CRISPRessoShared.find_wrong_nt(curr_guide_seq)
                        if wrong_nt:
                            raise CRISPRessoShared.NTException('The sgRNA sequence in row %d (%s) contains incorrect characters:%s' % (idx + 1, curr_guide_seq, ' '.join(wrong_nt)))
                    guide_mismatches = [[]] * len(guides)
                    guide_names = [""] * len(guides)
                    guide_qw_centers = CRISPRessoShared.set_guide_array(row.quantification_window_center, guides, 'guide quantification center')
                    guide_qw_sizes = CRISPRessoShared.set_guide_array(row.quantification_window_size, guides, 'guide quantification size')
                    guide_plot_cut_points = [1] * len(guides)
                    discard_guide_positions_overhanging_amplicon_edge = False
                    if 'discard_guide_positions_overhanging_amplicon_edge' in row:
                        discard_guide_positions_overhanging_amplicon_edge = row.discard_guide_positions_overhanging_amplicon_edge
                    (this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_sgRNA_names, this_sgRNA_include_idxs, this_include_idxs,
                        this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(curr_amplicon_seq, guides, guide_mismatches, guide_names, guide_qw_centers,
                        guide_qw_sizes, curr_amplicon_quant_window_coordinates, row.exclude_bp_from_left, row.exclude_bp_from_right, row.plot_window_size, guide_plot_cut_points, discard_guide_positions_overhanging_amplicon_edge)
                    for guide_seq in this_sgRNA_sequences:
                        guides_are_in_amplicon[guide_seq] = 1

            for guide_seq in guides_are_in_amplicon:
                if guides_are_in_amplicon[guide_seq] != 1:
                    raise CRISPRessoShared.BadParameterException('The guide sequence provided on row %d (%s) is not present in any amplicon sequence:%s! \nNOTE: The guide will be ignored for the analysis. Please check your input!' % (idx + 1, row.guide_seq, curr_amplicon_seq))

        crispresso_cmds = []
        batch_names_arr = []
        batch_input_names = {}
        for idx, row in batch_params.iterrows():

            batch_name = CRISPRessoShared.slugify(row["name"])
            batch_names_arr.append(batch_name)
            batch_input_names[batch_name] = row["name"]

            crispresso_cmd = args.crispresso_command + ' -o "%s" --name %s' % (OUTPUT_DIRECTORY, batch_name)
            crispresso_cmd = CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, crispresso_options_for_batch, batch_params, idx)
            if re.match(r'(?i)^$|^nan?$', str(row.amplicon_seq)) is not None:
                crispresso_cmd += ' --auto '
                crispresso_cmd = re.sub(r'--amplicon_seq\s+[^ ]+\s*', '', crispresso_cmd)
            if re.match(r'(?i)^$|^nan?$', str(row.guide_seq)) is not None:
                crispresso_cmd = re.sub(r'--guide_seq\s+[^ ]+\s*', '', crispresso_cmd)
            crispresso_cmds.append(crispresso_cmd)

        crispresso2_info['results']['batch_names_arr'] = batch_names_arr
        crispresso2_info['results']['batch_input_names'] = batch_input_names
        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, n_processes_for_batch, 'batch', args.skip_failed, start_end_percent=[10, 90])

        run_datas = []  # crispresso2 info from each row

        all_amplicons = set()
        amplicon_names = {}
        amplicon_counts = {}
        completed_batch_arr = []
        failed_batch_arr = []
        failed_batch_arr_desc = []
        for idx, row in batch_params.iterrows():
            batch_name = CRISPRessoShared.slugify(row["name"])
            folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
            # check if run failed
            failed_run_bool, failed_status_string = CRISPRessoShared.check_if_failed_run(folder_name, info)
            if failed_run_bool:
                failed_batch_arr.append(batch_name)
                failed_batch_arr_desc.append(failed_status_string)
                run_datas.append(None)
                continue

            run_data = CRISPRessoShared.load_crispresso_info(folder_name)
            run_datas.append(run_data)
            for ref_name in run_data['results']['ref_names']:
                ref_seq = run_data['results']['refs'][ref_name]['sequence']
                all_amplicons.add(ref_seq)
                # if this amplicon is called something else in another sample, just call it the amplicon
                if ref_name in amplicon_names and amplicon_names[ref_seq] != ref_name:
                    amplicon_names[ref_seq] = ref_seq
                else:
                    amplicon_names[ref_seq] = ref_name
                if ref_seq not in amplicon_counts:
                    amplicon_counts[ref_seq] = 0
                amplicon_counts[ref_seq] += 1

            completed_batch_arr.append(batch_name)

        crispresso2_info['results']['failed_batch_arr'] = failed_batch_arr
        crispresso2_info['results']['failed_batch_arr_desc'] = failed_batch_arr_desc
        crispresso2_info['results']['completed_batch_arr'] = completed_batch_arr

        # make sure amplicon names aren't super long
        for amplicon in all_amplicons:
            if len(amplicon_names[amplicon]) > 21:
                amplicon_names[amplicon] = amplicon_names[amplicon][0:21]

        # make sure no duplicate names (same name for the different amplicons)
        seen_names = {}
        for amplicon in all_amplicons:
            suffix_counter = 2
            orig_name = amplicon_names[amplicon]
            while amplicon_names[amplicon] in seen_names:
                amplicon_names[amplicon] = orig_name + "_" + str(suffix_counter)
                suffix_counter += 1
            seen_names[amplicon_names[amplicon]] = 1

        save_png = True
        if args.suppress_report:
            save_png = False

        crispresso2_info['results']['general_plots']['summary_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_datas'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_names'] = []
        crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_paths'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_datas'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_divs'] = {}

        crispresso2_info['results']['general_plots']['allele_modification_line_plot_names'] = []
        crispresso2_info['results']['general_plots']['allele_modification_line_plot_paths'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_line_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_line_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_line_plot_datas'] = {}
        crispresso2_info['results']['general_plots']['allele_modification_line_plot_divs'] = {}

        # =====================================================================
        # Pass 1: Aggregation — accumulate per-amplicon data and write CSVs
        # No plotting calls in this pass.
        # =====================================================================
        all_nuc_freq_dfs = {}
        all_nuc_pct_dfs = {}
        all_mod_freq_dfs = {}
        all_mod_pct_dfs = {}
        all_consensus_guides = {}
        all_consensus_include_idxs = {}
        all_consensus_sgRNA_intervals = {}
        all_consensus_sgRNA_plot_idxs = {}
        all_guides_all_same = {}
        all_summary_filenames = {}
        batch_amplicon_names = []

        percent_complete_start, percent_complete_end = 90, 95
        if all_amplicons:
            percent_complete_step = (percent_complete_end - percent_complete_start) / len(all_amplicons)
        else:
            percent_complete_step = 0

        for amplicon_index, amplicon_seq in enumerate(all_amplicons):
            if amplicon_counts[amplicon_seq] < 1:
                continue

            amplicon_name = amplicon_names[amplicon_seq]
            percent_complete = percent_complete_start + (amplicon_index * percent_complete_step)
            info('Reporting summary for amplicon: "' + amplicon_name + '"', {'percent_complete': percent_complete})

            consensus_sequence = ""
            nucleotide_frequency_summary = []
            nucleotide_percentage_summary = []
            modification_frequency_summary = []
            modification_percentage_summary = []

            amp_found_count = 0
            consensus_guides = []
            consensus_include_idxs = []
            consensus_sgRNA_plot_idxs = []
            consensus_sgRNA_intervals = []
            guides_all_same = True
            batches_with_this_amplicon = []
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
                run_data = run_datas[idx]
                if run_data is None:
                    continue
                batch_has_amplicon = False
                batch_amplicon_name = ''
                for ref_name in run_data['results']['ref_names']:
                    if amplicon_seq == run_data['results']['refs'][ref_name]['sequence']:
                        batch_has_amplicon = True
                        batch_amplicon_name = ref_name
                if not batch_has_amplicon:
                    continue
                batches_with_this_amplicon.append(idx)

                if consensus_guides == []:
                    consensus_guides = run_data['results']['refs'][batch_amplicon_name]['sgRNA_sequences']
                    consensus_include_idxs = run_data['results']['refs'][batch_amplicon_name]['include_idxs']
                    consensus_sgRNA_intervals = run_data['results']['refs'][batch_amplicon_name]['sgRNA_intervals']
                    consensus_sgRNA_plot_idxs = run_data['results']['refs'][batch_amplicon_name]['sgRNA_plot_idxs']

                if run_data['results']['refs'][batch_amplicon_name]['sgRNA_sequences'] != consensus_guides:
                    guides_all_same = False
                if set(run_data['results']['refs'][batch_amplicon_name]['include_idxs']) != set(consensus_include_idxs):
                    guides_all_same = False

                if 'nuc_freq_filename' not in run_data['results']['refs'][batch_amplicon_name]:
                    info("Skipping the amplicon '%s' in folder '%s'. Cannot find nucleotide information." % (batch_amplicon_name, folder_name))
                    continue

                nucleotide_frequency_file = os.path.join(folder_name, run_data['results']['refs'][batch_amplicon_name]['nuc_freq_filename'])
                ampSeq_nf, nuc_freqs = CRISPRessoShared.parse_count_file(nucleotide_frequency_file)

                nucleotide_pct_file = os.path.join(folder_name, run_data['results']['refs'][batch_amplicon_name]['nuc_pct_filename'])
                ampSeq_np, nuc_pcts = CRISPRessoShared.parse_count_file(nucleotide_pct_file)

                count_file = os.path.join(folder_name, run_data['results']['refs'][batch_amplicon_name]['mod_count_filename'])
                ampSeq_cf, mod_freqs = CRISPRessoShared.parse_count_file(count_file)

                if ampSeq_nf is None or ampSeq_np is None or ampSeq_cf is None:
                    info("Skipping the amplicon '%s' in folder '%s'. Could not parse batch output." % (batch_amplicon_name, folder_name))
                    info("Nucleotide frequency amplicon: '%s', Nucleotide percentage amplicon: '%s', Count vectors amplicon: '%s'" % (ampSeq_nf, ampSeq_np, ampSeq_cf))
                    continue
                if ampSeq_nf != ampSeq_np or ampSeq_np != ampSeq_cf:
                    warn("Skipping the amplicon '%s' in folder '%s'. Parsed amplicon sequences do not match\nnf:%s\nnp:%s\ncf:%s\nrf:%s" % (batch_amplicon_name, folder_name, ampSeq_nf, ampSeq_np, ampSeq_cf, amplicon_seq))
                    continue
                if consensus_sequence == "":
                    consensus_sequence = ampSeq_nf
                if ampSeq_nf != consensus_sequence:
                    info("Skipping the amplicon '%s' in folder '%s'. Amplicon sequences do not match." % (batch_amplicon_name, folder_name))
                    continue
                if 'Total' not in mod_freqs:
                    info("Skipping the amplicon '%s' in folder '%s'. Processing did not complete." % (batch_amplicon_name, folder_name))
                    continue
                if mod_freqs['Total'][0] == 0 or mod_freqs['Total'][0] == "0":
                    info("Skipping the amplicon '%s' in folder '%s'. Got no reads for amplicon." % (batch_amplicon_name, folder_name))
                    continue
                this_amp_total_reads = run_data['results']['alignment_stats']['counts_total'][batch_amplicon_name]
                if this_amp_total_reads < args.min_reads_for_inclusion:
                    info("Skipping the amplicon '%s' in folder '%s'. Got %s reads (min_reads_for_inclusion is %d)." % (batch_amplicon_name, folder_name, str(this_amp_total_reads), args.min_reads_for_inclusion))
                    continue

                mod_pcts = {}
                for key in mod_freqs:
                    mod_pcts[key] = np.array(mod_freqs[key]).astype(float) / float(mod_freqs['Total'][0])

                amp_found_count += 1

                for nuc in ['A', 'T', 'C', 'G', 'N', '-']:
                    row = [batch_name, nuc]
                    row.extend(nuc_freqs[nuc])
                    nucleotide_frequency_summary.append(row)

                    pct_row = [batch_name, nuc]
                    pct_row.extend(nuc_pcts[nuc])
                    nucleotide_percentage_summary.append(pct_row)

                for mod in ['Insertions', 'Insertions_Left', 'Deletions', 'Substitutions', 'All_modifications']:
                    row = [batch_name, mod]
                    row.extend(mod_freqs[mod])
                    modification_frequency_summary.append(row)

                    pct_row = [batch_name, mod]
                    pct_row.extend(mod_pcts[mod])
                    modification_percentage_summary.append(pct_row)

            if amp_found_count == 0:
                info("Couldn't find any data for amplicon '%s'. Not compiling results." % amplicon_name)
                continue

            amplicon_plot_name = amplicon_name + "."
            if len(amplicon_names) == 1 and amplicon_name == "Reference":
                amplicon_plot_name = ""

            # Build summary DataFrames
            colnames = ['Batch', 'Nucleotide']
            colnames.extend(list(consensus_sequence))
            nucleotide_frequency_summary_df = pd.DataFrame(nucleotide_frequency_summary, columns=colnames)
            nucleotide_frequency_summary_df = pd.concat([nucleotide_frequency_summary_df.iloc[:, 0:2],
                                                        nucleotide_frequency_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
            nucleotide_frequency_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_frequency_summary.txt')
            nucleotide_frequency_summary_df.to_csv(nucleotide_frequency_summary_filename, sep='\t', index=None)

            nucleotide_percentage_summary_df = pd.DataFrame(nucleotide_percentage_summary, columns=colnames)
            nucleotide_percentage_summary_df = pd.concat([nucleotide_percentage_summary_df.iloc[:, 0:2],
                                                    nucleotide_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
            nucleotide_percentage_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_percentage_summary.txt')
            nucleotide_percentage_summary_df.to_csv(nucleotide_percentage_summary_filename, sep='\t', index=None)

            colnames = ['Batch', 'Modification']
            colnames.extend(list(consensus_sequence))
            modification_frequency_summary_df = pd.DataFrame(modification_frequency_summary, columns=colnames)
            modification_frequency_summary_df = pd.concat([modification_frequency_summary_df.iloc[:, 0:2],
                                                        modification_frequency_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
            modification_frequency_summary_filename = _jp(amplicon_plot_name + 'MODIFICATION_FREQUENCY_SUMMARY.txt')
            modification_frequency_summary_df.to_csv(modification_frequency_summary_filename, sep='\t', index=None)

            modification_percentage_summary_df = pd.DataFrame(modification_percentage_summary, columns=colnames)
            modification_percentage_summary_df = pd.concat([modification_percentage_summary_df.iloc[:, 0:2],
                                                    modification_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
            modification_percentage_summary_filename = _jp(amplicon_plot_name + 'MODIFICATION_PERCENTAGE_SUMMARY.txt')
            modification_percentage_summary_df.to_csv(modification_percentage_summary_filename, sep='\t', index=None)

            crispresso2_info['results']['nucleotide_frequency_summary_filename'] = os.path.basename(nucleotide_frequency_summary_filename)
            crispresso2_info['results']['nucleotide_percentage_summary_filename'] = os.path.basename(nucleotide_percentage_summary_filename)
            crispresso2_info['results']['modification_frequency_summary_filename'] = os.path.basename(modification_frequency_summary_filename)
            crispresso2_info['results']['modification_percentage_summary_filename'] = os.path.basename(modification_percentage_summary_filename)

            # Per-sgRNA CSV writes (data output, not plotting)
            if guides_all_same and consensus_guides != []:
                info("All guides are equal. Performing comparison of batches for amplicon '%s'" % amplicon_name)
                for idx, sgRNA in enumerate(consensus_guides):
                    sgRNA_plot_idxs = consensus_sgRNA_plot_idxs[idx]
                    plot_idxs_flat = [0, 1]
                    plot_idxs_flat.extend([plot_idx + 2 for plot_idx in sgRNA_plot_idxs])

                    sub_nucleotide_frequency_summary_df = nucleotide_frequency_summary_df.iloc[:, plot_idxs_flat]
                    sub_nucleotide_frequency_summary_df = pd.concat([sub_nucleotide_frequency_summary_df.iloc[:, 0:2],
                                                                sub_nucleotide_frequency_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                    sub_nucleotide_frequency_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_frequency_summary_around_sgRNA_' + sgRNA + '.txt')
                    sub_nucleotide_frequency_summary_df.to_csv(sub_nucleotide_frequency_summary_filename, sep='\t', index=None)

                    sub_nucleotide_percentage_summary_df = nucleotide_percentage_summary_df.iloc[:, plot_idxs_flat]
                    sub_nucleotide_percentage_summary_df = pd.concat([sub_nucleotide_percentage_summary_df.iloc[:, 0:2],
                                                                sub_nucleotide_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                    sub_nucleotide_percentage_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_percentage_summary_around_sgRNA_' + sgRNA + '.txt')
                    sub_nucleotide_percentage_summary_df.to_csv(sub_nucleotide_percentage_summary_filename, sep='\t', index=None)

                    sub_modification_percentage_summary_df = modification_percentage_summary_df.iloc[:, plot_idxs_flat]
                    sub_modification_percentage_summary_df = pd.concat([sub_modification_percentage_summary_df.iloc[:, 0:2],
                                                                sub_modification_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                    sub_modification_percentage_summary_filename = _jp(amplicon_plot_name + 'Modification_percentage_summary_around_sgRNA_' + sgRNA + '.txt')
                    sub_modification_percentage_summary_df.to_csv(sub_modification_percentage_summary_filename, sep='\t', index=None)

            # Store per-amplicon data for context construction
            batch_amplicon_names.append(amplicon_name)
            all_nuc_freq_dfs[amplicon_name] = nucleotide_frequency_summary_df
            all_nuc_pct_dfs[amplicon_name] = nucleotide_percentage_summary_df
            all_mod_freq_dfs[amplicon_name] = modification_frequency_summary_df
            all_mod_pct_dfs[amplicon_name] = modification_percentage_summary_df
            all_consensus_guides[amplicon_name] = consensus_guides
            all_consensus_include_idxs[amplicon_name] = consensus_include_idxs
            all_consensus_sgRNA_intervals[amplicon_name] = consensus_sgRNA_intervals
            all_consensus_sgRNA_plot_idxs[amplicon_name] = consensus_sgRNA_plot_idxs
            all_guides_all_same[amplicon_name] = guides_all_same
            all_summary_filenames[amplicon_name] = {
                'nucleotide_frequency': nucleotide_frequency_summary_filename,
                'modification_frequency': modification_frequency_summary_filename,
            }
        # end pass 1

        # =====================================================================
        # Pass 2: Context construction + plotting
        # =====================================================================
        from CRISPResso2.plots.plot_context import BatchPlotContext

        batch_plot_context = BatchPlotContext(
            args=args,
            run_data=crispresso2_info,
            output_directory=OUTPUT_DIRECTORY,
            save_png=save_png,
            _jp=_jp,
            custom_config=custom_config,
            amplicon_names=batch_amplicon_names,
            nucleotide_frequency_summary_dfs=all_nuc_freq_dfs,
            nucleotide_percentage_summary_dfs=all_nuc_pct_dfs,
            modification_frequency_summary_dfs=all_mod_freq_dfs,
            modification_percentage_summary_dfs=all_mod_pct_dfs,
            consensus_guides=all_consensus_guides,
            consensus_include_idxs=all_consensus_include_idxs,
            consensus_sgRNA_intervals=all_consensus_sgRNA_intervals,
            consensus_sgRNA_plot_idxs=all_consensus_sgRNA_plot_idxs,
            guides_all_same=all_guides_all_same,
            all_summary_filenames=all_summary_filenames,
            sub_nucleotide_frequency_summary_filename=locals().get(
                'sub_nucleotide_frequency_summary_filename', ''
            ),
            sub_nucleotide_percentage_summary_filename=locals().get(
                'sub_nucleotide_percentage_summary_filename', ''
            ),
        )
        # Stash n_processes_for_batch on args so Pro's plot_runners
        # (via the on_batch_plots_complete hook) can read it.  CORE's
        # elif branch below uses the local n_processes_for_batch
        # directly.
        args.n_processes_for_batch = n_processes_for_batch

        pro_plots_ran = C2PRO_INSTALLED and CRISPRessoShared.run_C2Pro_hook('on_batch_plots_complete', batch_plot_context, logger)
        if not pro_plots_ran and not args.suppress_plots and not args.suppress_batch_summary_plots:
            # Built-in matplotlib plot iteration for the non-Pro path.
            n_processes = int(n_processes_for_batch)

            if n_processes > 1:
                process_pool = ProcessPoolExecutor(n_processes)
                process_futures = {}
            else:
                process_pool = None
                process_futures = None

            plot = partial(
                CRISPRessoMultiProcessing.run_plot,
                num_processes=n_processes,
                process_futures=process_futures,
                process_pool=process_pool,
                halt_on_plot_fail=args.halt_on_plot_fail,
            )

            general_plots = crispresso2_info['results']['general_plots']
            general_plots.setdefault('summary_plot_names', [])
            general_plots.setdefault('summary_plot_titles', {})
            general_plots.setdefault('summary_plot_labels', {})
            general_plots.setdefault('summary_plot_datas', {})

            window_nuc_pct_quilt_plot_names = []
            nuc_pct_quilt_plot_names = []
            window_nuc_conv_plot_names = []
            nuc_conv_plot_names = []

            for amplicon_name in batch_plot_context.amplicon_names:
                batch_plot_context.amplicon_name = amplicon_name
                nuc_freq_filename = batch_plot_context.all_summary_filenames[amplicon_name]['nucleotide_frequency']
                mod_freq_filename = batch_plot_context.all_summary_filenames[amplicon_name]['modification_frequency']
                guides_all_same = batch_plot_context.guides_all_same[amplicon_name]
                consensus_guides = batch_plot_context.consensus_guides[amplicon_name]

                if guides_all_same and consensus_guides:
                    for sgRNA_ind, sgRNA in enumerate(consensus_guides):
                        batch_plot_context.sgRNA_ind = sgRNA_ind

                        if should_plot_large_plots(
                            batch_plot_context.nucleotide_percentage_summary_dfs[amplicon_name].shape[0],
                            False, False,
                        ):
                            quilt_input = prep_batch_nuc_quilt_around_sgRNA(batch_plot_context)
                            debug(f'Plotting nucleotide percentage quilt for amplicon {amplicon_name}, sgRNA {sgRNA}')
                            plot(CRISPRessoPlot.plot_nucleotide_quilt, quilt_input)
                            plot_name = os.path.basename(quilt_input['fig_filename_root'])
                            window_nuc_pct_quilt_plot_names.append(plot_name)
                            general_plots['summary_plot_titles'][plot_name] = f'sgRNA: {sgRNA} Amplicon: {amplicon_name}'
                            general_plots['summary_plot_labels'][plot_name] = f'Composition of each base around the guide {sgRNA} for the amplicon {amplicon_name}'
                            general_plots['summary_plot_datas'][plot_name] = [
                                ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                                ('Modification frequencies', os.path.basename(mod_freq_filename)),
                            ]

                            if args.base_editor_output:
                                conv_input = prep_batch_conversion_map_around_sgRNA(batch_plot_context)
                                debug(f'Plotting nucleotide conversion map for amplicon {amplicon_name}, sgRNA {sgRNA}')
                                plot(CRISPRessoPlot.plot_conversion_map, conv_input)
                                plot_name = os.path.basename(conv_input['fig_filename_root'])
                                window_nuc_conv_plot_names.append(plot_name)
                                title = f'sgRNA: {sgRNA} Amplicon: {amplicon_name}'
                                if len(consensus_guides) == 1:
                                    title = ''
                                general_plots['summary_plot_titles'][plot_name] = title
                                general_plots['summary_plot_labels'][plot_name] = (
                                    f'{args.conversion_nuc_from}->{args.conversion_nuc_to} '
                                    f'conversion rates around the guide {sgRNA} for the amplicon {amplicon_name}'
                                )
                                general_plots['summary_plot_datas'][plot_name] = [
                                    ('Nucleotide frequencies around sgRNA', os.path.basename(batch_plot_context.sub_nucleotide_frequency_summary_filename)),
                                    ('Nucleotide percentages around sgRNA', os.path.basename(batch_plot_context.sub_nucleotide_percentage_summary_filename)),
                                ]

                    # Whole-region quilt (when guides_all_same)
                    if should_plot_large_plots(
                        batch_plot_context.nucleotide_percentage_summary_dfs[amplicon_name].shape[0],
                        False, False,
                    ):
                        quilt_input = prep_batch_nuc_quilt(batch_plot_context)
                        debug(f'Plotting nucleotide quilt for {amplicon_name}')
                        plot(CRISPRessoPlot.plot_nucleotide_quilt, quilt_input)
                        plot_name = os.path.basename(quilt_input['fig_filename_root'])
                        nuc_pct_quilt_plot_names.append(plot_name)
                        title = f'Amplicon: {amplicon_name}'
                        if len(batch_plot_context.amplicon_names) == 1:
                            title = ''
                        general_plots['summary_plot_titles'][plot_name] = title
                        general_plots['summary_plot_labels'][plot_name] = f'Composition of each base for the amplicon {amplicon_name}'
                        general_plots['summary_plot_datas'][plot_name] = [
                            ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                            ('Modification frequencies', os.path.basename(mod_freq_filename)),
                        ]

                        if args.base_editor_output:
                            conv_input = prep_batch_conversion_map(batch_plot_context)
                            debug(f'Plotting nucleotide conversion map for {amplicon_name}')
                            plot(CRISPRessoPlot.plot_conversion_map, conv_input)
                            plot_name = os.path.basename(conv_input['fig_filename_root'])
                            nuc_conv_plot_names.append(plot_name)
                            general_plots['summary_plot_titles'][plot_name] = ''
                            general_plots['summary_plot_labels'][plot_name] = (
                                f'{args.conversion_nuc_from}->{args.conversion_nuc_to} '
                                f'conversion rates for the amplicon {amplicon_name}'
                            )
                            general_plots['summary_plot_datas'][plot_name] = [
                                ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                                ('Modification frequencies', os.path.basename(mod_freq_filename)),
                            ]

                # Guides not all same — whole-region quilt only, no sgRNA data on plot
                elif should_plot_large_plots(
                    batch_plot_context.nucleotide_percentage_summary_dfs[amplicon_name].shape[0],
                    False, False,
                ):
                    quilt_input = prep_batch_nuc_quilt(batch_plot_context)
                    debug(f'Plotting nucleotide quilt for {amplicon_name}')
                    plot(CRISPRessoPlot.plot_nucleotide_quilt, quilt_input)
                    plot_name = os.path.basename(quilt_input['fig_filename_root'])
                    nuc_pct_quilt_plot_names.append(plot_name)
                    general_plots['summary_plot_labels'][plot_name] = f'Composition of each base for the amplicon {amplicon_name}'
                    general_plots['summary_plot_datas'][plot_name] = [
                        ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                        ('Modification frequencies', os.path.basename(mod_freq_filename)),
                    ]

                    if args.base_editor_output:
                        conv_input = prep_batch_conversion_map(batch_plot_context)
                        debug(f'Plotting BE nucleotide conversion map for {amplicon_name}')
                        plot(CRISPRessoPlot.plot_conversion_map, conv_input)
                        plot_name = os.path.basename(conv_input['fig_filename_root'])
                        nuc_conv_plot_names.append(plot_name)
                        general_plots['summary_plot_labels'][plot_name] = (
                            f'{args.conversion_nuc_from}->{args.conversion_nuc_to} '
                            f'conversion rates for the amplicon {amplicon_name}'
                        )
                        general_plots['summary_plot_datas'][plot_name] = [
                            ('Nucleotide frequencies', os.path.basename(nuc_freq_filename)),
                            ('Modification frequencies', os.path.basename(mod_freq_filename)),
                        ]

            general_plots['window_nuc_pct_quilt_plot_names'] = window_nuc_pct_quilt_plot_names
            general_plots['nuc_pct_quilt_plot_names'] = nuc_pct_quilt_plot_names
            general_plots['window_nuc_conv_plot_names'] = window_nuc_conv_plot_names
            general_plots['nuc_conv_plot_names'] = nuc_conv_plot_names

            if process_pool is not None:
                wait(process_futures)
                for future in process_futures:
                    try:
                        future.result()
                    except Exception as e:
                        warn(f'Error in plot pool: {e}')
                        debug(traceback.format_exc())
                process_pool.shutdown()

        # summarize amplicon modifications
        with open(_jp('CRISPRessoBatch_quantification_of_editing_frequency.txt'), 'w') as outfile:
            wrote_header = False
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
                run_data = run_datas[idx]
                if run_data is None:
                    continue

                amplicon_modification_file = os.path.join(folder_name, run_data['running_info']['quant_of_editing_freq_filename'])
                with open(amplicon_modification_file, 'r') as infile:
                    file_head = infile.readline()
                    if not wrote_header:
                        outfile.write('Batch\t' + file_head)
                        wrote_header = True
                    for line in infile:
                        outfile.write(batch_name + "\t" + line)

        # summarize frameshift and splicing site mods
        with open(_jp('CRISPRessoBatch_quantification_of_frameshift_splicing.txt'), 'w') as outfile:
            outfile.write("\t".join(['Batch', 'Reference', 'total_reads', 'modified_frameshift', 'modified_non_frameshift', 'non_modified_non_frameshift', 'splicing_sites_modified', 'splice_sites_unmodified']) + "\n")
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
                run_data = run_datas[idx]
                if run_data is None:
                    continue

                for ref_name in run_data['results']['ref_names']:
                    count_total = run_data['results']['alignment_stats']['counts_total'][ref_name]
                    count_modified_frameshift = run_data['results']['alignment_stats']['counts_modified_frameshift'][ref_name]
                    count_modified_non_frameshift = run_data['results']['alignment_stats']['counts_modified_non_frameshift'][ref_name]
                    count_non_modified_non_frameshift = run_data['results']['alignment_stats']['counts_non_modified_non_frameshift'][ref_name]
                    count_splicing_sites_modified = run_data['results']['alignment_stats']['counts_splicing_sites_modified'][ref_name]
                    count_splicing_sites_unmodified = count_total - count_splicing_sites_modified
                    outfile.write("\t".join([str(x) for x in [batch_name, ref_name, count_total, count_modified_frameshift, count_modified_non_frameshift, count_non_modified_non_frameshift, count_splicing_sites_modified, count_splicing_sites_unmodified]]) + "\n")

        # summarize alignment
        with open(_jp('CRISPRessoBatch_mapping_statistics.txt'), 'w') as outfile:
            wrote_header = False
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)

                run_data = run_datas[idx]
                if run_data is None:
                    continue
                amplicon_modification_file = os.path.join(folder_name, run_data['running_info']['mapping_stats_filename'])
                with open(amplicon_modification_file, 'r') as infile:
                    file_head = infile.readline()
                    if not wrote_header:
                        outfile.write('Batch\t' + file_head)
                        wrote_header = True
                    for line in infile:
                        outfile.write(batch_name + "\t" + line)

        # (Built-in plot multiprocessing teardown is handled inside
        # ``run_builtin_batch_plots``; nothing to do here.)

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Batch_report.html")
            else:
                report_name = OUTPUT_DIRECTORY + '.html'
            pro_report = CRISPRessoShared.get_C2Pro_hook('make_batch_report') if C2PRO_INSTALLED else None
            if pro_report:
                pro_report(crispresso2_info, report_name, OUTPUT_DIRECTORY, _ROOT, logger, batch_plot_context)
            else:
                CRISPRessoReport.make_batch_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT, logger)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        end_time = datetime.now()
        end_time_string = end_time.strftime('%Y-%m-%d %H:%M:%S')
        running_time = end_time - start_time
        running_time_string = str(running_time)

        crispresso2_info['running_info']['start_time'] = start_time
        crispresso2_info['running_info']['start_time_string'] = start_time_string
        crispresso2_info['running_info']['end_time'] = end_time
        crispresso2_info['running_info']['end_time_string'] = end_time_string
        crispresso2_info['running_info']['running_time'] = running_time
        crispresso2_info['running_info']['running_time_string'] = running_time_string

        CRISPRessoShared.write_crispresso_info(
            crispresso2Batch_info_file,
            crispresso2_info,
        )
        info('Analysis Complete!', {'percent_complete': 100})
        if args.zip_output:
            if args.output_folder == "":
                path_value = os.path.split(OUTPUT_DIRECTORY)
                CRISPRessoShared.zip_results(path_value[1])
            else:
                CRISPRessoShared.zip_results(OUTPUT_DIRECTORY)
        info(CRISPRessoShared.get_crispresso_footer())
        sys.exit(0)

    except Exception as e:
        debug_flag = False
        if 'args' in vars() and 'debug' in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)

        error('\n\nERROR: %s' % e)
        sys.exit(-1)


if __name__ == '__main__':
    main()
