# -*- coding: utf-8 -*-
"""CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
"""
import os
from copy import deepcopy
import sys
import traceback
from CRISPResso2 import CRISPRessoShared
from CRISPResso2.CRISPRessoReports import CRISPRessoReport
from CRISPResso2.plots.data_prep import (
    prep_compare_allele_table,
    prep_compare_editing_barchart,
    prep_compare_modification_positions,
)

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error = logger.critical
warn = logger.warning
debug = logger.debug
info = logger.info


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoCompare!' % library_name)
                sys.exit(1)


def parse_profile(profile_file):
    return np.loadtxt(profile_file, skiprows=1)


# EXCEPTIONS############################

class MixedRunningModeException(Exception):
    pass


class DifferentAmpliconLengthException(Exception):
    pass
############################


np = check_library('numpy')
pd = check_library('pandas')


_ROOT = os.path.abspath(os.path.dirname(__file__))


def normalize_name(name, output_folder_1, output_folder_2):
    get_name_from_folder = lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPResso_on_', '')
    if not name:
        return '{0}_VS_{1}'.format(
            get_name_from_folder(output_folder_1),
            get_name_from_folder(output_folder_2),
        )
    else:
        return name


def get_matching_allele_files(run_info_1, run_info_2):
    def get_amplicon_info(run_info):
        return {
            amplicon['sequence']: {
                'name': amplicon_name,
                'guides': amplicon['sgRNA_orig_sequences'],
                'cut_points': amplicon['sgRNA_cut_points'],
                'allele_files': amplicon['allele_frequency_files'],
            }
            for amplicon_name, amplicon in run_info['results']['refs'].items()
        }
    amplicons_1 = get_amplicon_info(run_info_1)
    amplicons_2 = get_amplicon_info(run_info_2)
    matching_allele_files = []
    for sequence_1 in amplicons_1:
        if sequence_1 in amplicons_2:
            if amplicons_1[sequence_1]['cut_points'] != amplicons_2[sequence_1]['cut_points']:
                warn(f'Report 1 has different cut points than report 2 for amplicon {amplicons_1[sequence_1]["name"]}, skipping comparison')
                continue
            guides_1 = set(amplicons_1[sequence_1]['guides'])
            guides_2 = set(amplicons_2[sequence_1]['guides'])
            if not guides_1 & guides_2:
                warn(f'Report 1 has no shared guides with report 2 for amplicon {amplicons_1[sequence_1]["name"]}, skipping comparison')
                continue
            matching_allele_files.extend((f_1, f_2) for f_1, f_2 in zip(amplicons_1[sequence_1]['allele_files'], amplicons_2[sequence_1]['allele_files']))

    return matching_allele_files


def main():
    try:
        parser = CRISPRessoShared.getCRISPRessoArgParser("Compare", parser_title='CRISPRessoCompare Parameters')

        args = parser.parse_args()

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)
        # Validate Pro-only inline config input early.
        CRISPRessoShared.check_custom_config(args)

        description = ['~~~CRISPRessoCompare~~~', '-Comparison of two CRISPResso analyses-']
        compare_header = r'''
 ___________________________
| __ __      __      __  __ |
|/  /  \|\/||__) /\ |__)|_  |
|\__\__/|  ||   /--\| \ |__ |
|___________________________|
        '''
        compare_header = CRISPRessoShared.get_crispresso_header(description, compare_header)
        info(compare_header)

        # CORE always uses matplotlib; when Pro is installed the hook
        # below skips this path entirely and Pro owns plotting decisions.
        from CRISPResso2.plots import CRISPRessoPlot

        if args.zip_output and not args.place_report_in_output_folder:
            warn('Invalid argument combination: If zip_output is True then place_report_in_output_folder must also be True. Setting place_report_in_output_folder to True.')
            args.place_report_in_output_folder = True
        # check that the CRISPResso output is present and fill amplicon_info
        quantification_file_1, amplicon_names_1, amplicon_info_1 = CRISPRessoShared.check_output_folder(args.crispresso_output_folder_1)
        quantification_file_2, amplicon_names_2, amplicon_info_2 = CRISPRessoShared.check_output_folder(args.crispresso_output_folder_2)

        run_info_1 = CRISPRessoShared.load_crispresso_info(args.crispresso_output_folder_1)

        run_info_2 = CRISPRessoShared.load_crispresso_info(args.crispresso_output_folder_2)

        sample_1_name = args.sample_1_name
        if args.sample_1_name is None:
            sample_1_name = "Sample 1"
            if 'running_info' in run_info_1 and 'name' in run_info_1['running_info'] and run_info_1['running_info']['name']:
                sample_1_name = run_info_1['running_info']['name']

        sample_2_name = args.sample_2_name
        if args.sample_2_name is None:
            sample_2_name = "Sample 2"
            if 'running_info' in run_info_2 and 'name' in run_info_2['running_info'] and run_info_2['running_info']['name']:
                sample_2_name = run_info_2['running_info']['name']

        if sample_1_name == sample_2_name:
            sample_2_name += '_2'

        OUTPUT_DIRECTORY = 'CRISPRessoCompare_on_{0}'.format(normalize_name(
            args.name, args.crispresso_output_folder_1, args.crispresso_output_folder_2,
        ))

        if args.output_folder:
            OUTPUT_DIRECTORY = os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename)  # handy function to put a file in the output directory
        log_filename = _jp('CRISPRessoCompare_RUNNING_LOG.txt')

        try:
            info('Creating Folder %s' % OUTPUT_DIRECTORY, {'percent_complete': 0})
            os.makedirs(OUTPUT_DIRECTORY)
            info('Done!')
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename = _jp('CRISPRessoCompare_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        logger.addHandler(CRISPRessoShared.StatusHandler(os.path.join(OUTPUT_DIRECTORY, 'CRISPRessoCompare_status.json')))

        with open(log_filename, 'w+') as outfile:
            outfile.write('[Command used]:\nCRISPRessoCompare %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Compare_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2Compare_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}}  # keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        crispresso2_info['results']['general_plots']['summary_plot_names'] = []
        crispresso2_info['results']['general_plots']['summary_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_datas'] = {}

        save_png = True
        if args.suppress_report:
            save_png = False

        # LOAD DATA
        amplicon_names_in_both = [amplicon_name for amplicon_name in amplicon_names_1 if amplicon_name in amplicon_names_2]

        # --- Build ComparePlotContext ---
        from CRISPResso2.plots.plot_context import ComparePlotContext

        plot_context = ComparePlotContext(
            args=args,
            run_data=crispresso2_info,
            output_directory=OUTPUT_DIRECTORY,
            save_png=save_png,
            _jp=_jp,
            custom_config={},
            amplicon_names=amplicon_names_in_both,
            sample_1_name=sample_1_name,
            sample_2_name=sample_2_name,
            run_info_1=run_info_1,
            run_info_2=run_info_2,
            amplicon_info_1=amplicon_info_1,
            amplicon_info_2=amplicon_info_2,
        )

        C2PRO_INSTALLED = CRISPRessoShared.is_C2Pro_installed()
        pro_plots_ran = C2PRO_INSTALLED and CRISPRessoShared.run_C2Pro_hook('on_compare_plots_complete', plot_context, logger)
        if not pro_plots_ran:
            # Built-in matplotlib plot iteration for the non-Pro path.
            general_plots = crispresso2_info['results']['general_plots']
            general_plots.setdefault('summary_plot_names', [])
            general_plots.setdefault('summary_plot_titles', {})
            general_plots.setdefault('summary_plot_labels', {})
            general_plots.setdefault('summary_plot_datas', {})

            sig_counts = {}
            sig_counts_quant_window = {}

            percent_complete_start, percent_complete_end = 10, 90
            if amplicon_names_in_both:
                percent_complete_step = (
                    (percent_complete_end - percent_complete_start)
                    / len(amplicon_names_in_both)
                )
            else:
                percent_complete_step = 0

            for amplicon_name in amplicon_names_in_both:
                percent_complete = (
                    percent_complete_start
                    + percent_complete_step * amplicon_names_in_both.index(amplicon_name)
                )
                info(
                    f'Loading data for amplicon {amplicon_name}',
                    extra={'percent_complete': percent_complete},
                )

                plot_context.amplicon_name = amplicon_name
                plot_context.profile_1 = parse_profile(
                    amplicon_info_1[amplicon_name]['quantification_file']
                )
                plot_context.profile_2 = parse_profile(
                    amplicon_info_2[amplicon_name]['quantification_file']
                )

                sig_counts[amplicon_name] = {}
                sig_counts_quant_window[amplicon_name] = {}

                try:
                    assert np.all(
                        plot_context.profile_1[:, 0]
                        == plot_context.profile_2[:, 0]
                    )
                except AssertionError:
                    raise DifferentAmpliconLengthException(
                        'Different amplicon lengths for the two amplicons.'
                    )

                plot_context.cut_points = (
                    run_info_1['results']['refs'][amplicon_name]['sgRNA_cut_points']
                )
                plot_context.sgRNA_intervals = (
                    run_info_1['results']['refs'][amplicon_name]['sgRNA_intervals']
                )

                # Plot 1: Editing comparison barchart
                barchart_input = prep_compare_editing_barchart(plot_context)
                CRISPRessoPlot.plot_quantification_comparison_barchart(**barchart_input)
                plot_name = os.path.basename(barchart_input['plot_path'])
                general_plots['summary_plot_names'].append(plot_name)
                general_plots['summary_plot_titles'][plot_name] = 'Editing efficiency comparison'
                general_plots['summary_plot_labels'][plot_name] = (
                    f'Figure 1: Comparison for amplicon {amplicon_name}; '
                    'Left: Percentage of modified and unmodified reads in each sample; '
                    'Right: relative percentage of modified and unmodified reads'
                )
                output_1 = os.path.join(
                    args.crispresso_output_folder_1,
                    run_info_1['running_info']['report_filename'],
                )
                output_2 = os.path.join(
                    args.crispresso_output_folder_2,
                    run_info_2['running_info']['report_filename'],
                )
                general_plots['summary_plot_datas'][plot_name] = []
                if os.path.isfile(output_1):
                    general_plots['summary_plot_datas'][plot_name].append(
                        (f'{sample_1_name} output', os.path.relpath(output_1, OUTPUT_DIRECTORY))
                    )
                if os.path.isfile(output_2):
                    general_plots['summary_plot_datas'][plot_name].append(
                        (f'{sample_2_name} output', os.path.relpath(output_2, OUTPUT_DIRECTORY))
                    )

                # Load modification count data
                mod_file_1 = amplicon_info_1[amplicon_name]['modification_count_file']
                amp_seq_1, mod_freqs_1 = CRISPRessoShared.parse_count_file(mod_file_1)
                mod_file_2 = amplicon_info_2[amplicon_name]['modification_count_file']
                amp_seq_2, mod_freqs_2 = CRISPRessoShared.parse_count_file(mod_file_2)
                if amp_seq_2 != amp_seq_1:
                    raise DifferentAmpliconLengthException(
                        'Different amplicon lengths for the two amplicons.'
                    )

                plot_context.mod_freqs_1 = mod_freqs_1
                plot_context.mod_freqs_2 = mod_freqs_2
                plot_context.consensus_sequence = amp_seq_1
                plot_context.quant_windows_1 = (
                    run_info_1['results']['refs'][amplicon_name]['include_idxs']
                )
                plot_context.quant_windows_2 = (
                    run_info_2['results']['refs'][amplicon_name]['include_idxs']
                )

                amplicon_plot_name = f'{amplicon_name}.'
                if len(amplicon_names_in_both) == 1 and amplicon_name == 'Reference':
                    amplicon_plot_name = ''

                # Plot 2: Modification positions (x4 mod types)
                for mod in ['Insertions', 'Deletions', 'Substitutions', 'All_modifications']:
                    plot_context.mod_type = mod
                    positions_data = prep_compare_modification_positions(plot_context)

                    mod_filename = _jp(f'{amplicon_plot_name}{mod}_quantification.txt')
                    positions_data['mod_df'].to_csv(mod_filename, sep='\t', index=None)

                    CRISPRessoPlot.plot_quantification_positions(**positions_data['plot_kwargs'])
                    plot_name = os.path.basename(positions_data['plot_kwargs']['plot_path'])
                    general_plots['summary_plot_names'].append(plot_name)
                    general_plots['summary_plot_titles'][plot_name] = (
                        f"{positions_data['mod_name']} locations"
                    )
                    general_plots['summary_plot_labels'][plot_name] = (
                        f"{positions_data['mod_name']} location comparison for amplicon "
                        f'{amplicon_name}; Top: percent difference; Bottom: p-value.'
                    )
                    general_plots['summary_plot_datas'][plot_name] = [
                        (f"{positions_data['mod_name']} quantification",
                         os.path.basename(mod_filename)),
                    ]

                    sig_counts[amplicon_name][mod] = positions_data['sig_count']
                    sig_counts_quant_window[amplicon_name][mod] = (
                        positions_data['sig_count_quant_window']
                    )

                # Plot 3: Allele table comparisons
                matching_allele_files = get_matching_allele_files(run_info_1, run_info_2)
                matching_allele_files.sort(key=lambda pair: 'base_edit' in pair[0])

                plot_context.allele_pairs = []
                for allele_file_1, allele_file_2 in matching_allele_files:
                    df1 = pd.read_csv(
                        os.path.join(args.crispresso_output_folder_1, allele_file_1),
                        sep='\t',
                    )
                    df2 = pd.read_csv(
                        os.path.join(args.crispresso_output_folder_2, allele_file_2),
                        sep='\t',
                    )
                    plot_context.allele_pairs.append(
                        (allele_file_1, allele_file_2, df1, df2)
                    )

                for pair_idx in range(len(plot_context.allele_pairs)):
                    allele_data = prep_compare_allele_table(plot_context, pair_idx)

                    allele_comparison_file = _jp(allele_data['file_root'] + '.txt')
                    allele_data['merged_df'].to_csv(
                        allele_comparison_file, sep='\t',
                    )

                    is_base_edit = allele_data['is_base_edit']
                    if is_base_edit:
                        title_prefix = 'Base edit comparison enriched in '
                        label_prefix = 'Base editing target nucleotide composition alleles.'
                    else:
                        title_prefix = 'Alleles enriched in '
                        label_prefix = 'Distribution comparison of alleles.'
                    label_suffix = (
                        ' Nucleotides are indicated by unique colors (A = green; '
                        'C = red; G = yellow; T = purple). Substitutions are shown in '
                        'bold font. Red rectangles highlight inserted sequences. '
                        'Horizontal dashed lines indicate deleted sequences. The '
                        'vertical dashed line indicates the predicted cleavage site. '
                        'The proportion and number of reads is shown for each sample '
                        f'on the right, with the values for {sample_1_name} followed '
                        f'by the values for {sample_2_name}.'
                    )

                    CRISPRessoPlot.plot_alleles_table_compare(**allele_data['plot_top_kwargs'])
                    plot_name = os.path.basename(allele_data['plot_top_kwargs']['fig_filename_root'])
                    general_plots['summary_plot_names'].append(plot_name)
                    general_plots['summary_plot_titles'][plot_name] = (
                        title_prefix + sample_1_name
                    )
                    general_plots['summary_plot_labels'][plot_name] = (
                        label_prefix + label_suffix
                        + f' Alleles are sorted for enrichment in {sample_1_name}.'
                    )
                    general_plots['summary_plot_datas'][plot_name] = [
                        ('Allele comparison table', os.path.basename(allele_comparison_file)),
                    ]

                    CRISPRessoPlot.plot_alleles_table_compare(**allele_data['plot_bottom_kwargs'])
                    plot_name = os.path.basename(allele_data['plot_bottom_kwargs']['fig_filename_root'])
                    general_plots['summary_plot_names'].append(plot_name)
                    general_plots['summary_plot_titles'][plot_name] = (
                        title_prefix + sample_2_name
                    )
                    general_plots['summary_plot_labels'][plot_name] = (
                        label_prefix + label_suffix
                        + f' Alleles are sorted for enrichment in {sample_2_name}.'
                    )
                    general_plots['summary_plot_datas'][plot_name] = [
                        ('Allele comparison table', os.path.basename(allele_comparison_file)),
                    ]

        # Pro's hook publishes Fisher/Bonferroni counts via
        # crispresso2_info['_sig_counts'] so the summary CSV below
        # can read them.  CORE's inline path sets local variables
        # directly (sig_counts / sig_counts_quant_window) in the
        # else-branch above; reuse them when available, otherwise
        # fall back to the Pro-published keys.
        if pro_plots_ran:
            sig_counts = crispresso2_info.get('_sig_counts', {})
            sig_counts_quant_window = crispresso2_info.get('_sig_counts_quant_window', {})

        debug('Calculating significant base counts...', {'percent_complete': 95})
        sig_counts_filename = _jp('CRISPRessoCompare_significant_base_counts.txt')
        with open(sig_counts_filename, 'w') as fout:
            fout.write('Amplicon\tModification\tsig_base_count\tsig_base_count_quant_window\n')
            for amplicon_name in amplicon_names_in_both:
                for mod in ['Insertions', 'Deletions', 'Substitutions', 'All_modifications']:
                    val = np.nan
                    if amplicon_name in sig_counts and mod in sig_counts[amplicon_name]:
                        val = sig_counts[amplicon_name][mod]

                    val_quant_window = np.nan
                    if amplicon_name in sig_counts_quant_window and mod in sig_counts_quant_window[amplicon_name]:
                        val_quant_window = sig_counts_quant_window[amplicon_name][mod]
                    line = "%s\t%s\t%s\t%s\n" % (amplicon_name, mod, val, val_quant_window)
                    fout.write(line)
        crispresso2_info['running_info']['sig_counts_report_location'] = sig_counts_filename

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Compare_report.html")
            else:
                report_name = OUTPUT_DIRECTORY + '.html'
            pro_report = CRISPRessoShared.get_C2Pro_hook('make_compare_report') if C2PRO_INSTALLED else None
            if pro_report:
                pro_report(crispresso2_info, report_name, OUTPUT_DIRECTORY, _ROOT, logger, plot_context)
            else:
                CRISPRessoReport.make_compare_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT, logger)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        CRISPRessoShared.write_crispresso_info(crispresso2Compare_info_file, crispresso2_info)

        if args.zip_output:
            CRISPRessoShared.zip_results(OUTPUT_DIRECTORY)

        info('Analysis Complete!', {'percent_complete': 100})
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
