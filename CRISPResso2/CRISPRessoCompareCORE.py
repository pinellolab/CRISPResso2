# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''
import os
from copy import deepcopy
import sys
import traceback
from CRISPResso2 import CRISPRessoShared
from CRISPResso2.CRISPRessoReports import CRISPRessoReport

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoCompare!' % library_name)
                sys.exit(1)


def parse_profile(profile_file):
    return np.loadtxt(profile_file, skiprows=1)


###EXCEPTIONS############################

class MixedRunningModeException(Exception):
    pass

class DifferentAmpliconLengthException(Exception):
    pass
############################

np=check_library('numpy')
pd=check_library('pandas')
import scipy.stats as stats


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
        description = ['~~~CRISPRessoCompare~~~', '-Comparison of two CRISPResso analyses-']
        compare_header = r'''
 ___________________________
| __ __      __      __  __ |
|/  /  \|\/||__) /\ |__)|_  |
|\__\__/|  ||   /--\| \ |__ |
|___________________________|
        '''
        compare_header = CRISPRessoShared.get_crispresso_header(description, compare_header)
        print(compare_header)

        parser = CRISPRessoShared.getCRISPRessoArgParser("Compare", parser_title = 'CRISPRessoCompare Parameters')

        args = parser.parse_args()

        if args.use_matplotlib or not CRISPRessoShared.is_C2Pro_installed():
            from CRISPResso2 import CRISPRessoPlot
        else:
            from CRISPRessoPro import plot as CRISPRessoPlot

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

        debug_flag = args.debug

        if args.zip_output and not args.place_report_in_output_folder:
            logger.warn('Invalid arguement combination: If zip_output is True then place_report_in_output_folder must also be True. Setting place_report_in_output_folder to True.')
            args.place_report_in_output_folder = True
        #check that the CRISPResso output is present and fill amplicon_info
        quantification_file_1, amplicon_names_1, amplicon_info_1=CRISPRessoShared.check_output_folder(args.crispresso_output_folder_1)
        quantification_file_2, amplicon_names_2, amplicon_info_2=CRISPRessoShared.check_output_folder(args.crispresso_output_folder_2)

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

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory
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
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be pickled and saved at the end of the run
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

        #LOAD DATA
        amplicon_names_in_both = [amplicon_name for amplicon_name in amplicon_names_1 if amplicon_name in amplicon_names_2]
        n_refs = len(amplicon_names_in_both)
        def get_plot_title_with_ref_name(plotTitle, refName):
            if n_refs > 1:
                return (plotTitle + ": " + refName)
            return plotTitle

        sig_counts = {}  # number of bp significantly modified (bonferonni corrected fisher pvalue)
        sig_counts_quant_window = {}
        percent_complete_start, percent_complete_end = 10, 90
        if amplicon_names_in_both:
            percent_complete_step = (percent_complete_end - percent_complete_start) / len(amplicon_names_in_both)
        for amplicon_name in amplicon_names_in_both:
            percent_complete = percent_complete_start + percent_complete_step * amplicon_names_in_both.index(amplicon_name)
            info('Loading data for amplicon %s' % amplicon_name, {'percent_complete': percent_complete})
            profile_1=parse_profile(amplicon_info_1[amplicon_name]['quantification_file'])
            profile_2=parse_profile(amplicon_info_2[amplicon_name]['quantification_file'])

            sig_counts[amplicon_name] = {}
            sig_counts_quant_window[amplicon_name] = {}

            amplicon_plot_name = amplicon_name+"."
            if len(amplicon_names_in_both) == 1 and amplicon_name == "Reference":
                amplicon_plot_name = ""

            try:
                assert np.all(profile_1[:, 0]==profile_2[:, 0])
            except:
                raise DifferentAmpliconLengthException('Different amplicon lengths for the two amplicons.')
            len_amplicon=profile_1.shape[0]
            effect_vector_any_1=profile_1[:, 1]
            effect_vector_any_2=profile_2[:, 1]

            cut_points = run_info_1['results']['refs'][amplicon_name]['sgRNA_cut_points']
            sgRNA_intervals = run_info_1['results']['refs'][amplicon_name]['sgRNA_intervals']


            #Quantification comparison barchart
            N_TOTAL_1 = float(amplicon_info_1[amplicon_name]['Reads_aligned'])
            N_UNMODIFIED_1 = float(amplicon_info_1[amplicon_name]['Unmodified'])
            N_MODIFIED_1 = float(amplicon_info_1[amplicon_name]['Modified'])

            N_TOTAL_2 = float(amplicon_info_2[amplicon_name]['Reads_aligned'])
            N_UNMODIFIED_2 = float(amplicon_info_2[amplicon_name]['Unmodified'])
            N_MODIFIED_2 = float(amplicon_info_2[amplicon_name]['Modified'])
            plot_titles = {
                'vs': get_plot_title_with_ref_name(
                    '%s VS %s' % (sample_1_name, sample_2_name),
                    amplicon_name,
                ),
                'diff': get_plot_title_with_ref_name(
                    '%s - %s' % (sample_1_name, sample_2_name),
                    amplicon_name,
                ),
            }
            plot_name = '1.'+amplicon_plot_name+'Editing_comparison'

            CRISPRessoPlot.plot_quantification_comparison_barchart(
                N_TOTAL_1,
                N_UNMODIFIED_1,
                N_MODIFIED_1,
                N_TOTAL_2,
                N_UNMODIFIED_2,
                N_MODIFIED_2,
                sample_1_name,
                sample_2_name,
                plot_titles,
                _jp(plot_name),
                save_png,
            )

            crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
            crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'Editing efficiency comparison'
            crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Figure 1: Comparison for amplicon ' + amplicon_name + '; Left: Percentage of modified and unmodified reads in each sample; Right: relative percentage of modified and unmodified reads'
            output_1 = os.path.join(args.crispresso_output_folder_1, run_info_1['running_info']['report_filename'])
            output_2 = os.path.join(args.crispresso_output_folder_2, run_info_2['running_info']['report_filename'])
            crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = []
            if os.path.isfile(output_1):
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name].append((sample_1_name +' output', os.path.relpath(output_1, OUTPUT_DIRECTORY)))
            if os.path.isfile(output_2):
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name].append((sample_2_name+ ' output', os.path.relpath(output_2, OUTPUT_DIRECTORY)))


            mod_file_1 = amplicon_info_1[amplicon_name]['modification_count_file']
            amp_seq_1, mod_freqs_1 = CRISPRessoShared.parse_count_file(mod_file_1)
            mod_file_2 = amplicon_info_2[amplicon_name]['modification_count_file']
            amp_seq_2, mod_freqs_2 = CRISPRessoShared.parse_count_file(mod_file_2)
            consensus_sequence = amp_seq_1
            if amp_seq_2 != consensus_sequence:
                raise DifferentAmpliconLengthException('Different amplicon lengths for the two amplicons.')

            quant_windows_1 = run_info_1['results']['refs'][amplicon_name]['include_idxs']
            quant_windows_2 = run_info_2['results']['refs'][amplicon_name]['include_idxs']
            quant_windows_are_equal = np.array_equal(quant_windows_1, quant_windows_2)

            for mod in ['Insertions', 'Deletions', 'Substitutions', 'All_modifications']:
                mod_name = mod
                if mod == "All_modifications":
                    mod_name = "Combined modifications (insertions, deletions and substitutions)"

                mod_counts_1 = np.array(mod_freqs_1[mod], dtype=float)
                tot_counts_1 = np.array(mod_freqs_1['Total'], dtype=float)
                unmod_counts_1 = tot_counts_1 - mod_counts_1

                mod_counts_2 = np.array(mod_freqs_2[mod], dtype=float)
                tot_counts_2 = np.array(mod_freqs_2['Total'], dtype=float)
                unmod_counts_2 = tot_counts_2 - mod_counts_2

                fisher_results = [stats.fisher_exact([[z[0], z[1]], [z[2], z[3]]]) if max(z) > 0 else [np.nan, 1.0] for z in zip(mod_counts_1, unmod_counts_1, mod_counts_2, unmod_counts_2)]
                oddsratios, pvalues = [a for a, b in fisher_results], [b for a, b in fisher_results]

                mod_df = []
                row = [sample_1_name+'_'+mod]
                row.extend(mod_counts_1)
                mod_df.append(row)

                row = [sample_1_name+'_total']
                row.extend(tot_counts_1)
                mod_df.append(row)

                row = [sample_2_name+'_'+mod]
                row.extend(mod_counts_2)
                mod_df.append(row)

                row = [sample_2_name+'_total']
                row.extend(tot_counts_2)
                mod_df.append(row)

                row = ['odds_ratios']
                row.extend(oddsratios)
                mod_df.append(row)

                row = ['pvalues']
                row.extend(pvalues)
                mod_df.append(row)

                m, pvals = len(pvalues), np.asarray(pvalues)
                qval_bonferroni = pvals * float(m)
                qval_bonferroni[np.where(qval_bonferroni > 1)] = 1
                row = ['qval_bonferroni']
                row.extend(qval_bonferroni)
                mod_df.append(row)

                sig_count = len(np.where(qval_bonferroni <= args.reported_qvalue_cutoff)[0])
                sig_counts[amplicon_name][mod] = sig_count

                sig_count_quant_window = np.nan
                if quant_windows_are_equal:
                    qvals_in_quant = np.take(qval_bonferroni, quant_windows_1)  # quant_windows_1 are the same as quant_windows_2
                    sig_count_quant_window = len(np.where(qvals_in_quant <= args.reported_qvalue_cutoff)[0])
                sig_counts_quant_window[amplicon_name][mod] = sig_count_quant_window

                colnames = ['Reference']
                colnames.extend(list(consensus_sequence))
                mod_df = pd.DataFrame(mod_df, columns=colnames)
#                mod_df = pd.concat([mod_df.iloc[:,0:2], mod_df.iloc[:,2:].apply(pd.to_numeric)],axis=1)
                #write to file
                mod_filename = _jp(amplicon_plot_name + mod + "_quantification.txt")
                mod_df.to_csv(mod_filename, sep='\t', index=None)

                #plot
                plot_title = get_plot_title_with_ref_name(
                    '%s: %s - %s' % (mod, sample_1_name, sample_2_name),
                    amplicon_name,
                )
                plot_name = '2.' + amplicon_plot_name + mod + '_quantification'
                CRISPRessoPlot.plot_quantification_positions(
                    mod_counts_1,
                    tot_counts_1,
                    mod_counts_2,
                    tot_counts_2,
                    len_amplicon,
                    pvalues,
                    len(consensus_sequence),
                    cut_points,
                    sgRNA_intervals,
                    plot_title,
                    _jp(plot_name),
                    save_png,
                )
                crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = mod_name +' locations'
                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = mod_name + ' location comparison for amplicon ' + amplicon_name + '; Top: percent difference; Bottom: p-value.'
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [(mod_name+' quantification', os.path.basename(mod_filename))]


            #create merged heatmaps for each cut site
            matching_allele_files = get_matching_allele_files(run_info_1, run_info_2)
            for allele_file_1, allele_file_2 in matching_allele_files:
                df1 = pd.read_csv(os.path.join(args.crispresso_output_folder_1, allele_file_1), sep="\t")
                df2 = pd.read_csv(os.path.join(args.crispresso_output_folder_2, allele_file_2), sep="\t")

                #find unmodified reference for comparison (if it exists)
                ref_seq_around_cut = ""
                if len(df1.loc[df1['Reference_Sequence'].str.contains('-')==False]) > 0:
                    ref_seq_around_cut = df1.loc[df1['Reference_Sequence'].str.contains('-')==False]['Reference_Sequence'].iloc[0]
                #otherwise figure out which sgRNA was used for this comparison
                elif len(df2.loc[df2['Reference_Sequence'].str.contains('-')==False]) > 0:
                    ref_seq_around_cut = df2.loc[df2['Reference_Sequence'].str.contains('-')==False]['Reference_Sequence'].iloc[0]
                else:
                    seq_len = df2[df2['Unedited']==True]['Reference_Sequence'].iloc[0]
                    for sgRNA_interval, cut_point in zip(sgRNA_intervals, cut_points):
                        sgRNA_seq = consensus_sequence[sgRNA_interval[0]:sgRNA_interval[1]]
                        if sgRNA_seq in allele_file_1:
                            this_sgRNA_seq = sgRNA_seq
                            this_cut_point = cut_point
                            ref_seq_around_cut=consensus_sequence[max(0, this_cut_point-args.offset_around_cut_to_plot+1):min(seq_len, cut_point+args.offset_around_cut_to_plot+1)]
                            break

                merged = pd.merge(df1, df2, on = ['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated'], suffixes=('_' + sample_1_name, '_'+sample_2_name), how='outer')
                quant_cols = ['#Reads_'+sample_1_name, '%Reads_'+sample_1_name, '#Reads_'+sample_2_name, '%Reads_'+sample_2_name]
                merged[quant_cols] = merged[quant_cols].fillna(0)
                lfc_error =0.1
                merged['each_LFC'] = np.log2(((merged['%Reads_'+sample_1_name]+lfc_error)/(merged['%Reads_'+sample_2_name]+lfc_error)).astype(float)).replace([np.inf, np.nan], 0)
                merged = merged.sort_values(['%Reads_'+sample_1_name, 'Reference_Sequence', 'n_deleted', 'n_inserted', 'n_mutated'], ascending=False)
                merged = merged.reset_index(drop=True).set_index('Aligned_Sequence')
                args.crispresso_output_folder_root = os.path.split(allele_file_1)[1].replace(".txt", "")
                allele_comparison_file = _jp(args.crispresso_output_folder_root+'.txt')
                merged.to_csv(allele_comparison_file, sep="\t", index=None)

                plot_name = '3.'+args.crispresso_output_folder_root+'_top'
                CRISPRessoPlot.plot_alleles_table_compare(ref_seq_around_cut, merged.sort_values(['each_LFC'], ascending=True), sample_1_name, sample_2_name, _jp(plot_name),
                            MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot, MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot, SAVE_ALSO_PNG=save_png)
                crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'Alleles enriched in ' + sample_1_name
                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Distribution comparison of alleles. Nucleotides are indicated by unique colors (A = green; C = red; G = yellow; T = purple). Substitutions are shown in bold font. Red rectangles highlight inserted sequences. Horizontal dashed lines indicate deleted sequences. The vertical dashed line indicates the predicted cleavage site. '+ \
                'The proportion and number of reads is shown for each sample on the right, with the values for ' + sample_1_name + ' followed by the values for ' + sample_2_name +'. Alleles are sorted for enrichment in ' + sample_1_name+'.'
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [('Allele comparison table', os.path.basename(allele_comparison_file))]

                plot_name = '3.'+args.crispresso_output_folder_root+'_bottom'
                CRISPRessoPlot.plot_alleles_table_compare(ref_seq_around_cut, merged.sort_values(['each_LFC'], ascending=False), sample_1_name, sample_2_name, _jp(plot_name),
                            MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot, MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot, SAVE_ALSO_PNG=save_png)
                crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'Alleles enriched in ' + sample_2_name
                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Distribution comparison of alleles. Nucleotides are indicated by unique colors (A = green; C = red; G = yellow; T = purple). Substitutions are shown in bold font. Red rectangles highlight inserted sequences. Horizontal dashed lines indicate deleted sequences. The vertical dashed line indicates the predicted cleavage site. '+ \
                'The proportion and number of reads is shown for each sample on the right, with the values for ' + sample_1_name + ' followed by the values for ' + sample_2_name +'. Alleles are sorted for enrichment in ' + sample_2_name+'.'
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [('Allele comparison table', os.path.basename(allele_comparison_file))]

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
                    line = "%s\t%s\t%s\t%s\n"%(amplicon_name, mod, val, val_quant_window)
                    fout.write(line)
        crispresso2_info['running_info']['sig_counts_report_location'] = sig_counts_filename

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Compare_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_compare_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT, logger)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        CRISPRessoShared.write_crispresso_info(crispresso2Compare_info_file, crispresso2_info)

        if args.zip_output:
            CRISPRessoShared.zip_results(OUTPUT_DIRECTORY)

        info('Analysis Complete!', {'percent_complete': 100})
        print(CRISPRessoShared.get_crispresso_footer())
        sys.exit(0)

    except Exception as e:
        debug_flag = False
        if 'args' in vars() and 'debug' in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)

        error('\n\nERROR: %s' % e)
        sys.exit(-1)
