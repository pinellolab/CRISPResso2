# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import os
import glob
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, wait
from functools import partial
import sys
import argparse
import numpy as np
import pandas as pd
import traceback
from datetime import datetime
from CRISPResso2 import CRISPRessoShared
from CRISPResso2.CRISPRessoReports import CRISPRessoReport
from CRISPResso2.CRISPRessoMultiProcessing import get_max_processes, run_plot

if CRISPRessoShared.is_C2Pro_installed():
    from CRISPRessoPro import __version__ as CRISPRessoProVersion
    C2PRO_INSTALLED = True
else:
    C2PRO_INSTALLED = False

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info

_ROOT = os.path.abspath(os.path.dirname(__file__))


def main():
    try:
        start_time =  datetime.now()
        start_time_string =  start_time.strftime('%Y-%m-%d %H:%M:%S')

        description = ['~~~CRISPRessoAggregate~~~', '-Aggregation of CRISPResso Run Data-']
        aggregate_string = r'''
___________________________________
|      __  __  _   _  __     ___ _ |
| /\  /__ /__ |_) |_ /__  /\  | |_ |
|/--\ \_| \_| | \ |_ \_| /--\ | |_ |
|__________________________________|
        '''
        info(CRISPRessoShared.get_crispresso_header(description, aggregate_string))

        parser = argparse.ArgumentParser(description="Aggregate CRISPResso2 Runs")
        parser.add_argument("-p", "--prefix", action='append', help="Prefix for CRISPResso folders to aggregate (may be specified multiple times)", default=[])
        parser.add_argument("-s", "--suffix", type=str, help="Suffix for CRISPResso folders to aggregate", default="")

        parser.add_argument("-n", "--name", type=str, help="Output name of the report", required=True)
        parser.add_argument('--min_reads_for_inclusion',  help='Minimum number of reads for a run to be included in the run summary', type=int, default=0)

        parser.add_argument('--place_report_in_output_folder',  help='If true, report will be written inside the CRISPResso output folder. By default, the report will be written one directory up from the report output.', action='store_true')
        parser.add_argument('--suppress_report',  help='Suppress output report', action='store_true')
        parser.add_argument('--suppress_plots',  help='Suppress output plots', action='store_true')
        parser.add_argument('--max_samples_per_summary_plot', type=int, help="Maximum number of samples on each page of the pdf report plots. If this number gets above ~150, they will be too big for matplotlib.", default=150)
        parser.add_argument('--n_processes', type=str, help='Specify the number of processes to use for analysis.\
        Please use with caution since increasing this parameter will significantly increase the memory required to run CRISPResso. Can be set to \'max\'.', default='1')

        parser.add_argument('--debug', help='Show debug messages', action='store_true')
        parser.add_argument('-v', '--verbosity', type=int, help='Verbosity level of output to the console (1-4), 4 is the most verbose', default=3)
        parser.add_argument('--halt_on_plot_fail', action="store_true", help="Halt execution if a plot fails to generate")

        # CRISPRessoPro params
        parser.add_argument('--use_matplotlib', action='store_true',
                        help='Use matplotlib for plotting instead of plotly/d3 when CRISPRessoPro is installed')

        args = parser.parse_args()

        if args.use_matplotlib or not CRISPRessoShared.is_C2Pro_installed():
            from CRISPResso2 import CRISPRessoPlot
        else:
            from CRISPRessoPro import plot as CRISPRessoPlot

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

        output_folder_name='CRISPRessoAggregate_on_%s' % args.name
        OUTPUT_DIRECTORY=os.path.abspath(output_folder_name)

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory

        try:
             info('Creating Folder %s' % OUTPUT_DIRECTORY)
             os.makedirs(OUTPUT_DIRECTORY)
        except:
             warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoAggregate_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        logger.addHandler(CRISPRessoShared.StatusHandler(os.path.join(OUTPUT_DIRECTORY, 'CRISPRessoAggregate_status.json')))

        with open(log_filename, 'w+') as outfile:
              outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Aggregate_info_file = os.path.join(
            OUTPUT_DIRECTORY, 'CRISPResso2Aggregate_info.json',
        )
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)
        crispresso2_info['running_info']['command_used'] = ' '.join(sys.argv)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        n_processes = 1
        if args.n_processes == 'max':
            n_processes = get_max_processes()
        else:
            n_processes = int(args.n_processes)

        if n_processes > 1:
            process_pool = ProcessPoolExecutor(n_processes)
            process_futures = {}
        else:
            process_pool = None
            process_futures = None

        plot = partial(
            run_plot,
            num_processes=n_processes,
            process_pool=process_pool,
            process_futures=process_futures,
            halt_on_plot_fail=args.halt_on_plot_fail,
        )

        #glob returns paths including the original prefix
        all_files = []
        for prefix in args.prefix:
            all_files.extend(glob.glob(prefix + '*'+args.suffix))
            if args.prefix != "":
                all_files.extend(glob.glob(prefix + '/*'+args.suffix)) #if a folder is given, add all subfolders

        seen_folders = {}
        crispresso2_folder_infos = {} #file_loc->crispresso_info; these are only CRISPResso runs -- this bit unrolls batch, pooled, and wgs runs
        successfully_imported_count = 0
        not_imported_count = 0
        for folder in all_files:
            if folder in seen_folders: #skip if we've seen this folder (glob could have added it twice)
                continue
            seen_folders[folder] = 1
            if os.path.isdir(folder) and str(folder).endswith(args.suffix):
                #first, try to import a plain CRISPResso2 run
                crispresso_info_file = os.path.join(folder, 'CRISPResso2_info.json')
                if os.path.exists(crispresso_info_file):
                    try:
                        run_data = CRISPRessoShared.load_crispresso_info(folder)
                        crispresso2_folder_infos[folder] = run_data
                        successfully_imported_count += 1
                    except Exception as e:
                        warn('Could not open CRISPResso2 info file in ' + folder)
                        not_imported_count += 1
                #second, check pooled
                pooled_info_file = os.path.join(folder, 'CRISPResso2Pooled_info.json')
                if os.path.exists(pooled_info_file):
                    pooled_data = CRISPRessoShared.load_crispresso_info(
                        folder, 'CRISPResso2Pooled_info.json',
                    )
                    if 'good_region_names' in pooled_data['results']:
                        run_names = pooled_data['results']['good_region_names']
                        for run_name in run_names:
                            run_folder_loc = os.path.join(folder, 'CRISPResso_on_%s'%run_name)
                            try:
                                run_data = CRISPRessoShared.load_crispresso_info(run_folder_loc)
                                crispresso2_folder_infos[run_folder_loc] = run_data
                                successfully_imported_count += 1
                            except Exception as e:
                                warn('Could not open CRISPResso2 info file in ' + run_folder_loc)
                                not_imported_count += 1
                    else:
                        warn('Could not process pooled folder ' + folder)
                        not_imported_count += 1
                #third, check batch
                batch_info_file = os.path.join(folder, 'CRISPResso2Batch_info.json')
                if os.path.exists(batch_info_file):
                    batch_data = CRISPRessoShared.load_crispresso_info(
                        folder, 'CRISPResso2Batch_info.json',
                    )
                    if 'completed_batch_arr' in batch_data['results']:
                        run_names = batch_data['results']['completed_batch_arr']
                        for run_name in run_names:
                            run_folder_loc = os.path.join(folder, 'CRISPResso_on_%s'%run_name)
                            try:
                                run_data = CRISPRessoShared.load_crispresso_info(run_folder_loc)
                                crispresso2_folder_infos[run_folder_loc] = run_data
                                successfully_imported_count += 1
                            except Exception as e:
                                warn('Could not open CRISPResso2 info file in ' + run_folder_loc)
                                not_imported_count += 1
                    else:
                        warn('Could not process batch folder ' + folder)
                        not_imported_count += 1
                #fourth, check WGS
                wgs_info_file = os.path.join(folder, 'CRISPResso2WGS_info.json')
                if os.path.exists(wgs_info_file):
                    wgs_data = CRISPRessoShared.load_crispresso_info(
                        folder, 'CRISPResso2WGS_info.json',
                    )
                    if 'good_region_folders' in wgs_data['results']:
                        run_names = wgs_data['results']['good_region_folders']
                        for run_name in run_names:
                            run_folder_loc = os.path.join(folder, 'CRISPResso_on_%s'%run_name)
                            try:
                                run_data = CRISPRessoShared.load_crispresso_info(run_folder_loc)
                                crispresso2_folder_infos[run_folder_loc] = run_data
                                successfully_imported_count += 1
                            except Exception:
                                warn('Could not open CRISPResso2 info file in ' + run_folder_loc)
                                not_imported_count += 1
                    else:
                        warn('Could not process WGS folder ' + folder)
                        not_imported_count += 1

        info('Read ' + str(successfully_imported_count) + ' folders (' + str(not_imported_count) + ' not imported)', {'percent_complete': 10})

        save_png = True
        if args.suppress_report:
            save_png = False

        if successfully_imported_count > 0:

            crispresso2_folders = list(sorted(crispresso2_folder_infos.keys()))
            crispresso2_folder_names = {}
            crispresso2_folder_htmls = {}#file_loc->html folder loc
            quilt_plots_to_show = {}  # name->{'href':path to report, 'img': png}
            for crispresso2_folder in crispresso2_folders:
                this_folder_name = CRISPRessoShared.slugify(crispresso2_folder)
                crispresso2_folder_names[crispresso2_folder] = this_folder_name
                this_sub_html_file = crispresso2_folder+".html"
                if crispresso2_folder_infos[crispresso2_folder]['running_info']['args'].place_report_in_output_folder:
                    this_sub_html_file = os.path.join(crispresso2_folder, crispresso2_folder_infos[crispresso2_folder]['running_info']['report_filename'])
                crispresso2_folder_htmls[crispresso2_folder] = os.path.abspath(this_sub_html_file)

                run_data = crispresso2_folder_infos[crispresso2_folder]
                for ref_name in run_data['results']['ref_names']:
                    if 'plot_2a_root' in run_data['results']['refs'][ref_name]:
                        plot_root = run_data['results']['refs'][ref_name]['plot_2a_root']
                        quilt_plots_to_show[this_folder_name+" "+ref_name] = {'href': os.path.abspath(this_sub_html_file),
                                'img': os.path.abspath(os.path.join(crispresso2_folder,plot_root+".png"))}

            all_amplicons = set()
            amplicon_names = {}  # sequence -> ref name (to check for amplicons with the same name but different sequences)
            amplicon_counts = {}
            amplicon_sources = {}
            for crispresso2_folder in crispresso2_folders:
                run_data = crispresso2_folder_infos[crispresso2_folder]
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
                        amplicon_sources[ref_seq] = []
                    amplicon_counts[ref_seq]+= 1
                    amplicon_sources[ref_seq].append(crispresso2_folder+'(' + ref_name + ')')

            # make sure amplicon names aren't super long
            for amplicon in all_amplicons:
                if len(amplicon_names[amplicon]) > 21:
                    amplicon_names[amplicon] = amplicon_names[amplicon][0:21]

            #make sure no duplicate amplicon names (same name for the different amplicons)
            seen_names = []
            for amplicon in all_amplicons:
                suffix_counter = 2
                orig_name = amplicon_names[amplicon]
                while amplicon_names[amplicon] in seen_names:
                    amplicon_names[amplicon] = orig_name+"_"+str(suffix_counter)
                    suffix_counter += 1
                seen_names.append(amplicon_names[amplicon])

            crispresso2_info['results']['ref_names'] = seen_names
            crispresso2_info['results']['refs'] = {}
            crispresso2_info['results']['general_plots']['summary_plot_names'] = []
            crispresso2_info['results']['general_plots']['summary_plot_titles'] = {}
            crispresso2_info['results']['general_plots']['summary_plot_labels'] = {}
            crispresso2_info['results']['general_plots']['summary_plot_datas'] = {}

            with open(_jp('CRISPRessoAggregate_amplicon_information.txt'), 'w') as outfile:
                outfile.write("\t".join(['Amplicon Name', 'Number of sources', 'Amplicon sources', 'Amplicon sequence']) + "\n")
                for amplicon in all_amplicons:
                    outfile.write("\t".join([amplicon_names[amplicon], str(amplicon_counts[amplicon]), ';'.join(amplicon_sources[amplicon]), amplicon]) + "\n")

            window_nuc_pct_quilt_plot_names = []
            nuc_pct_quilt_plot_names = []
            window_nuc_conv_plot_names = []
            nuc_conv_plot_names = []

            percent_complete_start, percent_complete_end = 11, 90
            percent_complete_step = (percent_complete_end - percent_complete_start) / len(all_amplicons)
            #report for amplicons that appear multiple times
            for amplicon_index, amplicon_seq in enumerate(all_amplicons):
                amplicon_name = amplicon_names[amplicon_seq]
                crispresso2_info['results']['refs'][amplicon_name] = {}
                #only perform comparison if amplicon seen in more than one sample
                if amplicon_counts[amplicon_seq] < 2:
                    continue

                percent_complete = percent_complete_start + (amplicon_index * percent_complete_step)
                info('Reporting summary for amplicon: "' + amplicon_name + '"', {'percent_complete': percent_complete})

                consensus_sequence = ""
                nucleotide_frequency_summary = []
                nucleotide_percentage_summary = []
                modification_frequency_summary = []
                modification_percentage_summary = []

                amp_found_count = 0 #how many folders had information for this amplicon
                consensus_guides = []
                consensus_include_idxs = []
                consensus_sgRNA_plot_idxs = []
                consensus_sgRNA_intervals = []
                guides_all_same = True
                runs_with_this_amplicon = []
                for crispresso2_folder in crispresso2_folders:
                    run_data = crispresso2_folder_infos[crispresso2_folder]
                    run_has_amplicon = False
                    run_amplicon_name = ''
                    for ref_name in run_data['results']['ref_names']:
                        if amplicon_seq == run_data['results']['refs'][ref_name]['sequence']:
                            run_has_amplicon = True
                            run_amplicon_name = ref_name
                    if not run_has_amplicon:
                        continue
                    runs_with_this_amplicon.append(crispresso2_folder)

                    if consensus_guides == []:
                        consensus_guides = run_data['results']['refs'][run_amplicon_name]['sgRNA_sequences']
                        consensus_include_idxs = run_data['results']['refs'][run_amplicon_name]['include_idxs']
                        consensus_sgRNA_intervals = run_data['results']['refs'][run_amplicon_name]['sgRNA_intervals']
                        consensus_sgRNA_plot_idxs = run_data['results']['refs'][run_amplicon_name]['sgRNA_plot_idxs']

                    if run_data['results']['refs'][run_amplicon_name]['sgRNA_sequences'] != consensus_guides:
                        guides_all_same = False
                    if set(run_data['results']['refs'][run_amplicon_name]['include_idxs']) != set(consensus_include_idxs):
                        guides_all_same = False

                    if 'nuc_freq_filename' not in run_data['results']['refs'][run_amplicon_name]:
                        info("Skipping the amplicon '%s' in folder '%s'. Cannot find nucleotide information."%(run_amplicon_name, crispresso2_folder))
                        continue

                    nucleotide_frequency_file = os.path.join(crispresso2_folder, run_data['results']['refs'][run_amplicon_name]['nuc_freq_filename'])
                    ampSeq_nf, nuc_freqs = CRISPRessoShared.parse_count_file(nucleotide_frequency_file)

                    nucleotide_pct_file = os.path.join(crispresso2_folder, run_data['results']['refs'][run_amplicon_name]['nuc_pct_filename'])
                    ampSeq_np, nuc_pcts = CRISPRessoShared.parse_count_file(nucleotide_pct_file)

                    count_file = os.path.join(crispresso2_folder, run_data['results']['refs'][run_amplicon_name]['mod_count_filename'])
                    ampSeq_cf, mod_freqs = CRISPRessoShared.parse_count_file(count_file)

                    if ampSeq_nf is None or ampSeq_np is None or ampSeq_cf is None:
                        info("Skipping the amplicon '%s' in folder '%s'. Could not parse run output."%(run_amplicon_name, crispresso2_folder))
                        info("Nucleotide frequency amplicon: '%s', Nucleotide percentage amplicon: '%s', Count vectors amplicon: '%s'"%(ampSeq_nf, ampSeq_np, ampSeq_cf))
                        continue
                    if ampSeq_nf != ampSeq_np or ampSeq_np != ampSeq_cf:
                        warn("Skipping the amplicon '%s' in folder '%s'. Parsed amplicon sequences do not match\nnf:%s\nnp:%s\ncf:%s\nrf:%s"%(run_amplicon_name, crispresso2_folder, ampSeq_nf, ampSeq_np, ampSeq_cf, amplicon_seq))
                        continue
                    if consensus_sequence == "":
                        consensus_sequence = ampSeq_nf
                    if ampSeq_nf != consensus_sequence:
                        info("Skipping the amplicon '%s' in folder '%s'. Amplicon sequences do not match."%(run_amplicon_name, crispresso2_folder))
                        continue
                    if 'Total' not in mod_freqs:
                        info("Skipping the amplicon '%s' in folder '%s'. Processing did not complete."%(run_amplicon_name, crispresso2_folder))
                        continue
                    if mod_freqs['Total'][0] == 0 or mod_freqs['Total'][0] == "0":
                        info("Skipping the amplicon '%s' in folder '%s'. Got no reads for amplicon."%(run_amplicon_name, crispresso2_folder))
                        continue
                    this_amp_total_reads = run_data['results']['alignment_stats']['counts_total'][run_amplicon_name]
                    if this_amp_total_reads < args.min_reads_for_inclusion:
                        info("Skipping the amplicon '%s' in folder '%s'. Got %s reads (min_reads_for_inclusion is %d)."%(run_amplicon_name, crispresso2_folder, str(this_amp_total_reads), args.min_reads_for_inclusion))
                        continue

                    mod_pcts = {}
                    for key in mod_freqs:
                        mod_pcts[key] = np.array(
                            mod_freqs[key],
                        ).astype(float) / float(this_amp_total_reads)

                    amp_found_count += 1

                    run_name = crispresso2_folder_names[crispresso2_folder]

                    for nuc in ['A', 'T', 'C', 'G', 'N', '-']:
                        row = [run_name, nuc]
                        row.extend(nuc_freqs[nuc])
                        nucleotide_frequency_summary.append(row)

                        pct_row = [run_name, nuc]
                        pct_row.extend(nuc_pcts[nuc])
                        nucleotide_percentage_summary.append(pct_row)

                    for mod in ['Insertions', 'Insertions_Left', 'Deletions', 'Substitutions', 'All_modifications']:
                        row = [run_name, mod]
                        row.extend(mod_freqs[mod])
                        modification_frequency_summary.append(row)

                        pct_row = [run_name, mod]
                        pct_row.extend(mod_pcts[mod])
                        modification_percentage_summary.append(pct_row)

                if amp_found_count == 0:
                    info("Couldn't find any data for amplicon '%s'. Not compiling results."%amplicon_name)
                else:
                    amplicon_plot_name = amplicon_name+"."
                    if len(amplicon_names) == 1 and amplicon_name == "Reference":
                        amplicon_plot_name = ""

                    colnames = ['Folder', 'Nucleotide']
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

                    colnames = ['Folder', 'Modification']
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

                    crispresso2_info['results']['refs'][amplicon_name]['nucleotide_frequency_summary_filename'] = os.path.basename(nucleotide_frequency_summary_filename)
                    crispresso2_info['results']['refs'][amplicon_name]['nucleotide_percentage_summary_filename'] = os.path.basename(nucleotide_percentage_summary_filename)

                    crispresso2_info['results']['refs'][amplicon_name]['modification_frequency_summary_filename'] = os.path.basename(modification_frequency_summary_filename)
                    crispresso2_info['results']['refs'][amplicon_name]['modification_percentage_summary_filename'] = os.path.basename(modification_percentage_summary_filename)

                    this_number_samples = len(pd.unique(nucleotide_percentage_summary_df['Folder']))

                    #if guides are all the same, merge substitutions and perform base editor comparison at guide quantification window
                    if guides_all_same and consensus_guides != []:
                        info("All guides are equal. Performing comparison of runs for amplicon '%s'"% amplicon_name)
                        include_idxs = consensus_include_idxs #include indexes are the same for all guides
                        for idx, sgRNA in enumerate(consensus_guides):
                            sgRNA_plot_idxs = consensus_sgRNA_plot_idxs[idx]
                            plot_idxs_flat = [0, 1] # guide, nucleotide
                            plot_idxs_flat.extend([plot_idx + 2 for plot_idx in sgRNA_plot_idxs])

                            sub_nucleotide_frequency_summary_df = nucleotide_frequency_summary_df.iloc[:, plot_idxs_flat]
                            sub_nucleotide_frequency_summary_df = pd.concat([sub_nucleotide_frequency_summary_df.iloc[:, 0:2],
                                                                        sub_nucleotide_frequency_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                            sub_nucleotide_frequency_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_frequency_summary_around_sgRNA_'+sgRNA+'.txt')
                            sub_nucleotide_frequency_summary_df.to_csv(sub_nucleotide_frequency_summary_filename, sep='\t', index=None)

                            sub_nucleotide_percentage_summary_df = nucleotide_percentage_summary_df.iloc[:, plot_idxs_flat]
                            sub_nucleotide_percentage_summary_df = pd.concat([sub_nucleotide_percentage_summary_df.iloc[:, 0:2],
                                                                        sub_nucleotide_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                            sub_nucleotide_percentage_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_percentage_summary_around_sgRNA_'+sgRNA+'.txt')
                            sub_nucleotide_percentage_summary_df.to_csv(sub_nucleotide_percentage_summary_filename, sep='\t', index=None)

                            sub_modification_percentage_summary_df = modification_percentage_summary_df.iloc[:, plot_idxs_flat]
                            sub_modification_percentage_summary_df = pd.concat([sub_modification_percentage_summary_df.iloc[:, 0:2],
                                                                        sub_modification_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                            sub_modification_percentage_summary_filename = _jp(amplicon_plot_name + 'Modification_percentage_summary_around_sgRNA_'+sgRNA+'.txt')
                            sub_modification_percentage_summary_df.to_csv(sub_modification_percentage_summary_filename, sep='\t', index=None)

                            if not args.suppress_plots and this_number_samples < args.max_samples_per_summary_plot:
                                # plot for each guide
                                # show all sgRNA's on the plot
                                sub_sgRNA_intervals = []
                                for sgRNA_interval in consensus_sgRNA_intervals:
                                    newstart = None
                                    newend = None
                                    for idx, i in enumerate(sgRNA_plot_idxs):
                                        if i <= sgRNA_interval[0]:
                                            newstart = idx
                                        if newend is None and i >= sgRNA_interval[1]:
                                            newend = idx

                                    #if guide doesn't overlap with plot idxs
                                    if newend == 0 or newstart == len(sgRNA_plot_idxs):
                                        continue
                                    #otherwise, correct partial overlaps
                                    elif newstart == None and newend == None:
                                        newstart = 0
                                        newend = len(include_idxs) -1
                                    elif newstart == None:
                                        newstart = 0
                                    elif newend == None:
                                        newend = len(include_idxs) -1
                                    #and add it to the list
                                    sub_sgRNA_intervals.append((newstart, newend))

                                this_window_nuc_pct_quilt_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt_around_sgRNA_'+sgRNA)
                                nucleotide_quilt_input = {
                                    'nuc_pct_df': sub_nucleotide_percentage_summary_df,
                                    'mod_pct_df': sub_modification_percentage_summary_df,
                                    'fig_filename_root': this_window_nuc_pct_quilt_plot_name,
                                    'save_also_png': save_png,
                                    'sgRNA_intervals': sub_sgRNA_intervals,
                                    'sgRNA_sequences': consensus_guides,
                                    'quantification_window_idxs': include_idxs,
                                    'group_column': 'Folder',
                                    'custom_colors': None,
                                }
                                plot(
                                    CRISPRessoPlot.plot_nucleotide_quilt,
                                    nucleotide_quilt_input,
                                )

                                plot_name = os.path.basename(this_window_nuc_pct_quilt_plot_name)
                                window_nuc_pct_quilt_plot_names.append(plot_name)
                                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'sgRNA: ' + sgRNA + ' Amplicon: ' + amplicon_name
                                if len(consensus_guides) == 1:
                                    crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = ''
                                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Composition of each base around the guide ' + sgRNA + ' for the amplicon ' + amplicon_name
                                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [(amplicon_name + ' nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), (amplicon_name + ' modification frequencies', os.path.basename(modification_frequency_summary_filename))]
                        # done with per-sgRNA plots

                        if not args.suppress_plots: # and this_number_samples < 500: # plot the whole region
                            this_plot_suffix = "" # in case we have a lot of regions, split them up and add a suffix here
                            this_plot_suffix_int = 1
                            nrow_per_sample_nucs = nucleotide_percentage_summary_df.shape[0] / this_number_samples # calculate number of rows per sample for subsetting the tables
                            nrow_per_sample_mods = modification_percentage_summary_df.shape[0] / this_number_samples
                            for sample_start_ind in range(0, this_number_samples, args.max_samples_per_summary_plot):
                                sample_end_ind = min(sample_start_ind + args.max_samples_per_summary_plot, this_number_samples)
                                this_nuc_pct_quilt_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt' + this_plot_suffix)
                                this_nuc_start_ind = int(sample_start_ind*nrow_per_sample_nucs)
                                this_nuc_end_ind = int((sample_end_ind+1)*nrow_per_sample_nucs - 1)
                                this_mod_start_ind = int(sample_start_ind*nrow_per_sample_mods)
                                this_mod_end_ind = int((sample_end_ind+1)*nrow_per_sample_mods - 1)
                                nucleotide_quilt_input = {
                                    'nuc_pct_df': nucleotide_percentage_summary_df.iloc[this_nuc_start_ind:this_nuc_end_ind, :],
                                    'mod_pct_df': modification_percentage_summary_df.iloc[this_mod_start_ind:this_mod_end_ind, :],
                                    'fig_filename_root': this_nuc_pct_quilt_plot_name,
                                    'save_also_png': save_png,
                                    'sgRNA_intervals': consensus_sgRNA_intervals,
                                    'sgRNA_sequences': consensus_guides,
                                    'quantification_window_idxs': include_idxs,
                                    'group_column': 'Folder',
                                    'custom_colors': None,
                                }
                                plot(
                                    CRISPRessoPlot.plot_nucleotide_quilt,
                                    nucleotide_quilt_input,
                                )

                                plot_name = os.path.basename(this_nuc_pct_quilt_plot_name)
                                nuc_pct_quilt_plot_names.append(plot_name)
                                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'Amplicon: ' + amplicon_name + this_plot_suffix
                                if len(amplicon_names) == 1:
                                    crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = ''
                                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Composition of each base for the amplicon ' + amplicon_name
                                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [(amplicon_name + ' nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), (amplicon_name + ' modification frequencies', os.path.basename(modification_frequency_summary_filename))]

                                this_plot_suffix_int += 1
                                this_plot_suffix = "_" + str(this_plot_suffix_int)

                    else: # guides are not the same
                        if not args.suppress_plots: # and this_number_samples < 150:
                            this_plot_suffix = "" # in case we have a lot of regions, split them up and add a suffix here
                            this_plot_suffix_int = 1
                            nrow_per_sample_nucs = nucleotide_percentage_summary_df.shape[0] / this_number_samples # calculate number of rows per sample for subsetting the tables
                            nrow_per_sample_mods = modification_percentage_summary_df.shape[0] / this_number_samples
                            for sample_start_ind in range(0,this_number_samples,args.max_samples_per_summary_plot):
                                sample_end_ind = min(sample_start_ind + args.max_samples_per_summary_plot, this_number_samples)
                                this_nuc_pct_quilt_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt' + this_plot_suffix)
                                this_nuc_start_ind = int(sample_start_ind*nrow_per_sample_nucs)
                                this_nuc_end_ind = int((sample_end_ind+1)*nrow_per_sample_nucs - 1)
                                this_mod_start_ind = int(sample_start_ind*nrow_per_sample_mods)
                                this_mod_end_ind = int((sample_end_ind+1)*nrow_per_sample_mods - 1)

                                nucleotide_quilt_input = {
                                    'nuc_pct_df': nucleotide_percentage_summary_df.iloc[this_nuc_start_ind:this_nuc_end_ind, :],
                                    'mod_pct_df': modification_percentage_summary_df.iloc[this_mod_start_ind:this_mod_end_ind, :],
                                    'fig_filename_root': this_nuc_pct_quilt_plot_name,
                                    'save_also_png': save_png,
                                    'sgRNA_intervals': consensus_sgRNA_intervals,
                                    'sgRNA_sequences': consensus_guides,
                                    'quantification_window_idxs': consensus_include_idxs,
                                    'group_column': 'Folder',
                                    'custom_colors': None,
                                }
                                plot(
                                    CRISPRessoPlot.plot_nucleotide_quilt,
                                    nucleotide_quilt_input,
                                )

                                plot_name = os.path.basename(this_nuc_pct_quilt_plot_name)
                                nuc_pct_quilt_plot_names.append(plot_name)
                                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'Amplicon: ' + amplicon_name + this_plot_suffix
                                if len(amplicon_names) == 1:
                                    crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = ''
                                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Composition of each base for the amplicon ' + amplicon_name
                                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [(amplicon_name + ' nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), (amplicon_name + ' modification frequencies', os.path.basename(modification_frequency_summary_filename))]

                                this_plot_suffix_int += 1
                                this_plot_suffix = "_" + str(this_plot_suffix_int)

                    if C2PRO_INSTALLED and not args.use_matplotlib and not args.suppress_plots:
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
                        if guides_all_same:
                            sgRNA_intervals = [consensus_sgRNA_intervals] * modification_frequency_summary_df.shape[0]
                        else:
                            sgRNA_intervals = [consensus_sgRNA_intervals]
                        for modification_type in ['Insertions', 'Deletions', 'Substitutions']:
                            modification_df = modification_frequency_summary_df[
                                modification_frequency_summary_df['Modification'] == modification_type
                            ]
                            modification_df.index = [
                                '{0} ({1})'.format(folder, folder_index)
                                for folder_index, folder in enumerate(
                                    modification_df['Folder'], 1,
                                )
                            ]
                            modification_df = modification_df.drop(
                                ['Modification', 'Folder'], axis=1,
                            )
                            modification_df.columns = [
                                '{0} ({1})'.format(column, position)
                                for position, column in
                                enumerate(modification_df.columns, 1)
                            ]
                            plot_name = 'CRISPRessoAggregate_percentage_of_{0}_across_alleles_{1}_heatmap'.format(modification_type.lower(), amplicon_name)
                            plot_path = '{0}.html'.format(_jp(plot_name))

                            heatmap_div_id = '{0}-allele-modification-heatmap-{1}'.format(amplicon_name.lower(), modification_type.lower())
                            allele_modification_heatmap_input = {
                                'sample_values': modification_df,
                                'sample_sgRNA_intervals': sgRNA_intervals,
                                'plot_path': plot_path,
                                'title': modification_type,
                                'div_id': heatmap_div_id,
                                'amplicon_name': amplicon_name,
                            }
                            plot(
                                CRISPRessoPlot.plot_allele_modification_heatmap,
                                allele_modification_heatmap_input,
                            )

                            crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_names'].append(plot_name)
                            crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_paths'][plot_name] = plot_path
                            crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_titles'][plot_name] = 'CRISPRessoAggregate {0} Across Samples for {1}'.format(
                                modification_type,
                                amplicon_name,
                            )
                            crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_labels'][plot_name] = 'Each row is a sample and each column is a position in the amplicon sequence. Each cell shows the percentage of {0} for the sample at that position relative to the amplicon. Guides for each sample are identified by a black rectangle.'.format(modification_type.lower())
                            crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_datas'][plot_name] = [
                                (
                                    'CRISPRessoAggregate Modification Frequency Summary',
                                    os.path.basename(
                                        modification_frequency_summary_filename,
                                    ),
                                ),
                            ]
                            crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_divs'][plot_name] = heatmap_div_id

                            plot_name = 'CRISPRessoAggregate_percentage_of_{0}_across_alleles_{1}_line'.format(modification_type.lower(), amplicon_name)
                            plot_path = '{0}.html'.format(_jp(plot_name))

                            line_div_id = '{0}-allele-modification-line-{1}'.format(amplicon_name.lower(), modification_type.lower())
                            allele_modification_line_input = {
                                'sample_values': modification_df,
                                'sample_sgRNA_intervals': sgRNA_intervals,
                                'plot_path': plot_path,
                                'title': modification_type,
                                'div_id': line_div_id,
                                'amplicon_name': amplicon_name,
                            }
                            plot(
                                CRISPRessoPlot.plot_allele_modification_line,
                                allele_modification_line_input,
                            )
                            crispresso2_info['results']['general_plots']['allele_modification_line_plot_names'].append(plot_name)
                            crispresso2_info['results']['general_plots']['allele_modification_line_plot_paths'][plot_name] = plot_path
                            crispresso2_info['results']['general_plots']['allele_modification_line_plot_titles'][plot_name] = 'CRISPRessoAggregate {0} Across Samples for {1}'.format(
                                modification_type,
                                amplicon_name,
                            )
                            crispresso2_info['results']['general_plots']['allele_modification_line_plot_labels'][plot_name] = 'Each line is a sample that indicates the percentage of {0} for the sample at that position relative to the amplicon. Guides are shown by a grey rectangle.'.format(modification_type.lower())
                            crispresso2_info['results']['general_plots']['allele_modification_line_plot_datas'][plot_name] = [
                                (
                                    'CRISPRessoAggregate Modification Frequency Summary',
                                    os.path.basename(
                                        modification_frequency_summary_filename,
                                    ),
                                ),
                            ]
                            crispresso2_info['results']['general_plots']['allele_modification_line_plot_divs'][plot_name] = line_div_id

            crispresso2_info['results']['general_plots']['window_nuc_pct_quilt_plot_names'] = window_nuc_pct_quilt_plot_names
            crispresso2_info['results']['general_plots']['nuc_pct_quilt_plot_names'] = nuc_pct_quilt_plot_names
            crispresso2_info['results']['general_plots']['window_nuc_conv_plot_names'] = window_nuc_conv_plot_names
            crispresso2_info['results']['general_plots']['nuc_conv_plot_names'] = nuc_conv_plot_names

            quantification_summary=[]
            #summarize amplicon modifications
            debug('Summarizing amplicon modifications...', {'percent_complete': 92})
            samples_quantification_summary_by_amplicon_filename = _jp('CRISPRessoAggregate_quantification_of_editing_frequency_by_amplicon.txt') #this file has separate lines for each amplicon in each run
            with open(samples_quantification_summary_by_amplicon_filename, 'w') as outfile:
                wrote_header = False
                for crispresso2_folder in crispresso2_folders:
                    run_data = crispresso2_folder_infos[crispresso2_folder]
                    run_name = crispresso2_folder_names[crispresso2_folder]
                    amplicon_modification_file=os.path.join(crispresso2_folder, run_data['running_info']['quant_of_editing_freq_filename'])
                    with open(amplicon_modification_file, 'r') as infile:
                        file_head = infile.readline()
                        if not wrote_header:
                            outfile.write('Folder\t' + file_head)
                            wrote_header = True
                        for line in infile:
                            outfile.write(crispresso2_folder + "\t" + line)

                    n_tot = run_data['running_info']['alignment_stats']['N_TOT_READS']
                    n_aligned = 0
                    n_unmod = 0
                    n_mod = 0
                    n_discarded = 0

                    n_insertion = 0
                    n_deletion = 0
                    n_substitution = 0
                    n_only_insertion = 0
                    n_only_deletion = 0
                    n_only_substitution = 0
                    n_insertion_and_deletion = 0
                    n_insertion_and_substitution = 0
                    n_deletion_and_substitution = 0
                    n_insertion_and_deletion_and_substitution = 0

                    for ref_name in run_data['results']['ref_names']: #multiple alleles could be provided
                        n_aligned += run_data['results']['alignment_stats']['counts_total'][ref_name]
                        n_unmod += run_data['results']['alignment_stats']['counts_unmodified'][ref_name]
                        n_mod += run_data['results']['alignment_stats']['counts_modified'][ref_name]
                        n_discarded += run_data['results']['alignment_stats']['counts_discarded'][ref_name]

                        n_insertion += run_data['results']['alignment_stats']['counts_insertion'][ref_name]
                        n_deletion += run_data['results']['alignment_stats']['counts_deletion'][ref_name]
                        n_substitution += run_data['results']['alignment_stats']['counts_substitution'][ref_name]
                        n_only_insertion += run_data['results']['alignment_stats']['counts_only_insertion'][ref_name]
                        n_only_deletion += run_data['results']['alignment_stats']['counts_only_deletion'][ref_name]
                        n_only_substitution += run_data['results']['alignment_stats']['counts_only_substitution'][ref_name]
                        n_insertion_and_deletion += run_data['results']['alignment_stats']['counts_insertion_and_deletion'][ref_name]
                        n_insertion_and_substitution += run_data['results']['alignment_stats']['counts_insertion_and_substitution'][ref_name]
                        n_deletion_and_substitution += run_data['results']['alignment_stats']['counts_deletion_and_substitution'][ref_name]
                        n_insertion_and_deletion_and_substitution += run_data['results']['alignment_stats']['counts_insertion_and_deletion_and_substitution'][ref_name]

                    unmod_pct = np.nan
                    mod_pct = np.nan
                    if n_aligned > 0:
                        unmod_pct = 100*n_unmod/float(n_aligned)
                        mod_pct = 100*n_mod/float(n_aligned)


                    vals = [run_name]
                    vals.extend([round(unmod_pct, 8), round(mod_pct, 8), n_aligned, n_tot, n_unmod, n_mod, n_discarded, n_insertion, n_deletion, n_substitution, n_only_insertion, n_only_deletion, n_only_substitution, n_insertion_and_deletion, n_insertion_and_substitution, n_deletion_and_substitution, n_insertion_and_deletion_and_substitution])
                    quantification_summary.append(vals)

            header = 'Name\tUnmodified%\tModified%\tReads_total\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions'
            header_els = header.split("\t")
            df_summary_quantification=pd.DataFrame(quantification_summary, columns=header_els).sort_values(by=['Name'])
            samples_quantification_summary_filename = _jp('CRISPRessoAggregate_quantification_of_editing_frequency.txt') #this file has one line for each run (sum of all amplicons)
            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename, sep='\t', index=None)
            crispresso2_info['results']['alignment_stats']['samples_quantification_summary_filename'] = os.path.basename(samples_quantification_summary_filename)
            crispresso2_info['results']['alignment_stats']['samples_quantification_summary_by_amplicon_filename'] = os.path.basename(samples_quantification_summary_by_amplicon_filename)
            df_summary_quantification.set_index('Name')

            if not args.suppress_plots:
                plot_root = _jp("CRISPRessoAggregate_reads_summary")
                debug('Plotting reads summary...', {'percent_complete': 94})
                reads_total_input = {
                    'fig_filename_root': plot_root,
                    'df_summary_quantification': df_summary_quantification,
                    'save_png': save_png,
                    'cutoff': args.min_reads_for_inclusion,
                }
                plot(CRISPRessoPlot.plot_reads_total, reads_total_input)

                plot_name = os.path.basename(plot_root)
                crispresso2_info['results']['general_plots']['summary_plot_root'] = plot_name
                crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'CRISPRessoAggregate Mapping Statistics Summary'
                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads in each sample. The vertical line shows the cutoff for analysis, set using the --min_reads_for_inclusion parameter.'
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [('CRISPRessoAggregate summary', os.path.basename(samples_quantification_summary_filename)), ('CRISPRessoAggregate summary by amplicon', os.path.basename(samples_quantification_summary_by_amplicon_filename))]

                plot_root = _jp("CRISPRessoAggregate_quantification_of_editing_frequency")

                unmod_mod_pcts_input = {
                    'fig_filename_root': plot_root,
                    'df_summary_quantification': df_summary_quantification,
                    'save_png': save_png,
                    'cutoff': args.min_reads_for_inclusion,
                }
                plot(CRISPRessoPlot.plot_unmod_mod_pcts, unmod_mod_pcts_input)

                plot_name = os.path.basename(plot_root)
                crispresso2_info['results']['general_plots']['summary_plot_root'] = plot_name
                crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
                crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'CRISPRessoAggregate Modification Summary'
                crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads aligned to each amplicon, divided into the reads that are modified and unmodified. The vertical line shows the cutoff for analysis, set using the --min_reads_for_inclusion parameter.'
                crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [('CRISPRessoAggregate summary', os.path.basename(samples_quantification_summary_filename)), ('CRISPRessoAggregate summary by amplicon', os.path.basename(samples_quantification_summary_by_amplicon_filename))]

            #summarize alignment
            debug('Summarizing alignment...', {'percent_complete': 96})
            with open(_jp('CRISPRessoAggregate_mapping_statistics.txt'), 'w') as outfile:
                wrote_header = False
                for crispresso2_folder in crispresso2_folders:
                    run_data = crispresso2_folder_infos[crispresso2_folder]
                    run_name = crispresso2_folder_names[crispresso2_folder]
                    mapping_file=os.path.join(crispresso2_folder, run_data['running_info']['mapping_stats_filename'])
                    with open(mapping_file, 'r') as infile:
                        file_head = infile.readline()
                        if not wrote_header:
                            outfile.write('Folder\t' + file_head)
                            wrote_header = True
                        for line in infile:
                            outfile.write(crispresso2_folder + "\t" + line)

            if not args.suppress_report:
                report_filename = OUTPUT_DIRECTORY+'.html'
                if (args.place_report_in_output_folder):
                    report_filename = _jp("CRISPResso2Aggregate_report.html")
                CRISPRessoReport.make_aggregate_report(
                    crispresso2_info,
                    args.name,
                    report_filename,
                    OUTPUT_DIRECTORY,
                    _ROOT,
                    crispresso2_folders,
                    crispresso2_folder_htmls,
                    logger,
                    compact_plots_to_show=quilt_plots_to_show,
                )
                crispresso2_info['running_info']['report_location'] = report_filename
                crispresso2_info['running_info']['report_filename'] = os.path.basename(report_filename)
        else: #no files successfully imported
            files_in_curr_dir = os.listdir('.')
            if len(files_in_curr_dir) > 15:
                files_in_curr_dir = files_in_curr_dir[0:15]
                files_in_curr_dir.append("(Complete listing truncated)")
            info('No CRISPResso runs could be imported.\nFiles in current directory:\n\t' + "\n\t".join(files_in_curr_dir))

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
            crispresso2Aggregate_info_file, crispresso2_info,
        )

        if n_processes > 1:
            wait(process_futures)
            if args.debug:
                debug('Plot pool results:')
                for future in process_futures:
                    debug('future: ' + str(future))
            for future in process_futures:
                try:
                    future.result()
                except Exception as e:
                    logger.warning('Error in plot pool: %s' % e)
                    logger.debug(traceback.format_exc())
            process_pool.shutdown()

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

if __name__ == '__main__':
    main()
