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
from datetime import datetime
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoPlot
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoReport

import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info

_ROOT = os.path.abspath(os.path.dirname(__file__))

####Support functions###
def propagate_options(cmd, options, params, paramInd):
####
# cmd - the command to run
# options - list of options to propagate e.g. crispresso options
# params - df from excel parser e.g. params['amplicon_name'] = ['name1','name2','name3'] where each item corresponds to a different run in the batch
# paramInd - index in dict - this is the run number in the batch

    for option in options :
        if option:
            if option in params:
                val = params.loc[paramInd, option]
                if val is None:
                    pass
                elif str(val) == "True":
                    cmd+=' --%s' % option
                elif str(val) =="False":
                    pass
                elif isinstance(val, str):
                    if val != "":
                        if " " in val or "-" in val:
                            cmd+=' --%s "%s"' % (option, str(val)) # quotes for options with spaces
                        else:
                            cmd+=' --%s %s' % (option, str(val))
                elif isinstance(val, bool):
                    if val:
                        cmd+=' --%s' % option
                else:
                    cmd+=' --%s %s' % (option, str(val))
#    print("cmd is " + str(cmd))
    return cmd

def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoBatch!' % library_name)
                sys.exit(1)



pd=check_library('pandas')
np=check_library('numpy')

def main():
    try:
        start_time =  datetime.now()
        start_time_string =  start_time.strftime('%Y-%m-%d %H:%M:%S')

        description = ['~~~CRISPRessoBatch~~~', '-Analysis of CRISPR/Cas9 outcomes from batch deep sequencing data-']
        batch_string = r'''
 _________________
| __    ___ __    |
||__) /\ | /  |__||
||__)/--\| \__|  ||
|_________________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description, batch_string))

        parser = CRISPRessoShared.getCRISPRessoArgParser(parserTitle = 'CRISPRessoBatch Parameters')

        #batch specific params
        parser.add_argument('-bs', '--batch_settings', type=str, help='Settings file for batch. Must be tab-separated text file. The header row contains CRISPResso parameters (e.g., fastq_r1, fastq_r2, amplicon_seq, and other optional parameters). Each following row sets parameters for an additional batch.', required=True)
        parser.add_argument('--skip_failed',  help='Continue with batch analysis even if one sample fails', action='store_true')
        parser.add_argument('--min_reads_for_inclusion',  help='Minimum number of reads for a batch to be included in the batch summary', type=int, default=0)
        parser.add_argument('-bo', '--batch_output_folder',  help='Directory where batch analysis output will be stored')
        parser.add_argument('--crispresso_command', help='CRISPResso command to call', default='CRISPResso')

        args = parser.parse_args()

        debug_flag = args.debug

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = {'name', 'output_folder'}
        crispresso_options_for_batch = list(crispresso_options-options_to_ignore)

        CRISPRessoShared.check_file(args.batch_settings)

        ##parse excel sheet
        batch_params=pd.read_csv(args.batch_settings, comment='#', sep='\t')
        #pandas either allows for auto-detect sep or for comment. not both
#        batch_params=pd.read_csv(args.batch_settings,sep=None,engine='python',error_bad_lines=False)
        batch_params.columns = batch_params.columns.str.strip(' -\xd0')

        #rename column "a" to "amplicon_seq", etc
        batch_params.rename(index=str, columns=CRISPRessoShared.get_crispresso_options_lookup(), inplace=True)
        batch_count = batch_params.shape[0]
        batch_params.index = range(batch_count)

        if 'fastq_r1' not in batch_params and 'bam_input' not in batch_params:
            raise CRISPRessoShared.BadParameterException("fastq_r1 must be specified in the batch settings file. Current headings are: "
                    + str(batch_params.columns.values))

        #add args from the command line to batch_params_df
        for arg in vars(args):
            if arg not in batch_params:
                batch_params[arg] = getattr(args, arg)
            else:
                if (getattr(args, arg) is not None):
                    batch_params[arg].fillna(value=getattr(args, arg), inplace=True)

        #assert that all names are unique
        #and clean names

        for i in range(batch_count):
            if batch_params.loc[i, 'name'] == '':
                batch_params.at[i, 'name'] = i
            batch_params.at[i, 'name'] = CRISPRessoShared.clean_filename(batch_params.loc[i, 'name'])

        if batch_params.drop_duplicates('name').shape[0] != batch_params.shape[0]:
            raise CRISPRessoShared.BadParameterException('Batch input names must be unique. The given names are not unique: ' + str(batch_params.loc[:, 'name']))

        #Check files
        batch_params["sgRNA_intervals"] = '' #create empty array for sgRNA intervals
        batch_params["sgRNA_intervals"] = batch_params["sgRNA_intervals"].apply(list)
        batch_params["cut_point_include_idx"] = '' #create empty array for cut point intervals for each batch based on sgRNA
        batch_params["cut_point_include_idx"] = batch_params["cut_point_include_idx"].apply(list)
        for idx, row in batch_params.iterrows():
            if 'fastq_r1' in row:
                if row.fastq_r1 is None:
                    raise CRISPRessoShared.BadParameterException("At least one fastq file must be given as a command line parameter or be specified in the batch settings file with the heading 'fastq_r1' (fastq_r1 on row %s '%s' is invalid)"%(int(idx)+1, row.fastq_r1))
                else:
                    CRISPRessoShared.check_file(row.fastq_r1)

            if 'fastq_r2' in row and row.fastq_r2 != "":
                CRISPRessoShared.check_file(row.fastq_r2)

            if 'input_bam' in row:
                if row.input_bam is None:
                    raise CRISPRessoShared.BadParameterException("At least one input file must be given as a command line parameter or be specified in the batch settings file with the heading 'fastq_r1' or 'input_bam' (input_bam on row %s '%s' is invalid)"%(int(idx)+1, row.input_bam))
                else:
                    CRISPRessoShared.check_file(row.input_bam)

            if args.auto:
                continue

            curr_amplicon_seq_str = row.amplicon_seq
            if curr_amplicon_seq_str is None:
                raise CRISPRessoShared.BadParameterException("Amplicon sequence must be given as a command line parameter or be specified in the batch settings file with the heading 'amplicon_seq' (Amplicon seq on row %s '%s' is invalid)"%(int(idx)+1, curr_amplicon_seq_str))

            guides_are_in_amplicon = {} #dict of whether a guide is in at least one amplicon sequence
            #iterate through amplicons
            for curr_amplicon_seq in str(curr_amplicon_seq_str).split(','):
                this_include_idxs=[] #mask for bp to include for this amplicon seq, as specified by sgRNA cut points
                this_sgRNA_intervals = []
                wrong_nt=CRISPRessoShared.find_wrong_nt(curr_amplicon_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The amplicon sequence in row %d (%s) contains incorrect characters:%s' % (idx+1, curr_amplicon_seq_str, ' '.join(wrong_nt)))

                #iterate through guides
                curr_guide_seq_string = row.guide_seq
                if curr_guide_seq_string is not None and curr_guide_seq_string != "":
                    guides = str(curr_guide_seq_string).strip().upper().split(',')
                    for curr_guide_seq in guides:
                        wrong_nt=CRISPRessoShared.find_wrong_nt(curr_guide_seq)
                        if wrong_nt:
                            raise CRISPRessoShared.NTException('The sgRNA sequence in row %d (%s) contains incorrect characters:%s'  % (idx+1, curr_guide_seq, ' '.join(wrong_nt)))
                    guide_mismatches = [[]]*len(guides)
                    guide_names = [""]*len(guides)
                    guide_qw_centers = CRISPRessoShared.set_guide_array(row.quantification_window_center, guides, 'guide quantification center')
                    guide_qw_sizes = CRISPRessoShared.set_guide_array(row.quantification_window_size, guides, 'guide quantification size')
                    guide_plot_cut_points = [1]*len(guides)
                    discard_guide_positions_overhanging_amplicon_edge = False
                    if 'discard_guide_positions_overhanging_amplicon_edge' in row:
                        discard_guide_positions_overhanging_amplicon_edge = row.discard_guide_positions_overhanging_amplicon_edge
                    (this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_sgRNA_names, this_include_idxs,
                        this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(curr_amplicon_seq, guides, guide_mismatches, guide_names, guide_qw_centers,
                        guide_qw_sizes, row.quantification_window_coordinates, row.exclude_bp_from_left, row.exclude_bp_from_right, row.plot_window_size, guide_plot_cut_points, discard_guide_positions_overhanging_amplicon_edge)
                    for guide_seq in this_sgRNA_sequences:
                        guides_are_in_amplicon[guide_seq] = 1

                batch_params.loc[idx, "cut_point_include_idx"].append(this_include_idxs)
                batch_params.loc[idx, "sgRNA_intervals"].append(this_sgRNA_intervals)

            for guide_seq in guides_are_in_amplicon:
                if guides_are_in_amplicon[guide_seq] != 1:
                    warn('\nThe guide sequence provided on row %d (%s) is not present in any amplicon sequence:%s! \nNOTE: The guide will be ignored for the analysis. Please check your input!' % (idx+1, row.guide_seq, curr_amplicon_seq))

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

        output_folder_name='CRISPRessoBatch_on_%s' % batch_folder_name
        OUTPUT_DIRECTORY=os.path.abspath(output_folder_name)

        if args.batch_output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.batch_output_folder), output_folder_name)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory

        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoBatch_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        with open(log_filename, 'w+') as outfile:
            outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Batch_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2Batch_info.json')
        crispresso2_info = {} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['version'] = CRISPRessoShared.__version__
        crispresso2_info['args'] = deepcopy(args)

        crispresso2_info['log_filename'] = os.path.basename(log_filename)

        #this is potentially the square-root of the input n_processes because sub-CRISPResso commands will take some processes as well
        args.n_processes = CRISPRessoShared.get_sub_n_processes(suppress_plots=args.suppress_plots, suppress_report=args.suppress_report, n_processes=args.n_processes)

        crispresso_cmd_to_write = ' '.join(sys.argv)
        if args.write_cleaned_report:
            cmd_copy = sys.argv[:]
            cmd_copy[0] = 'CRISPRessoBatch'
            for i in range(len(cmd_copy)):
                if os.sep in cmd_copy[i]:
                    cmd_copy[i] = os.path.basename(cmd_copy[i])

            crispresso_cmd_to_write = ' '.join(cmd_copy) #clean command doesn't show the absolute path to the executable or other files
        crispresso2_info['command_used'] = crispresso_cmd_to_write

        crispresso_cmds = []
        batch_names_arr = []
        batch_input_names = {}
        for idx, row in batch_params.iterrows():

            batch_name = CRISPRessoShared.slugify(row["name"])
            batch_names_arr.append(batch_name)
            batch_input_names[batch_name] = row["name"]

            crispresso_cmd= args.crispresso_command + ' -o %s --name %s' % (OUTPUT_DIRECTORY, batch_name)
            crispresso_cmd=propagate_options(crispresso_cmd, crispresso_options_for_batch, batch_params, idx)
            if row.amplicon_seq == "":
                crispresso_cmd += ' --auto '
            crispresso_cmds.append(crispresso_cmd)

        crispresso2_info['batch_names_arr'] = batch_names_arr
        crispresso2_info['batch_input_names'] = batch_input_names

        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, args.n_processes, 'batch', args.skip_failed)

        run_datas = [] #crispresso2 info from each row

        all_amplicons = set()
        amplicon_names = {}
        amplicon_counts = {}
        completed_batch_arr = []
        for idx, row in batch_params.iterrows():
            batch_name = CRISPRessoShared.slugify(row["name"])
            file_prefix = row['file_prefix']
            folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
            run_data_file = os.path.join(folder_name, 'CRISPResso2_info.json')
            if not os.path.isfile(run_data_file):
                info("Skipping folder '%s'. Cannot find run data at '%s'."%(folder_name, run_data_file))
                run_datas.append(None)
                continue

            run_data = CRISPRessoShared.load_crispresso_info(folder_name)
            run_datas.append(run_data)
            for ref_name in run_data['ref_names']:
                ref_seq = run_data['refs'][ref_name]['sequence']
                all_amplicons.add(ref_seq)
                #if this amplicon is called something else in another sample, just call it the amplicon
                if ref_name in amplicon_names and amplicon_names[ref_seq] != ref_name:
                    amplicon_names[ref_seq] = ref_seq
                else:
                    amplicon_names[ref_seq] = ref_name
                if ref_seq not in amplicon_counts:
                    amplicon_counts[ref_seq] = 0
                amplicon_counts[ref_seq]+= 1

            completed_batch_arr.append(batch_name)

        crispresso2_info['completed_batch_arr'] = completed_batch_arr

        #make sure amplicon names aren't super long
        for amplicon in all_amplicons:
            if len(amplicon_names[amplicon]) > 21:
                amplicon_names[amplicon] = amplicon_names[amplicon][0:21]

        #make sure no duplicate names (same name for the different amplicons)
        seen_names = {}
        for amplicon in all_amplicons:
            suffix_counter = 2
            orig_name = amplicon_names[amplicon]
            while amplicon_names[amplicon] in seen_names:
                amplicon_names[amplicon] = orig_name+"_"+str(suffix_counter)
                suffix_counter += 1
            seen_names[amplicon_names[amplicon]] = 1


        save_png = True
        if args.suppress_report:
            save_png = False

        window_nuc_pct_quilt_plot_names = []
        nuc_pct_quilt_plot_names = []
        window_nuc_conv_plot_names = []
        nuc_conv_plot_names = []

        #report for amplicons that appear multiple times
        for amplicon_index, amplicon_seq in enumerate(all_amplicons):
            #only perform comparison if amplicon seen in more than one sample
            if amplicon_counts[amplicon_seq] < 2:
                continue

            amplicon_name = amplicon_names[amplicon_seq]
            info('Reporting summary for amplicon: "' + amplicon_name + '"')

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
            batches_with_this_amplicon = []
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                file_prefix = row['file_prefix']
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
                run_data = run_datas[idx]
                if run_data is None:
                    continue
                batch_has_amplicon = False
                batch_amplicon_name = ''
                for ref_name in run_data['ref_names']:
                    if amplicon_seq == run_data['refs'][ref_name]['sequence']:
                        batch_has_amplicon = True
                        batch_amplicon_name = ref_name
                if not batch_has_amplicon:
                    continue
                batches_with_this_amplicon.append(idx)

                if consensus_guides == []:
                    consensus_guides = run_data['refs'][batch_amplicon_name]['sgRNA_sequences']
                    consensus_include_idxs = run_data['refs'][batch_amplicon_name]['include_idxs']
                    consensus_sgRNA_intervals = run_data['refs'][batch_amplicon_name]['sgRNA_intervals']
                    consensus_sgRNA_plot_idxs = run_data['refs'][batch_amplicon_name]['sgRNA_plot_idxs']

                if run_data['refs'][batch_amplicon_name]['sgRNA_sequences'] != consensus_guides:
                    guides_all_same = False
                if set(run_data['refs'][batch_amplicon_name]['include_idxs']) != set(consensus_include_idxs):
                    guides_all_same = False

                if 'nuc_freq_filename' not in run_data['refs'][batch_amplicon_name]:
                    info("Skipping the amplicon '%s' in folder '%s'. Cannot find nucleotide information."%(batch_amplicon_name, folder_name))
                    continue

                nucleotide_frequency_file = os.path.join(folder_name, run_data['refs'][batch_amplicon_name]['nuc_freq_filename'])
                ampSeq_nf, nuc_freqs = CRISPRessoShared.parse_count_file(nucleotide_frequency_file)

                nucleotide_pct_file = os.path.join(folder_name, run_data['refs'][batch_amplicon_name]['nuc_pct_filename'])
                ampSeq_np, nuc_pcts = CRISPRessoShared.parse_count_file(nucleotide_pct_file)

                count_file = os.path.join(folder_name, run_data['refs'][batch_amplicon_name]['mod_count_filename'])
                ampSeq_cf, mod_freqs = CRISPRessoShared.parse_count_file(count_file)

                if ampSeq_nf is None or ampSeq_np is None or ampSeq_cf is None:
                    info("Skipping the amplicon '%s' in folder '%s'. Could not parse batch output."%(batch_amplicon_name, folder_name))
                    info("Nucleotide frequency amplicon: '%s', Nucleotide percentage amplicon: '%s', Count vectors amplicon: '%s'"%(ampSeq_nf, ampSeq_np, ampSeq_cf))
                    continue
                if ampSeq_nf != ampSeq_np or ampSeq_np != ampSeq_cf:
                    warn("Skipping the amplicon '%s' in folder '%s'. Parsed amplicon sequences do not match\nnf:%s\nnp:%s\ncf:%s\nrf:%s"%(batch_amplicon_name, folder_name, ampSeq_nf, ampSeq_np, ampSeq_cf, amplicon_seq))
                    continue
                if consensus_sequence == "":
                    consensus_sequence = ampSeq_nf
                if ampSeq_nf != consensus_sequence:
                    info("Skipping the amplicon '%s' in folder '%s'. Amplicon sequences do not match."%(batch_amplicon_name, folder_name))
                    continue
                if 'Total' not in mod_freqs:
                    info("Skipping the amplicon '%s' in folder '%s'. Processing did not complete."%(batch_amplicon_name, folder_name))
                    continue
                if mod_freqs['Total'][0] == 0 or mod_freqs['Total'][0] == "0":
                    info("Skipping the amplicon '%s' in folder '%s'. Got no reads for amplicon."%(batch_amplicon_name, folder_name))
                    continue
                this_amp_total_reads = run_data['counts_total'][batch_amplicon_name]
                if this_amp_total_reads < args.min_reads_for_inclusion:
                    info("Skipping the amplicon '%s' in folder '%s'. Got %s reads (min_reads_for_inclusion is %d)."%(batch_amplicon_name, folder_name, str(this_amp_total_reads), args.min_reads_for_inclusion))
                    continue

                mod_pcts = {}
                for key in mod_freqs:
                    mod_pcts[key] = np.array(mod_freqs[key]).astype(np.float)/float(mod_freqs['Total'][0])

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
                info("Couldn't find any data for amplicon '%s'. Not compiling results."%amplicon_name)
            else:
                amplicon_plot_name = amplicon_name+"."
                if len(amplicon_names) == 1 and amplicon_name == "Reference":
                    amplicon_plot_name = ""

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

                crispresso2_info['nucleotide_frequency_summary_filename'] = os.path.basename(nucleotide_frequency_summary_filename)
                crispresso2_info['nucleotide_percentage_summary_filename'] = os.path.basename(nucleotide_percentage_summary_filename)

                crispresso2_info['modification_frequency_summary_filename'] = os.path.basename(modification_frequency_summary_filename)
                crispresso2_info['modification_percentage_summary_filename'] = os.path.basename(modification_percentage_summary_filename)

                crispresso2_info['summary_plot_titles'] = {}
                crispresso2_info['summary_plot_labels'] = {}
                crispresso2_info['summary_plot_datas'] = {}

                #if guides are all the same, merge substitutions and perform base editor comparison at guide quantification window
                if guides_all_same and consensus_guides != []:
                    info("All guides are equal. Performing comparison of batches for amplicon '%s'"% amplicon_name)
                    include_idxs = consensus_include_idxs #include indexes are the same for all guides
                    for idx, sgRNA in enumerate(consensus_guides):
                        sgRNA_intervals = consensus_sgRNA_intervals[idx]
                        sgRNA_plot_idxs = consensus_sgRNA_plot_idxs[idx]
                        plot_idxs_flat = [0, 1] # guide, nucleotide
                        plot_idxs_flat.extend([plot_idx + 2 for plot_idx in sgRNA_plot_idxs])
                        sub_nucleotide_frequency_summary_df = nucleotide_frequency_summary_df.iloc[:, plot_idxs_flat]
                        sub_nucleotide_percentage_summary_df = nucleotide_percentage_summary_df.iloc[:, plot_idxs_flat]
                        sub_modification_percentage_summary_df = modification_percentage_summary_df.iloc[:, plot_idxs_flat]

                        #show all sgRNA's on the plot
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

                        if not args.suppress_plots:
                            #plot for each guide
                            this_window_nuc_pct_quilt_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt_around_sgRNA_'+sgRNA)
                            CRISPRessoPlot.plot_nucleotide_quilt(sub_nucleotide_percentage_summary_df, sub_modification_percentage_summary_df, this_window_nuc_pct_quilt_plot_name, save_png, sgRNA_intervals=sub_sgRNA_intervals, quantification_window_idxs=include_idxs)
                            plot_name = os.path.basename(this_window_nuc_pct_quilt_plot_name)
                            window_nuc_pct_quilt_plot_names.append(plot_name)
                            crispresso2_info['summary_plot_titles'][plot_name] = 'sgRNA: ' + sgRNA + ' Amplicon: ' + amplicon_name
                            if len(consensus_guides) == 1:
                                crispresso2_info['summary_plot_titles'][plot_name] = ''
                            crispresso2_info['summary_plot_labels'][plot_name] = 'Composition of each base around the guide ' + sgRNA + ' for the amplicon ' + amplicon_name
                            crispresso2_info['summary_plot_datas'][plot_name] = [('Nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), ('Modification frequencies', os.path.basename(modification_frequency_summary_filename))]

                            sub_nucleotide_frequency_summary_df = pd.concat([sub_nucleotide_frequency_summary_df.iloc[:, 0:2],
                                                                        sub_nucleotide_frequency_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                            sub_nucleotide_frequency_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_frequency_summary_around_sgRNA_'+sgRNA+'.txt')
                            sub_nucleotide_frequency_summary_df.to_csv(sub_nucleotide_frequency_summary_filename, sep='\t', index=None)

                            sub_nucleotide_percentage_summary_df = pd.concat([sub_nucleotide_percentage_summary_df.iloc[:, 0:2],
                                                                        sub_nucleotide_percentage_summary_df.iloc[:, 2:].apply(pd.to_numeric)], axis=1)
                            sub_nucleotide_percentage_summary_filename = _jp(amplicon_plot_name + 'Nucleotide_percentage_summary_around_sgRNA_'+sgRNA+'.txt')
                            sub_nucleotide_percentage_summary_df.to_csv(sub_nucleotide_percentage_summary_filename, sep='\t', index=None)

                            if args.base_editor_output:
                                this_window_nuc_conv_plot_name = _jp(amplicon_plot_name + 'Nucleotide_conversion_map_around_sgRNA_'+sgRNA)
                                CRISPRessoPlot.plot_conversion_map(sub_nucleotide_percentage_summary_df, this_window_nuc_conv_plot_name, args.conversion_nuc_from, args.conversion_nuc_to, save_png, sgRNA_intervals=sub_sgRNA_intervals, quantification_window_idxs=include_idxs)
                                plot_name = os.path.basename(this_window_nuc_conv_plot_name)
                                window_nuc_conv_plot_names.append(plot_name)
                                crispresso2_info['summary_plot_titles'][plot_name] = 'sgRNA: ' + sgRNA + ' Amplicon: ' + amplicon_name
                                if len(consensus_guides) == 1:
                                    crispresso2_info['summary_plot_titles'][plot_name] = ''
                                crispresso2_info['summary_plot_labels'][plot_name] = args.conversion_nuc_from + '->' + args.conversion_nuc_to +' conversion rates around the guide ' + sgRNA + ' for the amplicon ' + amplicon_name
                                crispresso2_info['summary_plot_datas'][plot_name] = [('Nucleotide frequencies around sgRNA', os.path.basename(sub_nucleotide_frequency_summary_filename)),
                                                                                    ('Nucleotide percentages around sgRNA', os.path.basename(sub_nucleotide_percentage_summary_filename))
                                                                                ]

                    if not args.suppress_plots: # plot the whole region
                        this_nuc_pct_quilt_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt')
                        CRISPRessoPlot.plot_nucleotide_quilt(nucleotide_percentage_summary_df, modification_percentage_summary_df, this_nuc_pct_quilt_plot_name, save_png, sgRNA_intervals=consensus_sgRNA_intervals, quantification_window_idxs=include_idxs)
                        plot_name = os.path.basename(this_nuc_pct_quilt_plot_name)
                        nuc_pct_quilt_plot_names.append(plot_name)
                        crispresso2_info['summary_plot_titles'][plot_name] = 'Amplicon: ' + amplicon_name
                        if len(amplicon_names) == 1:
                            crispresso2_info['summary_plot_titles'][plot_name] = ''
                        crispresso2_info['summary_plot_labels'][plot_name] = 'Composition of each base for the amplicon ' + amplicon_name
                        crispresso2_info['summary_plot_datas'][plot_name] = [('Nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), ('Modification frequencies', os.path.basename(modification_frequency_summary_filename))]
                        if args.base_editor_output:
                            this_nuc_conv_plot_name = _jp(amplicon_plot_name + 'Nucleotide_conversion_map')
                            CRISPRessoPlot.plot_conversion_map(nucleotide_percentage_summary_df, this_nuc_conv_plot_name, args.conversion_nuc_from, args.conversion_nuc_to, save_png, sgRNA_intervals=consensus_sgRNA_intervals, quantification_window_idxs=include_idxs)
                            plot_name = os.path.basename(this_nuc_conv_plot_name)
                            nuc_conv_plot_names.append(plot_name)
                            crispresso2_info['summary_plot_titles'][plot_name] = 'Amplicon: ' + amplicon_name
                            if len(amplicon_names) == 1:
                                crispresso2_info['summary_plot_titles'][plot_name] = ''
                            crispresso2_info['summary_plot_titles'][plot_name] = ''
                            crispresso2_info['summary_plot_labels'][plot_name] = args.conversion_nuc_from + '->' + args.conversion_nuc_to +' conversion rates for the amplicon ' + amplicon_name
                            crispresso2_info['summary_plot_datas'][plot_name] = [('Nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), ('Modification frequencies', os.path.basename(modification_frequency_summary_filename))]

                else: #guides are not the same
                    if not args.suppress_plots:
                        this_nuc_pct_quilt_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt')
                        CRISPRessoPlot.plot_nucleotide_quilt(nucleotide_percentage_summary_df, modification_percentage_summary_df, this_nuc_pct_quilt_plot_name, save_png)
                        plot_name = os.path.basename(this_nuc_pct_quilt_plot_name)
                        nuc_pct_quilt_plot_names.append(plot_name)
                        crispresso2_info['summary_plot_labels'][plot_name] = 'Composition of each base for the amplicon ' + amplicon_name
                        crispresso2_info['summary_plot_datas'][plot_name] = [('Nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), ('Modification frequencies', os.path.basename(modification_frequency_summary_filename))]
                        if args.base_editor_output:
                            this_nuc_conv_plot_name = _jp(amplicon_plot_name + 'Nucleotide_percentage_quilt')
                            CRISPRessoPlot.plot_conversion_map(nucleotide_percentage_summary_df, this_nuc_conv_plot_name, args.conversion_nuc_from, args.conversion_nuc_to, save_png)
                            plot_name = os.path.basename(this_nuc_conv_plot_name)
                            nuc_conv_plot_names.append(plot_name)
                            crispresso2_info['summary_plot_labels'][plot_name] = args.conversion_nuc_from + '->' + args.conversion_nuc_to +' conversion rates for the amplicon ' + amplicon_name
                            crispresso2_info['summary_plot_datas'][plot_name] = [('Nucleotide frequencies', os.path.basename(nucleotide_frequency_summary_filename)), ('Modification frequencies', os.path.basename(modification_frequency_summary_filename))]


        crispresso2_info['window_nuc_pct_quilt_plot_names'] = window_nuc_pct_quilt_plot_names
        crispresso2_info['nuc_pct_quilt_plot_names'] = nuc_pct_quilt_plot_names
        crispresso2_info['window_nuc_conv_plot_names'] = window_nuc_conv_plot_names
        crispresso2_info['nuc_conv_plot_names'] = nuc_conv_plot_names

        #summarize amplicon modifications
        with open(_jp('CRISPRessoBatch_quantification_of_editing_frequency.txt'), 'w') as outfile:
            wrote_header = False
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                file_prefix = row['file_prefix']
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)
                run_data = run_datas[idx]
                if run_data is None:
                    continue

                amplicon_modification_file=os.path.join(folder_name, run_data['quant_of_editing_freq_filename'])
                with open(amplicon_modification_file, 'r') as infile:
                    file_head = infile.readline()
                    if not wrote_header:
                        outfile.write('Batch\t' + file_head)
                        wrote_header = True
                    for line in infile:
                        outfile.write(batch_name + "\t" + line)

        #summarize alignment
        with open(_jp('CRISPRessoBatch_mapping_statistics.txt'), 'w') as outfile:
            wrote_header = False
            for idx, row in batch_params.iterrows():
                batch_name = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % batch_name)

                run_data = run_datas[idx]
                if run_data is None:
                    continue
                amplicon_modification_file=os.path.join(folder_name, run_data['mapping_stats_filename'])
                with open(amplicon_modification_file, 'r') as infile:
                    file_head = infile.readline()
                    if not wrote_header:
                        outfile.write('Batch\t' + file_head)
                        wrote_header = True
                    for line in infile:
                        outfile.write(batch_name + "\t" + line)

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Batch_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_batch_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT)
            crispresso2_info['report_location'] = report_name
            crispresso2_info['report_filename'] = os.path.basename(report_name)

        end_time =  datetime.now()
        end_time_string =  end_time.strftime('%Y-%m-%d %H:%M:%S')
        running_time = end_time - start_time
        running_time_string =  str(running_time)

        crispresso2_info['end_time'] = end_time
        crispresso2_info['end_time_string'] = end_time_string
        crispresso2_info['running_time'] = running_time
        crispresso2_info['running_time_string'] = running_time_string

        CRISPRessoShared.write_crispresso_info(
            crispresso2Batch_info_file,
            crispresso2_info,
        )
        info('Analysis Complete!')
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

if __name__ == '__main__':
    main()
