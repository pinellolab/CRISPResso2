# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import os
from copy import deepcopy
import sys
import argparse
import re
import traceback
import json
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2.CRISPRessoReports import CRISPRessoReport


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info

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
                error('You need to install %s module to use CRISPRessoMeta!' % library_name)
                sys.exit(1)



pd=check_library('pandas')
np=check_library('numpy')

def main():
    try:
        description = ['~~~CRISPRessoMeta~~~', '-Analysis of CRISPR/Cas9 outcomes from deep sequencing data using a metadata file-']
        meta_string = r'''
 ________________________________________
|   _________   ______ _______  ______   |
|  | | | | | \ | |       | |   | |  | |  |
|  | | | | | | | |----   | |   | |__| |  |
|  |_| |_| |_| |_|____   |_|   |_|  |_|  |
|________________________________________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description, meta_string))

        parser = CRISPRessoShared.getCRISPRessoArgParser("Meta", parser_title = 'CRISPRessoMeta Parameters')

        #batch specific params
        parser.add_argument('--metadata', type=str, help='Metadata file according to NIST specification', required=True)
        parser.add_argument('-mo', '--meta_output_folder',  help='Directory where analysis output will be stored')
        parser.add_argument('--crispresso_command', help='CRISPResso command to call', default='CRISPResso')

        args = parser.parse_args()

        if args.use_matplotlib or not CRISPRessoShared.is_C2Pro_installed():
            from CRISPResso2 import CRISPRessoPlot
        else:
            from CRISPRessoPro import plot as CRISPRessoPlot

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

        debug_flag = args.debug

        crispresso_options = CRISPRessoShared.get_core_crispresso_options()
        options_to_ignore = {'name', 'output_folder'}
        crispresso_options_for_meta = list(crispresso_options-options_to_ignore)

        CRISPRessoShared.check_file(args.metadata)

        meta_params = pd.DataFrame(columns=['name', 'guide_seq', 'amplicon_seq'])
        with open(args.metadata) as metadata_file:
            metadata = json.load(metadata_file)

            exp = metadata['Experiment']
            for guide in data['Experiment']:
                print('Guide: ' + guide['name'])
                print('Sequence: ' + guide['sequence'])
                print('Amplicon: ' + guide['amplicon'])
                print('Fastq_R1: ' + guide['fastq_r1'])
                print('Fastq_R2: ' + guide['fastq_r2'])
                meta_params.append({'name':guide['name'],'guide_seq':guide['sequence'],'amplicon_seq':guide['amplicon'],'fastq_r1':guide['fastq_r1'],'fastq_r2':guide['fastq_r2']})


        print('table:')
        print(meta_params)
        #rename column "a" to "amplicon_seq", etc
        meta_params.rename(index=str, columns=CRISPRessoShared.get_crispresso_options_lookup("Core"), inplace=True)
        meta_count = meta_params.shape[0]
        meta_params.index = range(meta_count)

        if 'fastq_r1' not in meta_params:
            raise CRISPRessoShared.BadParameterException("fastq_r1 must be specified in the meta settings file. Current headings are: "
                    + str(meta_params.columns.values))

        #add args from the command line to meta_params
        for arg in vars(args):
            if arg not in meta_params:
                meta_params[arg] = getattr(args, arg)
            else:
                if (getattr(args, arg) is not None):
                    meta_params[arg].fillna(value=getattr(args, arg), inplace=True)

        #assert that all names are unique
        #and clean names

        for i in range(meta_count):
            if meta_params.loc[i, 'name'] == '':
                meta_params.at[i, 'name'] = i
            meta_params.at[i, 'name'] = CRISPRessoShared.clean_filename(meta_params.loc[i, 'name'])

        if meta_params.drop_duplicates('name').shape[0] != meta_params.shape[0]:
            raise CRISPRessoShared.BadParameterException('Sample input names must be unique. The given names are not unique: ' + str(meta_params.loc[:, 'name']))

        #Check files
        meta_params["sgRNA_intervals"] = '' #create empty array for sgRNA intervals
        meta_params["sgRNA_intervals"] = meta_params["sgRNA_intervals"].apply(list)
        meta_params["cut_point_include_idx"] = '' #create empty array for cut point intervals for each batch based on sgRNA
        meta_params["cut_point_include_idx"] = meta_params["cut_point_include_idx"].apply(list)
        for idx, row in meta_params.iterrows():
            if row.fastq_r1 is None:
                raise CRISPRessoShared.BadParameterException("At least one fastq file must be given as a command line parameter or be specified in the meta settings file with the heading 'fastq_r1' (fastq_r1 on row %s '%s' is invalid)"%(int(idx)+1, row.fastq_r1))
            CRISPRessoShared.check_file(row.fastq_r1)

            if row.fastq_r2 != "":
                CRISPRessoShared.check_file(row.fastq_r2)

            if args.auto:
                continue

            curr_amplicon_seq_str = row.amplicon_seq
            if curr_amplicon_seq_str is None:
                raise CRISPRessoShared.BadParameterException("Amplicon sequence must be given as a command line parameter or be specified in the meta settings file with the heading 'amplicon_seq' (Amplicon seq on row %s '%s' is invalid)"%(int(idx)+1, curr_amplicon_seq_str))

            guides_are_in_amplicon = {} #dict of whether a guide is in at least one amplicon sequence
            #iterate through amplicons
            for curr_amplicon_seq in curr_amplicon_seq_str.split(','):
                this_include_idxs=[] #mask for bp to include for this amplicon seq, as specified by sgRNA cut points
                this_sgRNA_intervals = []
                wrong_nt=CRISPRessoShared.find_wrong_nt(curr_amplicon_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The amplicon sequence in row %d (%s) contains incorrect characters:%s' % (idx+1, curr_amplicon_seq_str, ' '.join(wrong_nt)))

                #iterate through guides
                curr_guide_seq_string = row.guide_seq
                if curr_guide_seq_string is not None and curr_guide_seq_string != "":
                    guides = curr_guide_seq_string.strip().upper().split(',')
                    for curr_guide_seq in guides:
                        wrong_nt=CRISPRessoShared.find_wrong_nt(curr_guide_seq)
                        if wrong_nt:
                            raise CRISPRessoShared.NTException('The sgRNA sequence in row %d (%s) contains incorrect characters:%s'  % (idx+1, curr_guide_seq, ' '.join(wrong_nt)))
                    guide_names = ['']*len(guides)
                    guide_mismatches = [[]]*len(guides)
                    guide_qw_centers = CRISPRessoShared.set_guide_array(row.quantification_window_center, guides, 'guide quantification center')
                    guide_qw_sizes = CRISPRessoShared.set_guide_array(row.quantification_window_size, guides, 'guide quantification size')
                    guide_plot_cut_points = [1]*len(guides)
                    discard_guide_positions_overhanging_amplicon_edge = False
                    if 'discard_guide_positions_overhanging_amplicon_edge' in row:
                        discard_guide_positions_overhanging_amplicon_edge = row.discard_guide_positions_overhanging_amplicon_edge

                    (this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_sgRNA_names, this_sgRNA_include_idxs, this_include_idxs,
                        this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(curr_amplicon_seq, guides, guide_mismatches, guide_names, guide_qw_centers,
                        guide_qw_sizes, row.quantification_window_coordinates, row.exclude_bp_from_left, row.exclude_bp_from_right, row.plot_window_size, guide_plot_cut_points, discard_guide_positions_overhanging_amplicon_edge)
                    for guide_seq in this_sgRNA_sequences:
                        guides_are_in_amplicon[guide_seq] = 1

                meta_params.ix[idx, "cut_point_include_idx"].append(this_include_idxs)
                meta_params.ix[idx, "sgRNA_intervals"].append(this_sgRNA_intervals)

            for guide_seq in guides_are_in_amplicon:
                if guides_are_in_amplicon[guide_seq] != 1:
                    warn('\nThe guide sequence provided on row %d (%s) is not present in any amplicon sequence:%s! \nNOTE: The guide will be ignored for the analysis. Please check your input!' % (idx+1, row.guide_seq, curr_amplicon_seq))

        meta_folder_name = os.path.splitext(os.path.basename(args.metadata))[0]
        if args.name and args.name != "":
            meta_folder_name = args.name

        output_folder_name='CRISPRessoMeta_on_%s' % meta_folder_name
        OUTPUT_DIRECTORY=os.path.abspath(output_folder_name)

        if args.meta_output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.meta_output_folder), output_folder_name)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory

        try:
            info('Creating Folder %s' % OUTPUT_DIRECTORY, {'percent_complete': 0})
            os.makedirs(OUTPUT_DIRECTORY)
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoMeta_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        logger.addHandler(CRISPRessoShared.StatusHandler('CRISPRessoMeta_status.json'))

        with open(log_filename, 'w+') as outfile:
            outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Meta_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2Meta_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        crispresso_cmds = []
        meta_names_arr = []
        meta_input_names = {}
        for idx, row in meta_params.iterrows():

            metaName = CRISPRessoShared.slugify(row["name"])
            meta_names_arr.append(metaName)
            meta_input_names[metaName] = row["name"]

            crispresso_cmd= args.crispresso_command + ' -o %s --name %s' % (OUTPUT_DIRECTORY, metaName)
            crispresso_cmd=propagate_options(crispresso_cmd, crispresso_options_for_meta, meta_params, idx)
            crispresso_cmds.append(crispresso_cmd)

        crispresso2_info['meta_names_arr'] = meta_names_arr
        crispresso2_info['meta_input_names'] = meta_input_names

        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, args.n_processes, 'meta', args.skip_failed, start_end_percent=(10, 90))

        run_datas = [] #crispresso2 info from each row

        all_amplicons = set()
        amplicon_names = {}
        amplicon_counts = {}
        completed_meta_arr = []
        for idx, row in meta_params.iterrows():
            metaName = CRISPRessoShared.slugify(row["name"])
            folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % metaName)
            run_data_file = os.path.join(folder_name, 'CRISPResso2_info.json')
            if os.path.isfile(run_data_file) is False:
                info("Skipping folder '%s'. Cannot find run data at '%s'."%(folder_name, run_data_file))
                run_datas.append(None)
                continue

            run_data = CRISPRessoShared.load_crispresso_info(folder_name)
            run_datas.append(run_data)
            for ref_name in run_data['results']['ref_names']:
                ref_seq = run_data['results']['refs'][ref_name]['sequence']
                all_amplicons.add(ref_seq)
                #if this amplicon is called something else in another sample, just call it the amplicon
                if ref_name in amplicon_names and amplicon_names[ref_seq] != ref_name:
                    amplicon_names[ref_seq] = ref_seq
                else:
                    amplicon_names[ref_seq] = ref_name
                if ref_seq not in amplicon_counts:
                    amplicon_counts[ref_seq] = 0
                amplicon_counts[ref_seq]+= 1

            completed_meta_arr.append(metaName)

        crispresso2_info['completed_meta_arr'] = completed_meta_arr

        #make sure amplicon names aren't super long
        for amplicon in all_amplicons:
            if len(amplicon_names[amplicon]) > 20:
                amplicon_names[amplicon] = amplicon_names[amplicon][0:20]

        #make sure no duplicate names (same name for the different amplicons)
        seen_names = {}
        for amplicon in all_amplicons:
            suffix_counter = 2
            while amplicon_names[amplicon] in seen_names:
                amplicon_names[amplicon] = amplicon_names[amplicon]+"_"+str(suffix_counter)
                suffix_counter += 1
            seen_names[amplicon_names[amplicon]] = 1


        save_png = True
        if args.suppress_report:
            save_png = False

        #summarize amplicon modifications
        with open(_jp('CRISPRessoBatch_quantification_of_editing_frequency.txt'), 'w') as outfile:
            wrote_header = False
            for idx, row in meta_params.iterrows():
                metaName = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % metaName)
                run_data = run_datas[idx]
                if run_data is None:
                    continue

                amplicon_modification_file=os.path.join(folder_name, run_data['running_info']['quant_of_editing_freq_filename'])
                with open(amplicon_modification_file, 'r') as infile:
                    file_head = infile.readline()
                    if not wrote_header:
                        outfile.write('Batch\t' + file_head)
                        wrote_header = True
                    for line in infile:
                        outfile.write(metaName + "\t" + line)

        #summarize alignment
        with open(_jp('CRISPRessoBatch_mapping_statistics.txt'), 'w') as outfile:
            wrote_header = False
            for idx, row in meta_params.iterrows():
                metaName = CRISPRessoShared.slugify(row["name"])
                folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % metaName)

                run_data = run_datas[idx]
                if run_data is None:
                    continue
                amplicon_modification_file=os.path.join(folder_name, run_data['running_info']['mapping_stats_filename'])
                with open(amplicon_modification_file, 'r') as infile:
                    file_head = infile.readline()
                    if not wrote_header:
                        outfile.write('Batch\t' + file_head)
                        wrote_header = True
                    for line in infile:
                        outfile.write(metaName + "\t" + line)

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Meta_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_meta_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT, logger)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        CRISPRessoShared.write_crispresso_info(
            crispresso2Meta_info_file, crispresso2_info,
        )
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

if __name__ == '__main__':
    main()
