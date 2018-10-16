# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import os
import errno
import sys
import subprocess as sb
import glob
import argparse
import unicodedata
import re
import string
import traceback
from CRISPResso import CRISPRessoShared
from CRISPResso import CRISPRessoPlot
from CRISPResso import CRISPRessoMultiProcessing
from CRISPResso import CRISPRessoReport

running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

if running_python3:
    import pickle as cp #python 3
else:
    import cPickle as cp #python 2.7


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
CRISPResso_to_call = os.path.join(os.path.dirname(_ROOT),'CRISPResso.py')

####Support functions###
def propagate_options(cmd,options,params,paramInd):
####
# cmd - the command to run
# options - list of options to propagate e.g. crispresso options
# params - df from excel parser e.g. params['amplicon_name'] = ['name1','name2','name3'] where each item corresponds to a different run in the batch
# paramInd - index in dict - this is the run number in the batch

    for option in options :
        if option:
            if option in params:
                val = params.loc[paramInd,option]
                if val is None:
                    pass
                elif str(val) == "True":
                    cmd+=' --%s' % option
                elif str(val) =="False":
                    pass
                elif type(val)==str:
                    if val != "":
                        if " " in val or "-" in val:
                            cmd+=' --%s "%s"' % (option,str(val)) # quotes for options with spaces
                        else:
                            cmd+=' --%s %s' % (option,str(val))
                elif type(val)==bool:
                    if val:
                        cmd+=' --%s' % option
                else:
                    cmd+=' --%s %s' % (option,str(val))
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
        description = ['~~~CRISPRessoBatch~~~','-Analysis of CRISPR/Cas9 outcomes from batch deep sequencing data-']
        batch_string = r'''
 _________________
| __    ___ __    |
||__) /\ | /  |__||
||__)/--\| \__|  ||
|_________________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description,batch_string))

        parser = CRISPRessoShared.getCRISPRessoArgParser(_ROOT, parserTitle = 'CRISPRessoBatch Parameters')

        #batch specific params
        parser.add_argument('-bs','--batch_settings', type=str, help='Settings file for batch. Must be tab-separated text file. The header row contains CRISPResso parameters (e.g., fastq_r1, fastq_r2, amplicon_seq, and other optional parameters). Each following row sets parameters for an additional batch.',required=True)
        parser.add_argument('--skip_failed',  help='Continue with batch analysis even if one sample fails',action='store_true')
        parser.add_argument('-p','--n_processes',type=int, help='Specify the number of processes to use for quantification.\
        Please use with caution since increasing this parameter will increase the memory required to run CRISPResso.',default=1)
        parser.add_argument('-bo','--batch_output_folder',  help='Directory where batch analysis output will be stored')

        args = parser.parse_args()

        debug_flag = args.debug

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = set(['output_folder'])
        crispresso_options_for_batch = list(crispresso_options-options_to_ignore)

        CRISPRessoShared.check_file(args.batch_settings)

        ##parse excel sheet
        batch_params=pd.read_csv(args.batch_settings,comment='#',sep='\t')
        #pandas either allows for auto-detect sep or for comment. not both
#        batch_params=pd.read_csv(args.batch_settings,sep=None,engine='python',error_bad_lines=False)
        batch_params.columns = batch_params.columns.str.strip(' -\xd0')

        #rename column "a" to "amplicon_seq", etc
        batch_params.rename(index=str,columns=CRISPRessoShared.get_crispresso_options_lookup(),inplace=True)
        batch_count = batch_params.shape[0]
        batch_params.index = range(batch_count)

        if 'fastq_r1' not in batch_params:
            raise CRISPRessoShared.BadParameterException("fastq_r1 must be specified in the batch settings file. Current headings are: "
                    + str(batch_params.columns.values))

        #add args from the command line to batch_params_df
        for arg in vars(args):
            if arg not in batch_params:
                batch_params[arg] = getattr(args,arg)
            else:
                if (getattr(args,arg) is not None):
                    batch_params[arg].fillna(value=getattr(args,arg), inplace=True)

        #assert that all names are unique

        for i in range(batch_count):
            if batch_params.loc[i,'name'] == '':
                batch_params.at[i,'name'] = i

        if batch_params.drop_duplicates('name').shape[0] != batch_params.shape[0]:
            raise CRISPRessoShared.BadParameterException('Batch input names must be unique. The given names are not unique: ' + str(batch_params.loc[:,'name']))

        quantification_window_center = -3 #default value
        if args.quantification_window_center is not None:
            quantification_window_center = args.quantification_window_center
        plot_window_size = 20 #default value
        if args.plot_window_size is not None:
            plot_window_size = args.plot_window_size


        #Check files
        batch_params["sgRNA_intervals"] = '' #create empty array for sgRNA intervals
        batch_params["sgRNA_intervals"] = batch_params["sgRNA_intervals"].apply(list)
        batch_params["cut_point_include_idx"] = '' #create empty array for cut point intervals for each batch based on sgRNA
        batch_params["cut_point_include_idx"] = batch_params["cut_point_include_idx"].apply(list)
        for idx,row in batch_params.iterrows():
            if row.fastq_r1 is None:
                raise CRISPRessoShared.BadParameterException("At least one fastq file must be given as a command line parameter or be specified in the batch settings file with the heading 'fastq_r1' (fastq_r1 on row %s '%s' is invalid)"%(int(idx)+1,row.fastq_r1))
            CRISPRessoShared.check_file(row.fastq_r1)

            if row.fastq_r2 != "":
                CRISPRessoShared.check_file(row.fastq_r2)

            curr_amplicon_seq_str = row.amplicon_seq
            if curr_amplicon_seq_str is None:
                raise CRISPRessoShared.BadParameterException("Amplicon sequence must be given as a command line parameter or be specified in the batch settings file with the heading 'amplicon_seq' (Amplicon seq on row %s '%s' is invalid)"%(int(idx)+1,curr_amplicon_seq_str))

            guides_are_in_amplicon = {} #dict of whether a guide is in at least one amplicon sequence
            #iterate through amplicons
            for curr_amplicon_seq in curr_amplicon_seq_str.split(','):
                this_include_idxs=[] #mask for bp to include for this amplicon seq, as specified by sgRNA cut points
                this_sgRNA_intervals = []
                wrong_nt=CRISPRessoShared.find_wrong_nt(curr_amplicon_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The amplicon sequence in row %d (%s) contains incorrect characters:%s' % (idx+1,curr_amplicon_seq_str,' '.join(wrong_nt)))

                #iterate through guides
                curr_guide_seq_string = row.guide_seq
                if curr_guide_seq_string is not None and curr_guide_seq_string != "":
                    guides = curr_guide_seq_string.strip().upper().split(',')
                    for curr_guide_seq in guides:
                        wrong_nt=CRISPRessoShared.find_wrong_nt(curr_guide_seq)
                        if wrong_nt:
                            raise CRISPRessoShared.NTException('The sgRNA sequence in row %d (%s) contains incorrect characters:%s'  % (idx+1,curr_guide_seq, ' '.join(wrong_nt)))
                    (this_sgRNA_sequences, this_sgRNA_intervals, this_cut_points, this_sgRNA_plot_offsets, this_include_idxs,
                        this_exclude_idxs, this_plot_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(curr_amplicon_seq,guides,row.quantification_window_center,
                        row.quantification_window_size,row.quantification_window_coordinates,row.exclude_bp_from_left,row.exclude_bp_from_right,row.plot_window_size)
                    for guide_seq in this_sgRNA_sequences:
                        guides_are_in_amplicon[guide_seq] = 1

                batch_params.ix[idx,"cut_point_include_idx"].append(this_include_idxs)
                batch_params.ix[idx,"sgRNA_intervals"].append(this_sgRNA_intervals)

            for guide_seq in guides_are_in_amplicon:
                if guides_are_in_amplicon[guide_seq] != 1:
                    warn('\nThe guide sequence provided on row %d (%s) is not present in any amplicon sequence:%s! \nNOTE: The guide will be ignored for the analysis. Please check your input!' % (idx+1,row.guide_seq,curr_amplicon_seq))

        batch_folder_name = os.path.splitext(os.path.basename(args.batch_settings))[0]

        OUTPUT_DIRECTORY='CRISPRessoBatch_on_%s' % batch_folder_name

        if args.batch_output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.batch_output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory

        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoBatch_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\nCRISPRessoBatch %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso_cmds = []
        for idx,row in batch_params.iterrows():
            batchName = row["name"]
            crispresso_cmd=CRISPResso_to_call + ' -o %s --name %s' % (OUTPUT_DIRECTORY,batchName)
            crispresso_cmd=propagate_options(crispresso_cmd,crispresso_options_for_batch,batch_params,idx)
            crispresso_cmds.append(crispresso_cmd)

        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,args.n_processes,'batch',args.skip_failed)

        run_datas = [] #crispresso2 info from each row

        for idx,row in batch_params.iterrows():
            batchName = row["name"]
            file_prefix = row['file_prefix']
            folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)
            run_data = cp.load(open(os.path.join(folder_name,'CRISPResso2_info.pickle'),'rb'))
            run_datas.append(run_data)

        save_png = True
        if args.suppress_report:
            save_png = False

        #if amplicons are all the same, merge substitutions and perform base editor comparison
        if batch_params.drop_duplicates('amplicon_seq').shape[0]  == 1:
            info("All amplicons are equal. Performing comparison of batches")

            def report_nucleotide_summary(amplicon_seq,amplicon_name,amplicon_index):
                consensus_sequence = ""
                nucleotide_frequency_summary = []
                nucleotide_percentage_summary = []
                modification_frequency_summary = []
                modification_percentage_summary = []

                amp_found_count = 0 #how many folders had information for this amplicon
                for idx,row in batch_params.iterrows():
                    batchName = row["name"]
                    file_prefix = row['file_prefix']
                    folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)
                    run_data = run_datas[idx]

                    if 'nuc_freq_filename' not in run_data['refs'][amplicon_name]:
                        info("Skipping the amplicon '%s' in folder '%s'. Cannot find nucleotide information."%(amplicon_name,folder_name))
                        continue

                    nucleotide_frequency_file = run_data['refs'][amplicon_name]['nuc_freq_filename']
                    ampSeq_nf,nuc_freqs = CRISPRessoShared.parse_count_file(nucleotide_frequency_file)

                    nucleotide_pct_file=run_data['refs'][amplicon_name]['nuc_pct_filename']
                    ampSeq_np,nuc_pcts = CRISPRessoShared.parse_count_file(nucleotide_pct_file)

                    count_file=run_data['refs'][amplicon_name]['mod_count_filename']
                    ampSeq_cf,mod_freqs = CRISPRessoShared.parse_count_file(count_file)

                    if ampSeq_nf is None or ampSeq_np is None or ampSeq_cf is None:
                        info("Skipping the amplicon '%s' in folder '%s'. Could not parse batch output."%(amplicon_name,folder_name))
                        info("Nucleotide frequency amplicon: '%s', Nucleotide percentage amplicon: '%s', Count vectors amplicon: '%s'"%(ampSeq_nf,ampSeq_np,ampSeq_cf))
                        continue
                    if ampSeq_nf != ampSeq_np or ampSeq_np != ampSeq_cf:
                        warn("Skipping the amplicon '%s' in folder '%s'. Parsed amplicon sequences do not match\nnf:%s\nnp:%s\ncf:%s\nrf:%s"%(amplicon_name,folder_name,ampSeq_nf,ampSeq_np,ampSeq_cf,amplicon_seq))
                        continue
                    if consensus_sequence == "":
                        consensus_sequence = ampSeq_nf
                    if ampSeq_nf != consensus_sequence:
                        info("Skipping the amplicon '%s' in folder '%s'. Amplicon sequences do not match."%(amplicon_name,folder_name))
                        continue
                    if 'Total' not in mod_freqs:
                        info("Skipping the amplicon '%s' in folder '%s'. Processing did not complete."%(amplicon_name,folder_name))
                        continue
                    if mod_freqs['Total'][0] == 0 or mod_freqs['Total'][0] == "0":
                        info("Skipping the amplicon '%s' in folder '%s'. Got no reads for amplicon."%(amplicon_name,folder_name))
                        continue

                    mod_pcts = {}
                    for key in mod_freqs:
                        mod_pcts[key] = np.array(mod_freqs[key]).astype(np.float)/float(mod_freqs['Total'][0])

                    amp_found_count += 1

                    for nuc in ['A','T','C','G','N','-']:
                        row = [batchName,nuc]
                        row.extend(nuc_freqs[nuc])
                        nucleotide_frequency_summary.append(row)

                        pct_row = [batchName,nuc]
                        pct_row.extend(nuc_pcts[nuc])
                        nucleotide_percentage_summary.append(pct_row)

                    for mod in ['Insertions','Insertions_Left','Deletions','Substitutions','All_modifications']:
                        row = [batchName,mod]
                        row.extend(mod_freqs[mod])
                        modification_frequency_summary.append(row)

                        pct_row = [batchName,mod]
                        pct_row.extend(mod_pcts[mod])
                        modification_percentage_summary.append(pct_row)

                if amp_found_count == 0:
                    info("Couldn't find any data for amplicon '%s'. Not compiling results."%amplicon_name)
                    return()

                colnames = ['Batch','Nucleotide']
                colnames.extend(list(consensus_sequence))
                nucleotide_frequency_summary_df = pd.DataFrame(nucleotide_frequency_summary,columns=colnames)
                nucleotide_frequency_summary_df = pd.concat([nucleotide_frequency_summary_df.iloc[:,0:2],
                                                            nucleotide_frequency_summary_df.iloc[:,2:].apply(pd.to_numeric)],axis=1)
                nucleotide_frequency_summary_df.to_csv(_jp(amplicon_name + '.NUCLEOTIDE_FREQUENCY_SUMMARY.txt'),sep='\t',index=None)

                nucleotide_percentage_summary_df = pd.DataFrame(nucleotide_percentage_summary,columns=colnames)
                nucleotide_percentage_summary_df = pd.concat([nucleotide_percentage_summary_df.iloc[:,0:2],
                                                        nucleotide_percentage_summary_df.iloc[:,2:].apply(pd.to_numeric)],axis=1)
                nucleotide_percentage_summary_df.to_csv(_jp(amplicon_name + '.NUCLEOTIDE_PERCENTAGE_SUMMARY.txt'),sep='\t',index=None)

                colnames = ['Batch','Modification']
                colnames.extend(list(consensus_sequence))
                modification_frequency_summary_df = pd.DataFrame(modification_frequency_summary,columns=colnames)
                modification_frequency_summary_df = pd.concat([modification_frequency_summary_df.iloc[:,0:2],
                                                            modification_frequency_summary_df.iloc[:,2:].apply(pd.to_numeric)],axis=1)
                modification_frequency_summary_df.to_csv(_jp(amplicon_name + '.MODIFICATION_FREQUENCY_SUMMARY.txt'),sep='\t',index=None)

                modification_percentage_summary_df = pd.DataFrame(modification_percentage_summary,columns=colnames)
                modification_percentage_summary_df = pd.concat([modification_percentage_summary_df.iloc[:,0:2],
                                                        modification_percentage_summary_df.iloc[:,2:].apply(pd.to_numeric)],axis=1)
                modification_percentage_summary_df.to_csv(_jp(amplicon_name + '.MODIFICATION_PERCENTAGE_SUMMARY.txt'),sep='\t',index=None)



                #if guides are all the same, merge substitutions and perform base editor comparison at guide quantification window
                if batch_params.drop_duplicates('guide_seq').shape[0]  == 1 and batch_params.ix[0,"cut_point_include_idx"][amplicon_index]:
                    include_idxs = batch_params.ix[0,"cut_point_include_idx"][amplicon_index]
                    sgRNA_intervals = batch_params.ix[0,"sgRNA_intervals"][amplicon_index]
                    info("All guides are equal. Performing comparison of batches for amplicon '%s'"% amplicon_name)
                    include_idxs_flat = [0,1] # guide, nucleotide
                    include_idxs_flat.extend([cutidx + 2 for cutidx in include_idxs])
                    sub_nucleotide_frequency_summary_df = nucleotide_frequency_summary_df.iloc[:,include_idxs_flat]
                    sub_nucleotide_percentage_summary_df = nucleotide_percentage_summary_df.iloc[:,include_idxs_flat]
                    sub_modification_percentage_summary_df = modification_percentage_summary_df.iloc[:,include_idxs_flat]
                    sub_sgRNA_intervals = []
                    for sgRNA_interval in sgRNA_intervals:
                        newstart = None
                        newend = None
                        for idx,i in enumerate(include_idxs):
                            if i <= sgRNA_interval[0]:
                                newstart = idx
                            if newend is None and i >= sgRNA_interval[1]:
                                newend = idx

                        #if guide doesn't overlap with include indexes
                        if newend == 0 or newstart == len(include_idxs):
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
                        sub_sgRNA_intervals.append((newstart,newend))

                    CRISPRessoPlot.plot_nucleotide_quilt(sub_nucleotide_percentage_summary_df,sub_modification_percentage_summary_df,_jp(amplicon_name
                        + '.Quantification_Window_Nucleotide_Percentage_Quilt'),save_png,sgRNA_intervals=sub_sgRNA_intervals)
                    if args.base_editor_output:
                        CRISPRessoPlot.plot_conversion_map(sub_nucleotide_percentage_summary_df,_jp(amplicon_name + '.Quantification_Window_Nucleotide_Conversion'),args.conversion_nuc_from,args.conversion_nuc_to,save_png,sgRNA_intervals=sub_sgRNA_intervals)

                    CRISPRessoPlot.plot_nucleotide_quilt(nucleotide_percentage_summary_df,modification_percentage_summary_df,_jp(amplicon_name + '.Nucleotide_Percentage_Quilt'),save_png,sgRNA_intervals=sgRNA_intervals,quantification_window_idxs=include_idxs)
                    if args.base_editor_output:
                        CRISPRessoPlot.plot_conversion_map(nucleotide_percentage_summary_df,_jp(amplicon_name + '.Nucleotide_Conversion'),args.conversion_nuc_from,args.conversion_nuc_to,save_png,sgRNA_intervals=sgRNA_intervals)
                else: #guides are not the same
                    CRISPRessoPlot.plot_nucleotide_quilt(nucleotide_percentage_summary_df,modification_percentage_summary_df,_jp(amplicon_name + '.Nucleotide_Percentage_Quilt'),save_png)
                    if args.base_editor_output:
                        CRISPRessoPlot.plot_conversion_map(nucleotide_percentage_summary_df,_jp(amplicon_name + '.Nucleotide_Conversion'),args.conversion_nuc_from,args.conversion_nuc_to,save_png)

            amplicon_seqs = row.amplicon_seq.split(',')
            amplicon_names = row.amplicon_name.split(',')
            for i in range(len(amplicon_seqs)):
                this_amplicon_seq = amplicon_seqs[i]
                this_amplicon_name = amplicon_names[i]
                info('Plotting for amplicon ' + str(i+1) + ": '" + this_amplicon_name + "'" )
                report_nucleotide_summary(this_amplicon_seq,this_amplicon_name,i)


            #summarize amplicon modifications
            with open(_jp('CRISPRessoBatch_quantification_of_editing_frequency.txt'),'w') as outfile:
                wrote_header = False
                for idx,row in batch_params.iterrows():
                    batchName = row["name"]
                    file_prefix = row['file_prefix']
                    folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)
                    run_data = run_datas[idx]

                    amplicon_modification_file=run_data['quant_of_editing_freq_filename']
                    with open(amplicon_modification_file,'r') as infile:
                        file_head = infile.readline()
                        if not wrote_header:
                            outfile.write('Batch\t' + file_head)
                            wrote_header = True
                        for line in infile:
                            outfile.write(batchName + "\t" + line)

            #summarize alignment
            with open(_jp('CRISPRessoBatch_mapping_statistics.txt'),'w') as outfile:
                wrote_header = False
                for idx,row in batch_params.iterrows():
                    batchName = row["name"]
                    folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)

                    run_data = run_datas[idx]
                    amplicon_modification_file=run_data['mapping_stats_filename']
                    with open(amplicon_modification_file,'r') as infile:
                        file_head = infile.readline()
                        if not wrote_header:
                            outfile.write('Batch\t' + file_head)
                            wrote_header = True
                        for line in infile:
                            outfile.write(batchName + "\t" + line)

        if not args.suppress_report:
            report_name = _jp('CRISPResso2Batch_report.html')
            CRISPRessoReport.make_batch_report_from_folder(report_name,OUTPUT_DIRECTORY,_ROOT)

        info('All Done!')
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
