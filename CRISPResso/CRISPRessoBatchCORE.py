# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2017 The General Hospital Corporation. All Rights Reserved.
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
import multiprocessing as mp
import signal
from functools import partial
import CRISPRessoShared
import CRISPRessoPlot


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


def runCrispresso(OUTPUT_DIRECTORY, selected_options,batch_params,idx):
    batchName = batch_params.loc[idx,'name']

    testing = True
    if (testing):
        curr_file_loc =  str(os.path.abspath(__file__))
        if '2017_07_CRISPRESSO_DOS' in curr_file_loc:
            warn('Need to update this with installed Crispresso')
            crispresso_cmd='python /data/pinello/PROJECTS/2017_07_CRISPRESSO_DOS/src/CRISPResso.py -o %s --name %s' % (OUTPUT_DIRECTORY,batchName)
        else:
            crispresso_cmd='/opt/conda/bin/python /CRISPResso/CRISPResso.py -o %s --name %s' % (OUTPUT_DIRECTORY,batchName)
    else:
        crispresso_cmd='CRISPResso -o %s --name %s' % (OUTPUT_DIRECTORY,batchName)

    crispresso_cmd=propagate_options(crispresso_cmd,selected_options,batch_params,idx)
    info('Running CRISPResso batch #%d: %s' % (idx,crispresso_cmd))

    #don't actually process -- just use already processed files (for debug)
    runFake = False

    if (runFake):
        return_value = sb.call("sleep 1",shell=True)
    else:
        return_value = sb.call(crispresso_cmd,shell=True)


    if return_value != 0:
        warn('CRISPResso command failed (return value ' + str(return_value) + ') on batch #' + str(idx) + ": " + crispresso_cmd)
    else:
        info('Finished CRISPResso batch #%d' % idx)
    return return_value

def get_data(path):
        return os.path.join(_ROOT, 'data', path)

def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoBatch!' % library_name)
                sys.exit(1)



pd=check_library('pandas')
np=check_library('numpy')

###EXCEPTIONS############################
class FlashException(Exception):
    pass

class TrimmomaticException(Exception):
    pass

class Bowtie2Exception(Exception):
    pass

class AmpliconsNotUniqueException(Exception):
    pass

class AmpliconsNamesNotUniqueException(Exception):
    pass

class NoReadsAlignedException(Exception):
    pass

class DonorSequenceException(Exception):
    pass

class AmpliconEqualDonorException(Exception):
    pass

class SgRNASequenceException(Exception):
    pass

class NTException(Exception):
    pass

class ExonSequenceException(Exception):
    pass

class BadParameterException(Exception):
    pass


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
        parser.add_argument('--batch_settings', type=str, help='Settings file for batch. Must be tab-separated text file. The header row contains CRISPResso parameters (e.g., fastq_r1, fastq_r2, amplicon_seq, and other optional parameters). Each following row sets parameters for an additional batch.',required=True)
        parser.add_argument('-p','--n_processes',type=int, help='Specify the number of processes to use for the quantification.\
        Please use with caution since increasing this parameter will increase the memory required to run CRISPResso.',default=1)
        parser.add_argument('-bo','--batch_output_folder',  help='Directory where batch analysis output will be stored')

        args = parser.parse_args()

        debug_flag = args.debug

        crispresso_options_for_batch=['fastq_r1','fastq_r2','amplicon_seq',
    		'amplicon_name', 'amplicon_min_alignment_score',
    		'amplicon_min_unmodified_score',
    		'guide_seq', 'coding_seq',
    		'min_average_read_quality', 'min_single_bp_quality','min_bp_quality_or_N',
    		'name',
            #'output_folder', #disable setting of output folder
    		'split_paired_end', 'trim_sequences', 'trimmomatic_options_string',
    		'min_paired_end_reads_overlap', 'max_paired_end_reads_overlap',
    		'window_around_sgrna', 'cleavage_offset',
    		'exclude_bp_from_left', 'exclude_bp_from_right',
    		'ignore_substitutions', 'ignore_insertions', 'ignore_deletions',
    		'needleman_wunsch_gap_open', 'needleman_wunsch_gap_extend','needleman_wunsch_gap_incentive',
    		'keep_intermediate', 'dump', 'save_also_png',
    		'offset_around_cut_to_plot',
    		'min_frequency_alleles_around_cut_to_plot',
    		'max_rows_alleles_around_cut_to_plot',
    		'default_min_aln_score',
            'conversion_nuc_from','conversion_nuc_to','base_editor_mode','analysis_window_coordinates','crispresso1_mode','debug'
        ]

        ##parse excel sheet
        batch_params=pd.read_csv(args.batch_settings,comment='#',sep='\t')
        #pandas either allows for auto-detect sep or for comment. not both
#        batch_params=pd.read_csv(args.batch_settings,sep=None,engine='python',error_bad_lines=False)
        batch_params.columns = batch_params.columns.str.strip(' -\xd0')

        #rename column "a" to "amplicon_seq", etc
        batch_params.rename(index=str,columns=CRISPRessoShared.crispresso_options_lookup,inplace=True)
        batch_count = batch_params.shape[0]
        batch_params.index = range(batch_count)

        if 'fastq_r1' not in batch_params:
            raise BadParameterException("fastq_r1 must be specified in the batch settings file. Current headings are: "
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
            raise Exception('Batch input names must be unique. The given names are not unique: ' + str(batch_params.loc[:,'name']))

        cleavage_offset = -3 #default value
        if args.cleavage_offset is not None:
            cleavage_offset = args.cleavage_offset
        offset_around_cut_to_plot = 20 #default value
        if args.offset_around_cut_to_plot is not None:
            offset_around_cut_to_plot = args.offset_around_cut_to_plot


        #Check files
        batch_params["sgRNA_intervals"] = '' #create empty array for sgRNA intervals
        batch_params["sgRNA_intervals"] = batch_params["sgRNA_intervals"].apply(list)
        batch_params["cut_point_include_idx"] = '' #create empty array for cut point intervals for each batch based on sgRNA
        batch_params["cut_point_include_idx"] = batch_params["cut_point_include_idx"].apply(list)
        for idx,row in batch_params.iterrows():
            if row.fastq_r1 is None:
                raise Exception("At least one fastq file must be given as a command line parameter or be specified in the batch settings file with the heading 'fastq_r1' (fastq_r1 on row %s '%s' is invalid)"%(int(idx)+1,row.fastq_r1))
            CRISPRessoShared.check_file(row.fastq_r1)

            if row.fastq_r2 != "":
                CRISPRessoShared.check_file(row.fastq_r2)

            curr_amplicon_seq_str = row.amplicon_seq
            if curr_amplicon_seq_str is None:
                raise Exception("Amplicon sequence must be given as a command line parameter or be specified in the batch settings file with the heading 'amplicon_seq' (Amplicon seq on row %s '%s' is invalid)"%(int(idx)+1,curr_amplicon_seq_str))

            guides_are_in_amplicon = {} #dict of whether a guide is in at least one amplicon sequence
            #iterate through amplicons
            for curr_amplicon_seq in curr_amplicon_seq_str.split(','):
                this_include_idxs=[] #mask for bp to include for this amplicon seq, as specified by sgRNA cut points
                wrong_nt=CRISPRessoShared.find_wrong_nt(curr_amplicon_seq)
                if wrong_nt:
                    raise NTException('The amplicon sequence in row %d (%s) contains incorrect characters:%s' % (idx+1,curr_amplicon_seq_str,' '.join(wrong_nt)))

                #iterate through guides
                curr_guide_seq_string = row.guide_seq
                sgRNA_intervals = []
                if curr_guide_seq_string is not None and curr_guide_seq_string != "":
                    for curr_guide_seq in curr_guide_seq_string.strip().upper().split(','):
                        wrong_nt=CRISPRessoShared.find_wrong_nt(curr_guide_seq)
                        if wrong_nt:
                            raise NTException('The sgRNA sequence in row %d (%s) contains incorrect characters:%s'  % (idx+1,curr_guide_seq, ' '.join(wrong_nt)))

                        offset_fw=cleavage_offset+len(curr_guide_seq)-1
                        offset_rc=(-cleavage_offset)-1

                        cut_points = [m.start() + offset_fw for m in re.finditer(curr_guide_seq, curr_amplicon_seq)] + \
                                    [m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse_complement(curr_guide_seq), curr_amplicon_seq)] + \
                                    [m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse(curr_guide_seq), curr_amplicon_seq)]
                        sgRNA_intervals += [(m.start(),m.start()+len(curr_guide_seq)-1) for m in re.finditer(curr_guide_seq, curr_amplicon_seq)]+\
                                      [(m.start(),m.start()+len(curr_guide_seq)-1) for m in re.finditer(CRISPRessoShared.reverse_complement(curr_guide_seq), curr_amplicon_seq)]+\
                                      [(m.start(),m.start()+len(curr_guide_seq)-1) for m in re.finditer(CRISPRessoShared.reverse(curr_guide_seq), curr_amplicon_seq)]

                        #create mask of positions in which to include/exclude indels
                        if cut_points and offset_around_cut_to_plot >0:
                            for cut_p in cut_points:
                                st=max(0,cut_p-offset_around_cut_to_plot+1)
                                en=min(len(curr_amplicon_seq)-1,cut_p+offset_around_cut_to_plot+1)
                                this_include_idxs.extend(range(st,en))
                    if (not cut_points) and not (curr_guide_seq in guides_are_in_amplicon):
                        guides_are_in_amplicon[curr_guide_seq] = 0
                    elif (cut_points):
                        guides_are_in_amplicon[curr_guide_seq] = 1

                this_exclude_idxs=[]

                if row.exclude_bp_from_left:
                   this_exclude_idxs+=range(row.exclude_bp_from_left)

                if row.exclude_bp_from_right:
                   this_exclude_idxs+=range(len(curr_amplicon_seq))[-row.exclude_bp_from_right:]

                this_include_idxs=set(np.setdiff1d(this_include_idxs,this_exclude_idxs))
                this_include_idxs = sorted(list(set(this_include_idxs)))

                batch_params.ix[idx,"cut_point_include_idx"].append(this_include_idxs)
                batch_params.ix[idx,"sgRNA_intervals"].append(sgRNA_intervals)
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

        info("Running CRISPResso with %d processes" % args.n_processes)
        pool = mp.Pool(processes = args.n_processes)
        idxs = range(batch_count)
#        print("running on odir: %s\nopt: %s\nargs: %s\nbatchParams: %s\n"%(str(OUTPUT_DIRECTORY),str(crispresso_options_for_batch),str(args),str(batch_params)))
        pFunc = partial(runCrispresso,OUTPUT_DIRECTORY,crispresso_options_for_batch,batch_params)

        #handle signals -- bug in python 2.7 (https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python)
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        try:
            res = pool.map_async(pFunc,idxs)
            res.get(60*60) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            pool.terminate()
            info('Caught SIGINT. Program Terminated')
            raise Exception('CRISPResso2 Terminated')
            exit (0)
        except Exception as e:
            print('CAUGHT EXCEPTION HERE!!!')
            raise e
        else:
            info("Finished all batches")
            pool.close()
        pool.join()

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
                    folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)

                    nucleotide_frequency_file=os.path.join(folder_name,amplicon_name+'.nucleotide_frequency_table.txt')
                    ampSeq_nf,nuc_freqs = CRISPRessoShared.parse_count_file(nucleotide_frequency_file)

                    nucleotide_pct_file=os.path.join(folder_name,amplicon_name+'.nucleotide_percentage_table.txt')
                    ampSeq_np,nuc_pcts = CRISPRessoShared.parse_count_file(nucleotide_pct_file)

                    count_file=os.path.join(folder_name,amplicon_name+'.modification_count_vectors.txt')
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



                #if guides are all the same, merge substitutions and perform base editor comparison at guide target region
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
                            newend = len(include_idxs)
                        elif newstart == None:
                            newstart = 0
                        elif newend == None:
                            newend = len(include_idxs)
                        #and add it to the list
                        sub_sgRNA_intervals.append((newstart,newend))

                    CRISPRessoPlot.plot_nucleotide_quilt(sub_nucleotide_percentage_summary_df,sub_modification_percentage_summary_df,_jp(amplicon_name + '.TARGET_NUCLEOTIDE_PERCENTAGE_QUILT_INDEL'),args.save_also_png,sgRNA_intervals=sub_sgRNA_intervals)
                    if args.base_editor_mode:
                        CRISPRessoPlot.plot_conversion_map(sub_nucleotide_percentage_summary_df,_jp(amplicon_name + '.TARGET_NUCLEOTIDE_CONVERSION'),args.conversion_nuc_from,args.conversion_nuc_to,args.save_also_png,sgRNA_intervals=sub_sgRNA_intervals)

                    CRISPRessoPlot.plot_nucleotide_quilt(nucleotide_percentage_summary_df,modification_percentage_summary_df,_jp(amplicon_name + '.NUCLEOTIDE_PERCENTAGE_QUILT_INDEL'),args.save_also_png,sgRNA_intervals=sgRNA_intervals)
                    if args.base_editor_mode:
                        CRISPRessoPlot.plot_conversion_map(nucleotide_percentage_summary_df,_jp(amplicon_name + '.NUCLEOTIDE_CONVERSION'),args.conversion_nuc_from,args.conversion_nuc_to,args.save_also_png,sgRNA_intervals=sgRNA_intervals)
                else: #guides are not the same
                    CRISPRessoPlot.plot_nucleotide_quilt(nucleotide_percentage_summary_df,modification_percentage_summary_df,_jp(amplicon_name + '.NUCLEOTIDE_PERCENTAGE_QUILT_INDEL'),args.save_also_png)
                    if args.base_editor_mode:
                        CRISPRessoPlot.plot_conversion_map(nucleotide_percentage_summary_df,_jp(amplicon_name + '.NUCLEOTIDE_CONVERSION'),args.conversion_nuc_from,args.conversion_nuc_to,args.save_also_png)

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
                    folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)
                    amplicon_modification_file=os.path.join(folder_name,'CRISPResso_quantification_of_editing_frequency.txt')
                    with open(amplicon_modification_file,'r') as infile:
                        file_head = infile.readline()
                        if not wrote_header:
                            outfile.write('Batch\t' + file_head)
                            wrote_header = True
                        line = infile.readline()
                        while(line):
                            outfile.write(batchName + "\t" + line)
                            line = infile.readline()

            #summarize alignment
            with open(_jp('CRISPRessoBatch_mapping_statistics.txt'),'w') as outfile:
                wrote_header = False
                for idx,row in batch_params.iterrows():
                    batchName = row["name"]
                    folder_name = os.path.join(OUTPUT_DIRECTORY,'CRISPResso_on_%s' % batchName)
                    amplicon_modification_file=os.path.join(folder_name,'CRISPResso_mapping_statistics.txt')
                    with open(amplicon_modification_file,'r') as infile:
                        file_head = infile.readline()
                        if not wrote_header:
                            outfile.write('Batch\t' + file_head)
                            wrote_header = True
                        line = infile.readline()
                        while(line):
                            outfile.write(batchName + "\t" + line)
                            line = infile.readline()

        info('All Done!')
        print CRISPRessoShared.get_crispresso_footer()
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
