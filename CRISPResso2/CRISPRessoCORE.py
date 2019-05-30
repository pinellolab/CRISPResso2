#!/usr/bin/env python
# -*- coding: utf8 -*-

'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''



import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import argparse
from collections import defaultdict
from copy import deepcopy
import errno
import gzip
import zipfile
import os
import re
import subprocess as sb
import traceback
import unicodedata

if running_python3:
    import pickle as cp #python 3
else:
    import cPickle as cp #python 2.7

from CRISPResso2 import CRISPRessoCOREResources
from CRISPResso2 import CRISPRessoReport
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoPlot
from CRISPResso2 import CRISPResso2Align

from datetime import datetime
present = datetime.now()
#d1 = datetime.strptime('21/07/2019','%d/%m/%Y')
#if present > d1:
#    print('\nYour version of CRISPResso2 is out of date. Please download a new version.\n')
#    sys.exit(1)

import logging
#from test._mock_backport import inplace
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

####Support functions###

_ROOT = os.path.abspath(os.path.dirname(__file__))

def check_library(library_name):
    try:
        return __import__(library_name)
    except:
        error('You need to install %s module to use CRISPResso!' % library_name)
        sys.exit(1)

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def check_program(binary_name,download_url=None):
    if not which(binary_name):
        error('You need to install and have the command #####%s##### in your PATH variable to use CRISPResso!\n Please read the documentation!' % binary_name)
        if download_url:
            error('You can download it here:%s' % download_url)
        sys.exit(1)



def get_avg_read_length_fastq(fastq_filename):
     cmd=('z' if fastq_filename.endswith('.gz') else '' ) +('cat < \"%s\"' % fastq_filename)+\
                  r''' | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' '''
     p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
     return int(p.communicate()[0].strip())

def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < \"%s\" | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     return int(float(p.communicate()[0])/4.0)

#import time
#start = time.time()
matplotlib=check_library('matplotlib')
#end = time.time()
#start = time.time()
from matplotlib import font_manager as fm
CRISPRessoPlot.setMatplotlibDefaults()
#end = time.time()

#start = time.time()
plt=check_library('pylab')
#end = time.time()

from matplotlib import font_manager as fm
import matplotlib.gridspec as gridspec

pd=check_library('pandas')
np=check_library('numpy')
check_program('flash')

#start = time.time()
sns=check_library('seaborn')
#end = time.time()
sns.set_context('poster')
sns.set(font_scale=2.2)
sns.set_style('white')

#########################################


def process_fastq(fastq_filename,variantCache,ref_names,refs,args):
    """process_fastq processes each of the reads contained in a fastq file, given a cache of pre-computed variants
        fastqIn: name of fastq (e.g. output of FLASH)
            This file can be gzipped or plain text

        variantCache: dict with keys: sequence
            dict with keys:
                'count' : number of time sequence was observed
                'aln_ref_names' : names of reference it was aligned to
                'aln_scores' : score of alignment to each reference
                'class_name' : string with class names it was aligned to
                for each reference, there is a key: variant_ref_name with a payload object
            # payload object:
            The payload is a dict with keys:
                ### from CRISPRessoCOREResources.find_indels_substitutions
                # 'all_insertion_positions' #arr with 1's where there are insertions (including those outside of include_idxs quantification window)
                # 'all_insertion_left_positions' #arr with 1's to the left of where the insertion occurs
                # 'insertion_positions' # arr with 1's where there are insertions (1bp before and 1bp after insertion) that overlap with include_idxs quantification window
                # 'insertion_coordinates' # one entry per insertion, tuple of (start,end)
                # 'insertion_sizes'
                # 'insertion_n'
                # 'all_deletion_positions' #arr with 1's where there are insertions
                # 'deletion_positions' #arr with 1's where there are insertions that overlap the include_idxs quantification window
                # 'deletion_coordinates' # one entry per deletion
                # 'deletion_sizes' # correspond to entries in 'deletion_coordinates'
                # 'deletion_n'
                # 'all_substitution_positions'
                # 'substitution_positions'
                # 'substitution_n'
                # 'substitution_values'
                # 'ref_positions'
                ### added in this function
                # 'ref_name' # name of sequence that it most closely aligns to
                # 'classification' # MODIFIED or UNMODIFIED or AMBIGUOUS
                # 'aln_scores' # scores of alignment to each other reference sequence
                # 'aln_seq' # NW-aligned sequence
                # 'aln_ref' # NW-aligned sequence of corresponding reference (ref_name)

        refNameList: list of reference names
        refs: dictionary of sequences name>ref object
            ##ref object:
                # 'name'
                # 'sequence'
                # 'sequence_length'
                # 'min_aln_score' #sequence must align with at least this score
                # 'gap_incentive' #incentive for gaps at each position of the reference - to force gaps at the cut points, the indices of these cut points are set to 1  i.e. gap_incentive[4] = 1 would incentivise alignments with insertions to the right of the 4th character in the reference, or deletions of the 4th character in the reference.
                # 'contains_guide'
                # 'sgRNA_intervals'
                # 'sgRNA_sequences'
                # 'sgRNA_cut_points'
                # 'contains_coding_seq'
                # 'exon_positions'
                # 'exon_intervals'
                # 'exon_len_mods': the modification to the original exon length (if we copied the exon positions from another reference, this reference could introduce an indel, resulting in a non-zero length modification)
                # 'splicing_positions'
                # 'include_idxs' # sorted numpy array
                # 'exclude_idxs'
                # 'plot_idxs' #sorted numpy array
                # 'idx_cloned_from' #if this reference didn't contain a guide (or exon sequence), it was aligned to 'idx_cloned_from' reference, and cut_points, gap_incentive, sgRNA_intervals, inculde_idx, ane exon information were cloned from it (at the appropriate indices)
           Examples of these seqences can include:
           -the amplicon sequence
           -the repaired CRISPR expected output
           -allelic varaints if two variants are known to exist

        """

    N_TOT_READS = 0
    N_CACHED_ALN = 0 # read was found in cache
    N_CACHED_NOTALN = 0 #read was found in 'not aligned' cache
    N_COMPUTED_ALN = 0 # not in cache, aligned to at least 1 sequence with min cutoff
    N_COMPUTED_NOTALN = 0 #not in cache, not aligned to any sequence with min cutoff

    aln_matrix_loc = os.path.join(_ROOT,args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    if (args.needleman_wunsch_gap_open > 0):
        raise CRISPRessoShared.BadParameterException("Needleman Wunsch gap open penalty must be <= 0")
    if (args.needleman_wunsch_gap_extend > 0):
        raise CRISPRessoShared.BadParameterException("Needleman Wunsch gap extend penalty must be <= 0")


    not_aln = {} #cache for reads that don't align

    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename)
    else:
        fastq_handle=open(fastq_filename)

    count_seed_fw = 0
    count_seed_rv = 0
    count_seed_both = 0
    while(fastq_handle.readline()):

        #read through fastq in sets of 4
        fastq_seq = fastq_handle.readline().strip()
        fastq_strand = fastq_handle.readline()
        fastq_qual = fastq_handle.readline()

        if (N_TOT_READS % 10000 == 0):
            info("Processing reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS,N_COMPUTED_ALN,N_CACHED_ALN,N_COMPUTED_NOTALN,N_CACHED_NOTALN))

        N_TOT_READS+=1

        payload = []

        #if the sequence has been seen and can't be aligned, skip it
        if (fastq_seq in not_aln):
            N_CACHED_NOTALN += 1
            continue
        #if the sequence is already associated with a variant in the variant cache, pull it out
        if (fastq_seq in variantCache):
            N_CACHED_ALN+=1
            variantCache[fastq_seq]['count'] += 1

        #otherwise, create a new variant object, and put it in the cache
        else:
            aln_scores = []
            best_match_score = -1
            best_match_s1s = []
            best_match_s2s = []
            best_match_names = []
            for idx,ref_name in enumerate(ref_names):
                #get alignment and score from cython
                #score = 100 * #matchedBases / length(including gaps)
                seed_i = 0
                found_forward_count = 0
                found_reverse_count = 0
                while seed_i < args.aln_seed_count and seed_i < len(refs[ref_name]['fw_seeds']):
                    if refs[ref_name]['fw_seeds'][seed_i] in fastq_seq: #is forward
                        found_forward_count += 1
                    if refs[ref_name]['rc_seeds'][seed_i] in fastq_seq: #is rc
                        found_reverse_count += 1
                    seed_i += 1
                if found_forward_count > args.aln_seed_min and found_reverse_count == 0:
                    fws1,fws2,fwscore=CRISPResso2Align.global_align(fastq_seq, refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                    s1 = fws1
                    s2 = fws2
                    score = fwscore
                    count_seed_fw += 1
                elif found_forward_count == 0 and found_reverse_count > args.aln_seed_min:
                    rvs1,rvs2,rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                    s1 = rvs1
                    s2 = rvs2
                    score = rvscore
                    count_seed_rv += 1
                else:
                    count_seed_both += 1
                    fws1,fws2,fwscore=CRISPResso2Align.global_align(fastq_seq, refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                    rvs1,rvs2,rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                    s1 = fws1
                    s2 = fws2
                    score = fwscore
                    if (rvscore > fwscore):
                        s1 = rvs1
                        s2 = rvs2
                        score = rvscore

#                print "for " + ref_name + " got fws1: " + str(fws1) + " and fws2: " + str(fws2) + " score: " +str(fwscore)
                aln_scores.append(score)

                #reads are matched to the reference to which they best align. The 'min_aln_score' is calculated using only the changes in 'include_idxs'
                if score > best_match_score and score > refs[ref_name]['min_aln_score']:
                    best_match_score = score
                    best_match_s1s = [s1]
                    best_match_s2s = [s2]
                    best_match_names = [ref_name]
                elif score == best_match_score:
                    best_match_s1s.append(s1)
                    best_match_s2s.append(s2)
                    best_match_names.append(ref_name)

#            print "bestMatchscore: " + str(best_match_score) + " to " + str(best_match_names)

            if best_match_score > 0:
                N_COMPUTED_ALN+=1
                variantCache[fastq_seq] = {}
                variantCache[fastq_seq]['count'] = 1
                variantCache[fastq_seq]['aln_ref_names'] = best_match_names
                variantCache[fastq_seq]['aln_scores'] = aln_scores
                class_names = []

#                for idx, best_match_name in enumerate(best_match_names):
                for idx in range(len(best_match_names)):
                    best_match_name = best_match_names[idx]
                    payload=CRISPRessoCOREResources.find_indels_substitutions(best_match_s1s[idx],best_match_s2s[idx],refs[best_match_name]['include_idxs'])
                    payload['ref_name'] = best_match_name
                    payload['aln_scores'] = aln_scores

                    # If there is an insertion/deletion/substitution in the quantification window, the read is modified.
                    is_modified = False
                    if not args.ignore_deletions and payload['deletion_n'] > 0:
                        is_modified = True
                    elif not args.ignore_insertions and payload['insertion_n'] > 0:
                        is_modified = True
                    elif not args.ignore_substitutions and payload['substitution_n'] > 0:
                        is_modified = True

                    if is_modified:
                        class_names.append(best_match_name+"_MODIFIED")
                        payload['classification'] = 'MODIFIED'
                    else:
                        class_names.append(best_match_name+"_UNMODIFIED")
                        payload['classification'] = 'UNMODIFIED'

                    payload['aln_seq'] = best_match_s1s[idx]
                    payload['aln_ref'] = best_match_s2s[idx]

                    variantCache[fastq_seq]['variant_'+best_match_name] = payload

                variantCache[fastq_seq]['class_name'] = "&".join(class_names)

            else:
                N_COMPUTED_NOTALN+=1
                not_aln[fastq_seq] = 1


    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS,N_COMPUTED_ALN,N_CACHED_ALN,N_COMPUTED_NOTALN,N_CACHED_NOTALN))
    aln_stats = {"N_TOT_READS" : N_TOT_READS,
               "N_CACHED_ALN" : N_CACHED_ALN,
               "N_CACHED_NOTALN" : N_CACHED_NOTALN,
               "N_COMPUTED_ALN" : N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN" : N_COMPUTED_NOTALN,
               'count_seed_fw' : count_seed_fw,
               'count_seed_rv' : count_seed_rv,
               'count_seed_both' : count_seed_both
               }
    return(aln_stats)


def add_hist(hist_to_add,hist_global):
    for key,value in hist_to_add.iteritems():
        hist_global[key]+=value
    return hist_global



def split_paired_end_reads_single_file(fastq_filename,output_filename_r1,output_filename_r2):

    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename)
    else:
        fastq_handle=open(fastq_filename)

    #we cannot use with on gzip with python 2.6 :(
    try:
        fastq_splitted_outfile_r1=gzip.open(output_filename_r1,'w+')
        fastq_splitted_outfile_r2=gzip.open(output_filename_r2,'w+')
        [fastq_splitted_outfile_r1.write(line) if (i % 8 < 4) else fastq_splitted_outfile_r2.write(line) for i, line in enumerate(fastq_handle)]
    except:
        raise CRISPRessoShared.BadParameterException('Error in splitting read pairs from a single file')

    return output_filename_r1,output_filename_r2


def main():

    def print_stacktrace_if_debug():
        debug_flag = False
        if 'args' in vars() and 'debug' in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)
            error(traceback.format_exc())

    try:

        start_time =  datetime.now()
        start_time_string =  start_time.strftime('%Y-%m-%d %H:%M:%S')
        description = ['~~~CRISPResso 2~~~','-Analysis of genome editing outcomes from deep sequencing data-']
        header = CRISPRessoShared.get_crispresso_header(description=description,header_str=None)
        print(header)

        args = CRISPRessoShared.getCRISPRessoArgParser(requiredParams={'fastq_r1':True}).parse_args()

        aln_matrix_loc = os.path.join(_ROOT,"EDNAFULL")
        CRISPRessoShared.check_file(aln_matrix_loc)
        aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

        #check files
        CRISPRessoShared.check_file(args.fastq_r1)
        if args.fastq_r2:
            CRISPRessoShared.check_file(args.fastq_r2)


        if args.amplicon_seq is None and args.auto is False:
            raise CRISPRessoShared.BadParameterException('Please provide an amplicon sequence for analysis.')

        #create output directory
        get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')

        #normalize name and remove not allowed characters
        if not args.name:
            if args.fastq_r2!='':
                database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1),get_name_from_fasta(args.fastq_r2))
            else:
                database_id='%s' % get_name_from_fasta(args.fastq_r1)

        else:
            clean_name=CRISPRessoShared.slugify(args.name)
            if args.name!= clean_name:
                warn('The specified name %s contained invalid characters and was changed to: %s' % (args.name,clean_name))
            database_id=clean_name


        clean_file_prefix = ""
        if args.file_prefix != "":
            clean_file_prefix = CRISPRessoShared.slugify(args.file_prefix)
            if not clean_file_prefix.endswith("."):
                clean_file_prefix += "."

        OUTPUT_DIRECTORY='CRISPResso_on_%s' % database_id

        if args.output_folder:
            OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,clean_file_prefix + filename) #handy function to put a file in the output directory

        crispresso2_info_file = os.path.join(OUTPUT_DIRECTORY,'CRISPResso2_info.pickle')
        crispresso2_info = {} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['version'] = CRISPRessoShared.__version__
        crispresso2_info['args'] = deepcopy(args)

        log_filename=_jp('CRISPResso_RUNNING_LOG.txt')
        crispresso2_info['log_filename'] = os.path.basename(log_filename)

        crispresso2_info['name'] = database_id

        if args.no_rerun:
            if os.path.exists(crispresso2_info_file):
                previous_run_data = cp.load(open(crispresso2_info_file,'rb'))
                if previous_run_data['version'] == CRISPRessoShared.__version__:
                    args_are_same = True
                    for arg in vars(args):
                        if arg is "no_rerun":
                            continue
                        if arg not in vars(previous_run_data['args']):
                            info('Comparing current run to previous run: old run had argument ' + str(arg) + ' \nRerunning.')
                            args_are_same = False
                        elif str(getattr(previous_run_data['args'],arg)) != str(getattr(args,arg)):
                            info('Comparing current run to previous run:\n\told argument ' + str(arg) + ' = ' + str(getattr(previous_run_data['args'],arg)) + '\n\tnew argument: ' + str(arg) + ' = ' + str(getattr(args,arg)) + '\nRerunning.')
                            args_are_same = False

                    if args_are_same:
                        info('Analysis already completed on %s!'%previous_run_data['end_time_string'])
                        sys.exit(0)
                else:
                    info('The no_rerun flag is set, but this analysis will be rerun because the existing run was performed using an old version of CRISPResso (' + str(previous_run_data['version']) + ').')

        #### ASSERT GUIDE(S)
        guides = []
        if args.guide_seq:
            for current_guide_seq in args.guide_seq.split(','):
                wrong_nt=CRISPRessoShared.find_wrong_nt(current_guide_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The sgRNA sequence contains bad characters:%s'  % ' '.join(wrong_nt))
                guides.append(current_guide_seq)

        ###FRAMESHIFT SUPPORT###
        coding_seqs = []
        if args.coding_seq:
            for exon_seq in args.coding_seq.strip().upper().split(','):
                #check for wrong NT
                wrong_nt=CRISPRessoShared.find_wrong_nt(exon_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The coding sequence contains bad characters:%s' % ' '.join(wrong_nt))

                coding_seqs.append(exon_seq)

        ####SET REFERENCES TO COMPARE###
        ref_names = [] #ordered list of names
        refs = {} #dict of ref_name > ref object


#        #if we should automatically infer amplicon sequence, pull out the most frequent read and assign it to be the amplicon
        if args.auto:
            number_of_reads_to_consider = 1000 * 4 #1000 fastq sequences (4 lines each)

            amplicon_seq_arr = CRISPRessoShared.guess_amplicons(args.fastq_r1,args.fastq_r2,number_of_reads_to_consider,args.flash_command,args.max_paired_end_reads_overlap,args.min_paired_end_reads_overlap,
                aln_matrix,args.needleman_wunsch_gap_open,args.needleman_wunsch_gap_extend)
            amp_dummy = ['']
            amp_dummy.extend(list(range(2,len(amplicon_seq_arr)+1)))
            amplicon_name_arr = ['Amplicon'+str(x) for x in amp_dummy]
            if len(guides) == 0:
                for amplicon_seq in amplicon_seq_arr:
                    (potential_guide,is_base_editor) = CRISPRessoShared.guess_guides(amplicon_seq,args.fastq_r1,args.fastq_r2,number_of_reads_to_consider,args.flash_command,
                        args.max_paired_end_reads_overlap,args.min_paired_end_reads_overlap,aln_matrix,args.needleman_wunsch_gap_open,args.needleman_wunsch_gap_extend)
                    if potential_guide is not None and potential_guide not in guides:
                        guides.append(potential_guide)

            amplicon_min_alignment_score_arr = []
            plural_string = ""
            if len(amplicon_seq_arr) > 1:
                plural_string = "s"
            info("Auto-detected %d reference amplicon%s"%(len(amplicon_seq_arr),plural_string))

            if args.debug:
                for idx,seq in enumerate(amplicon_seq_arr):
                    info('Detected amplicon ' + str(idx) + ":" + str(seq))

            if len(guides) > 1:
                plural_string = "s"
            info("Auto-detected %d guide%s"%(len(guides),plural_string))
            if args.debug:
                for idx,seq in enumerate(guides):
                    info('Detected guide ' + str(idx) + ":" + str(seq))
            if amplicon_seq_arr == 0:
                raise BadParameterException("Cannot automatically infer amplicon sequence.")

        else: #not auto
            amplicon_seq_arr = args.amplicon_seq.split(",")
            amplicon_name_arr = args.amplicon_name.split(",")
            #split on commas, only accept empty values
            amplicon_min_alignment_score_arr = [float(x) for x in args.amplicon_min_alignment_score.split(",") if x]

        if args.expected_hdr_amplicon_seq != "":
            amplicon_seq_arr.append(args.expected_hdr_amplicon_seq)
            amplicon_name_arr.append('HDR')

        found_guide_seq = [False]*len(guides)
        found_coding_seq = [False]*len(coding_seqs)

        max_amplicon_len = 0 #for flash
        min_amplicon_len = 99**99 #for flash

        quant_window_coordinates_arr = []
        if args.quantification_window_coordinates is not None:
            quant_window_coordinates_arr = args.quantification_window_coordinates.split(",")

        for idx,seq in enumerate(amplicon_seq_arr):
            this_seq = seq.strip().upper()
            this_seq_length = len(this_seq)
            if this_seq_length > max_amplicon_len:
                max_amplicon_len = this_seq_length
            if this_seq_length < min_amplicon_len:
                min_amplicon_len = this_seq_length

            this_name = 'Amplicon'+str(idx)
            if idx < len(amplicon_name_arr):
                this_name = amplicon_name_arr[idx]

            wrong_nt=CRISPRessoShared.find_wrong_nt(this_seq)
            if wrong_nt:
                raise CRISPRessoShared.NTException('Reference amplicon sequence %d (%s) contains invalid characters:%s' % idx,this_name, ' '.join(wrong_nt))

            this_min_aln_score = args.default_min_aln_score
            if idx < len(amplicon_min_alignment_score_arr):
                this_min_aln_score = amplicon_min_alignment_score_arr[idx]

            this_quant_window_coordinates = None
            if idx < len(quant_window_coordinates_arr) and quant_window_coordinates_arr[idx] != "":
                this_quant_window_coordinates = quant_window_coordinates_arr[idx]


            # Calculate cut sites for this reference
            (this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_idxs, this_include_idxs,
                this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(this_seq,guides,args.quantification_window_center,
                args.quantification_window_size,this_quant_window_coordinates,args.exclude_bp_from_left,args.exclude_bp_from_right,args.plot_window_size)

            this_contains_guide = False
            if len(this_sgRNA_sequences) > 0:
                this_contains_guide = True

            for guide_idx, guide_seq in enumerate(guides):
                if guide_seq in this_sgRNA_sequences:
                    found_guide_seq[guide_idx] = True

            # Calculate coding sequence for this reference
            this_exon_positions = set()
            this_exon_intervals = []
            this_exon_len_mods = []
            this_splicing_positions = []
            this_contains_coding_seq = False
            for exon_idx, exon_seq in enumerate(coding_seqs):
                st_exon = this_seq.find(exon_seq)
                if st_exon >= 0:
                    found_coding_seq[exon_idx] = True
                    this_contains_coding_seq = True
                    en_exon = st_exon + len(exon_seq)  # this do not include the upper bound as usual in python
                    this_exon_len_mods.append(0)
                    this_exon_intervals.append((st_exon, en_exon))
                    this_exon_positions = this_exon_positions.union(set(range(st_exon, en_exon)))

                    # consider 2 base pairs before and after each exon
                    this_splicing_positions += [max(0, st_exon - 2), max(0, st_exon - 1), min(this_seq_length - 1, en_exon), min(this_seq_length - 1, en_exon + 1)]

            this_exon_positions = sorted(this_exon_positions)

            # protect from the wrong splitting of exons by the users to avoid false splicing sites
            this_splicing_positions = set(this_splicing_positions).difference(this_exon_positions)

            this_gap_incentive = np.zeros(this_seq_length+1,dtype=np.int)
            for cut_point in this_sgRNA_cut_points:
                this_gap_incentive[cut_point+1] = args.needleman_wunsch_gap_incentive

            seq_rc = CRISPRessoShared.reverse_complement(this_seq)
            seeds = []
            rc_seeds = []
            seedStarts = list(range(args.exclude_bp_from_left,this_seq_length-args.exclude_bp_from_right-args.aln_seed_len,args.aln_seed_count)) #define all possible seed starts
            for seedStart in seedStarts:
                attemptsToFindSeed = 0
                thisSeedStart = seedStart
                potentialSeed = this_seq[thisSeedStart:thisSeedStart+args.aln_seed_len]
                #seed shouldn't be in reverse complement of sequence, and it should also not be the same as another seed
                while potentialSeed in seq_rc or potentialSeed in seeds:
                    attemptsToFindSeed += 1
                    if attemptsToFindSeed > 100:
                        #raise CRISPRessoShared.AlignmentException("Can't find alignment seed that is unique to the forward sequence")
                        break
                    # if this seed would extend past the end of the sequence, reset to position 0 (possibly in the excluded region, but hey, we're desperate)
                    if (seedStart > this_seq_length - args.aln_seed_len):
                        thisSeedStart = 0
                    thisSeedStart += 1
                    potentialSeed = this_seq[thisSeedStart:thisSeedStart+args.aln_seed_len]
                seed_rc = CRISPRessoShared.reverse_complement(potentialSeed)
                if seed_rc in this_seq:
                    #raise CRISPRessoShared.AlignmentException("Reverse compliment of seed %s is in amplicon %s even though seed is not in reverse compliment of amplicon"%(seed,amplicon))
                    continue
                if potentialSeed not in seq_rc:
                    seeds.append(potentialSeed)
                    rc_seeds.append(seed_rc)


            refObj = {'name':this_name,
                   'sequence':this_seq,
                   'sequence_length':this_seq_length,
                   'min_aln_score':this_min_aln_score,
                   'gap_incentive':this_gap_incentive,
                   'sgRNA_cut_points':this_sgRNA_cut_points,
                   'sgRNA_intervals':this_sgRNA_intervals,
                   'sgRNA_sequences':this_sgRNA_sequences,
                   'sgRNA_plot_idxs':this_sgRNA_plot_idxs,
                   'contains_guide':this_contains_guide,
                   'contains_coding_seq':this_contains_coding_seq,
                   'exon_positions':this_exon_positions,
                   'exon_len_mods':this_exon_len_mods,
                   'exon_intervals':this_exon_intervals,
                   'splicing_positions':this_splicing_positions,
                   'include_idxs':this_include_idxs,
                   'exclude_idxs':this_exclude_idxs,
                   'idx_cloned_from':None,
                   'fw_seeds':seeds,
                   'rc_seeds':rc_seeds,
                   }
            ref_names.append(this_name)
            refs[this_name] = refObj

        #throw error if guides, or coding seqs not found in any reference
        if args.guide_seq:
            for idx, presence_bool in enumerate(found_guide_seq):
                if not presence_bool:
                    raise CRISPRessoShared.SgRNASequenceException('The guide sequence %d (%s) provided is not present in the amplicon sequences!\n\nPlease check your input!' % (idx, guides[idx]))

        if args.coding_seq:
            for idx, presence_bool in enumerate(found_coding_seq):
                if not presence_bool:
                    raise CRISPRessoShared.ExonSequenceException('The coding subsequence %d (%s) provided is not contained in any amplicon sequence!\n\nPlease check your input!' % (idx,coding_seqs[idx]))


        #clone cut points and include idx from first reference where those are set (also exons)
        clone_ref_name = None
        clone_has_cut_points = False
        clone_has_exons = False
        for ref_name in ref_names:
            cut_points = refs[ref_name]['sgRNA_cut_points']
            exon_positions = refs[ref_name]['exon_positions']
            if cut_points:
                if len(ref_names) > 1 and args.debug:
                    info("Using cut points from %s as template for other references"%ref_name)
                clone_ref_name = ref_name
                clone_has_cut_points = True
                if exon_positions:
                    clone_has_exons = True
                break
            if exon_positions:
                if len(ref_names) > 1 and args.debug:
                    info("Using exons positions from %s as template for other references"%ref_name)
                clone_ref_name = ref_name
                clone_has_exon_positions = True
                if cut_points:
                    clone_has_cut_points = True
                break

        if clone_ref_name is not None:
            for ref_name in ref_names:
                cut_points = refs[ref_name]['sgRNA_cut_points']
                sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
                exon_positions = refs[ref_name]['exon_positions']

                needs_cut_points = False
                needs_sgRNA_intervals = False
                needs_exon_positions = False

                if cut_points:
                    if len(ref_names) > 1 and args.debug:
                        info("Reference '%s' has cut points defined: %s. Not inferring."%(ref_name,cut_points))
                else:
                    needs_cut_points = True
                if sgRNA_intervals:
                    if len(ref_names) > 1 and args.debug:
                        info("Reference '%s' has sgRNA_intervals defined: %s. Not inferring."%(ref_name,sgRNA_intervals))
                else:
                    needs_sgRNA_intervals = True

                if exon_positions:
                    if len(exon_positions) > 1 and args.debug:
                        info("Reference '%s' has exon_positions defined: %s. Not inferring."%(ref_name,exon_positions))
                else:
                    needs_exon_positions = True

                if not needs_cut_points and not needs_sgRNA_intervals and not needs_exon_positions:
                    continue

                fws1,fws2,fwscore=CRISPResso2Align.global_align(refs[ref_name]['sequence'], refs[clone_ref_name]['sequence'],matrix=aln_matrix,gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,gap_incentive=refs[clone_ref_name]['gap_incentive'])
                if fwscore < 60:
                    continue

                if (needs_sgRNA_intervals or needs_cut_points) and clone_has_cut_points and args.debug:
                    info("Reference '%s' has NO cut points or sgRNA intervals idxs defined. Inferring from '%s'."%(ref_name,clone_ref_name))
                if needs_exon_positions and clone_has_exons and args.debug:
                    info("Reference '%s' has NO exon_positions defined. Inferring from '%s'."%(ref_name,clone_ref_name))
                #Create a list such that the nucleotide at ix in the old reference corresponds to s1inds[ix]
                s1inds = []
                s1ix = -1
                s2ix = -1
                for ix in range(len(fws1)):
                    if fws1[ix] != "-":
                        s1ix += 1
                    if fws2[ix] != "-":
                        s2ix += 1
                        s1inds.append(s1ix)
#                print("aln:\n%s\n%s"%(fws1,fws2))
#                print(str(s1inds))

                if (needs_cut_points or needs_sgRNA_intervals) and clone_has_cut_points:
                    this_cut_points = [s1inds[X] for X in refs[clone_ref_name]['sgRNA_cut_points']]
                    this_gap_incentive = np.zeros(refs[ref_name]['sequence_length']+1,dtype=np.int)
                    for cut_point in this_cut_points:
                        this_gap_incentive[cut_point + 1] = args.needleman_wunsch_gap_incentive

                    this_sgRNA_intervals = []
                    for (sgRNA_interval_start,sgRNA_interval_end) in refs[clone_ref_name]['sgRNA_intervals']:
                        this_sgRNA_intervals.append((s1inds[sgRNA_interval_start],s1inds[sgRNA_interval_end]))

                    this_sgRNA_plot_idxs = []
                    for plot_idx_list in refs[clone_ref_name]['sgRNA_plot_idxs']:
                        this_sgRNA_plot_idxs.append([s1inds[x] for x in plot_idx_list])

                    this_include_idxs = [s1inds[x] for x in refs[clone_ref_name]['include_idxs']]
                    #subtract any indices in 'exclude_idxs' -- e.g. in case some of the cloned include_idxs were near the read ends (excluded)
                    this_exclude_idxs = sorted(list(set(refs[ref_name]['exclude_idxs'])))
                    this_include_idxs = sorted(list(set(np.setdiff1d(this_include_idxs,this_exclude_idxs))))

                    refs[ref_name]['gap_incentive'] = this_gap_incentive
                    refs[ref_name]['sgRNA_cut_points'] = this_cut_points
                    refs[ref_name]['sgRNA_intervals'] = this_sgRNA_intervals
                    refs[ref_name]['sgRNA_sequences'] = refs[clone_ref_name]['sgRNA_sequences']
                    refs[ref_name]['sgRNA_plot_idxs'] = this_sgRNA_plot_idxs
                    refs[ref_name]['include_idxs'] = this_include_idxs
                    refs[ref_name]['contains_guide'] = True


                if needs_exon_positions and clone_has_exons:
                    this_exon_positions = set()
                    this_splicing_positions = []
                    this_exon_intervals = []
                    this_exon_len_mods = []
                    this_seq_length = refs[ref_name]['sequence_length']
                    for (exon_interval_start,exon_interval_end) in refs[clone_ref_name]['exon_intervals']:
                        this_exon_start = s1inds[exon_interval_start]
                        this_exon_end = s1inds[exon_interval_end]
                        this_exon_intervals.append((this_exon_start,this_exon_end))
                        this_exon_len_mods.append(((this_exon_end - this_exon_start)-(exon_interval_end - exon_interval_start)))
                        this_exon_positions = this_exon_positions.union(set(range(this_exon_start,this_exon_end)))
                        # consider 2 base pairs before and after each exon
                        this_splicing_positions += [max(0, this_exon_start - 2), max(0, this_exon_start - 1), min(this_seq_length - 1, this_exon_end), min(this_seq_length - 1, this_exon_end + 1)]

                    this_exon_positions = sorted(this_exon_positions)

                    # protect from the wrong splitting of exons by the users to avoid false splicing sites
                    this_splicing_positions = set(this_splicing_positions).difference(this_exon_positions)

                    refs[ref_name]['exon_positions'] = this_exon_positions
                    refs[ref_name]['exon_intervals'] = this_exon_intervals
                    refs[ref_name]['exon_len_mods'] = this_exon_len_mods
                    refs[ref_name]['splicing_positions'] = this_splicing_positions
                    refs[ref_name]['contains_coding_seq'] = True

                refs[ref_name]['idx_cloned_from'] = clone_ref_name


        crispresso_cmd_to_write = ' '.join(sys.argv)
        if args.write_cleaned_report:
            cmd_copy = sys.argv[:]
            cmd_copy[0] = 'CRISPResso'
            for i in range(len(cmd_copy)):
                if os.sep in cmd_copy[i]:
                    cmd_copy[i] = os.path.basename(cmd_copy[i])

            crispresso_cmd_to_write = ' '.join(cmd_copy) #clean command doesn't show the absolute path to the executable or other files
        crispresso2_info['command_used'] = crispresso_cmd_to_write


        try:
            os.makedirs(OUTPUT_DIRECTORY)
            info('Creating Folder %s' % OUTPUT_DIRECTORY)
#            info('Done!') #crispresso2 doesn't announce that the folder is created... save some electricity here
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        finally:
            logging.getLogger().addHandler(logging.FileHandler(log_filename))

            with open(log_filename,'w+') as outfile:
                outfile.write('CRISPResso version %s\n[Command used]:\n%s\n\n[Execution log]:\n' %(CRISPRessoShared.__version__,crispresso_cmd_to_write))

        N_READS_INPUT=get_n_reads_fastq(args.fastq_r1)


        args_string_arr = []
        clean_args_string_arr = []
        for arg in vars(args):
            val = str(getattr(args,arg))
            args_string_arr.append("%s: %s"%(str(arg),val))
            if os.sep in val:
                val = os.path.basename(val)
            clean_args_string_arr.append("%s: %s"%(str(arg),val))

        if args.write_cleaned_report:
            crispresso2_info['args_string'] = '\n'.join(sorted(clean_args_string_arr))
        else:
            crispresso2_info['args_string'] = '\n'.join(sorted(args_string_arr))

        crispresso2_info['start_time'] = start_time
        crispresso2_info['start_time_string'] = start_time_string

        if args.split_paired_end:
            if args.fastq_r2!='':
                raise CRISPRessoShared.BadParameterException('The option --split_paired_end is available only when a single fastq file is specified!')
            else:
                info('Splitting paired end single fastq file into two files...')
                args.fastq_r1,args.fastq_r2=split_paired_end_reads_single_file(args.fastq_r1,
                    output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_splitted_r1.fastq.gz'),
                    output_filename_r2=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_splitted_r2.fastq.gz'),)
                splitted_files_to_remove=[args.fastq_r1,args.fastq_r2]

                info('Done!')

        if args.min_average_read_quality>0 or args.min_single_bp_quality>0 or args.min_bp_quality_or_N>0:
            info('Filtering reads with average bp quality < %d and single bp quality < %d and replacing bases with quality < %d with N ...' % (args.min_average_read_quality,args.min_single_bp_quality,args.min_bp_quality_or_N))
            min_av_quality = None
            if args.min_average_read_quality > 0:
                min_av_quality = args.min_average_read_quality

            min_single_bp_quality = None
            if args.min_single_bp_quality > 0:
                min_single_bp_quality = args.min_single_bp_quality

            min_bp_quality_or_N = None
            if args.min_bp_quality_or_N > 0:
                min_bp_quality_or_N = args.min_bp_quality_or_N

            if args.fastq_r2!='':
                output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_filtered.fastq.gz')
                output_filename_r2=_jp(os.path.basename(args.fastq_r2.replace('.fastq','')).replace('.gz','')+'_filtered.fastq.gz')

                import filterFastqs
                filterFastqs.filterFastqs(fastq_r1=args.fastq_r1,fastq_r2=args.fastq_r2,fastq_r1_out=output_filename_r1,fastq_r2_out=output_filename_r2,min_bp_qual_in_read=min_single_bp_quality,min_av_read_qual=min_av_quality,min_bp_qual_or_N=min_bp_quality_or_N)

                args.fastq_r1 = output_filename_r1
                args.fastq_r2 = output_filename_r2

            else:
                output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq','')).replace('.gz','')+'_filtered.fastq.gz')

                import filterFastqs
                filterFastqs.filterFastqs(fastq_r1=args.fastq_r1,fastq_r1_out=output_filename_r1,min_bp_qual_in_read=min_single_bp_quality,min_av_read_qual=min_av_quality,min_bp_qual_or_N=min_bp_quality_or_N)

                args.fastq_r1 = output_filename_r1


        if args.fastq_r2=='': #single end reads
            #check if we need to trim
            if not args.trim_sequences:
                #create a symbolic link
                symlink_filename=_jp(os.path.basename(args.fastq_r1))
                CRISPRessoShared.force_symlink(os.path.abspath(args.fastq_r1),symlink_filename)
                output_forward_filename=symlink_filename
            else:
                output_forward_filename=_jp('reads.trimmed.fq.gz')
                #Trimming with trimmomatic
                cmd='%s SE -phred33 %s  %s %s >>%s 2>&1'\
                % (args.trimmomatic_command,args.fastq_r1,
                   output_forward_filename,
                   args.trimmomatic_options_string.replace('NexteraPE-PE.fa','TruSeq3-SE.fa'),
                   log_filename)
                #print cmd
                TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)

                if TRIMMOMATIC_STATUS:
                        raise CRISPRessoShared.TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')
                crispresso2_info['trimmomatic_command'] = cmd


            processed_output_filename=output_forward_filename

        else:#paired end reads case

            if not args.trim_sequences:
                output_forward_paired_filename=args.fastq_r1
                output_reverse_paired_filename=args.fastq_r2
            else:
                info('Trimming sequences with Trimmomatic...')
                output_forward_paired_filename=_jp('output_forward_paired.fq.gz')
                output_forward_unpaired_filename=_jp('output_forward_unpaired.fq.gz')
                output_reverse_paired_filename=_jp('output_reverse_paired.fq.gz')
                output_reverse_unpaired_filename=_jp('output_reverse_unpaired.fq.gz')

                #Trimming with trimmomatic
                cmd='%s PE -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
                    % (args.trimmomatic_command,
                        args.fastq_r1,args.fastq_r2,output_forward_paired_filename,
                        output_forward_unpaired_filename,output_reverse_paired_filename,
                        output_reverse_unpaired_filename,args.trimmomatic_options_string,log_filename)
                #print cmd
                TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)
                if TRIMMOMATIC_STATUS:
                    raise CRISPRessoShared.TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')
                crispresso2_info['trimmomatic_command'] = cmd

                info('Done!')


            info('Estimating average read length...')
            if args.debug:
                info('Checking average read length from ' + output_forward_paired_filename)
            if get_n_reads_fastq(output_forward_paired_filename):
                avg_read_length=get_avg_read_length_fastq(output_forward_paired_filename)
                if args.debug:
                    info('Average read length is ' + str(avg_read_length) + ' from ' + output_forward_paired_filename)
            else:
               raise CRISPRessoShared.NoReadsAfterQualityFilteringException('No reads survived the average or single bp quality filtering.')

            #Merging with Flash
            info('Merging paired sequences with Flash...')
            min_overlap = args.min_paired_end_reads_overlap
            max_overlap = args.max_paired_end_reads_overlap
            if args.stringent_flash_merging:
                expected_max_overlap=2*avg_read_length - min_amplicon_len
                expected_min_overlap=2*avg_read_length - max_amplicon_len
    #            print('avg read len: ' + str(avg_read_length))
    #            print('expected_max_overlap' + str(expected_max_overlap))
    #            print('expected_min_overlap' + str(expected_min_overlap))
    #            print('min amplicon len:' + str(min_amplicon_len))
    #            print('max amplicon len:' + str(max_amplicon_len))
                indel_overlap_tolerance = 10 # magic number bound on how many bp inserted/deleted in ~90% of reads (for flash)
                #max overlap is either the entire read (avg_read_length) or the expected amplicon length + indel tolerance
                max_overlap = max(10,min(avg_read_length, expected_max_overlap+indel_overlap_tolerance))
                #min overlap is either 4bp (as in crispresso1) or the expected amplicon length - indel tolerance
                min_overlap = max(4,expected_min_overlap-indel_overlap_tolerance)
    #            print('max_overlap: ' + str(max_overlap))
    #            print('min_overlap: ' + str(min_overlap))
                # if reads are longer than the amplicon, there is no way to tell flash to have them overlap like this..
                if avg_read_length > min_amplicon_len:
                    info('Warning: Reads are longer than amplicon.')
                    min_overlap = avg_read_length-10
                    max_overlap = 2*avg_read_length

            output_prefix = "out"
            if clean_file_prefix != "":
                output_prefix = clean_file_prefix + "out"
            cmd='%s %s %s --min-overlap %d --max-overlap %d --allow-outies -z -d %s -o %s >>%s 2>&1' %\
            (args.flash_command,
                 output_forward_paired_filename,
                 output_reverse_paired_filename,
                 min_overlap,
                 max_overlap,
                 OUTPUT_DIRECTORY,
                 output_prefix,
                 log_filename)
            #cmd='flash %s %s --min-overlap %d --max-overlap %d -f %d -z -d %s >>%s 2>&1' %\
            #(output_forward_paired_filename,
            #     output_reverse_paired_filename,
            #     args.min_paired_end_reads_overlap,
            #     max_amplicon_len,
            #     avg_read_length,
            #     OUTPUT_DIRECTORY,log_filename)
#            cmd='flash %s %s --min-overlap %d -f %d -r %d -s %d  -z -d %s >>%s 2>&1' %\
#            (output_forward_paired_filename,
#                 output_reverse_paired_filename,
#                 args.min_paired_end_reads_overlap,
#                 len_amplicon,avg_read_length,
#                 std_fragment_length,
#                 OUTPUT_DIRECTORY,log_filename)

            info('Running FLASH command: ' + cmd)
            crispresso2_info['flash_command'] = cmd
            FLASH_STATUS=sb.call(cmd,shell=True)
            if FLASH_STATUS:
                raise CRISPRessoShared.FlashException('Flash failed to run, please check the log file.')

            info('Done!')

            flash_hist_filename=_jp('out.hist')
            flash_histogram_filename=_jp('out.histogram')
            flash_not_combined_1_filename=_jp('out.notCombined_1.fastq.gz')
            flash_not_combined_2_filename=_jp('out.notCombined_2.fastq.gz')

            processed_output_filename=_jp('out.extendedFrags.fastq.gz')
            if os.path.isfile(processed_output_filename) is False:
                raise CRISPRessoShared.FlashException('Flash failed to produce merged reads file, please check the log file.')

        #count reads
        N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)
        if N_READS_AFTER_PREPROCESSING == 0:
            raise CRISPRessoShared.NoReadsAfterQualityFilteringException('No reads in input or no reads survived the average or single bp quality filtering.')

        info('Aligning sequences...')

        ####INITIALIZE CACHE####
        variantCache = {}
        #put empty sequence into cache
        cache_fastq_seq = ''
        variantCache[cache_fastq_seq] = {}
        variantCache[cache_fastq_seq]['count'] = 0

        #operates on variantCache
        aln_stats = process_fastq(processed_output_filename,variantCache,ref_names,refs,args)

        info('Done!')

        info('Quantifying indels/substitutions...')

        #ANALYZE ALIGNMENTS

        ###initialize
        N_TOTAL = 0
        N_DISCARDED = 0

        counts_total = {}
        counts_modified = {}
        counts_unmodified = {}
        counts_discarded = {}

        counts_insertion = {}
        counts_deletion = {}
        counts_substitution = {}

        counts_only_insertion = {}
        counts_only_deletion = {}
        counts_only_substitution = {}
        counts_insertion_and_deletion = {}
        counts_insertion_and_substitution = {}
        counts_deletion_and_substitution = {}
        counts_insertion_and_deletion_and_substitution = {}

        #we have insertions/deletions that change the concatenated exon sequence length and the difference between the final sequence
        #and the original sequence length is not a multiple of 3
        counts_modified_frameshift = {}

        #we have insertions/deletions that change the concatenated exon sequence length and the difference between the final sequence
        #and the original sequence length is a multiple of 3. We are in this case also when no indels are present but we have
        #substitutions
        counts_modified_non_frameshift = {}

        #read is modified but not at the exon
        counts_non_modified_non_frameshift = {}

        counts_splicing_sites_modified = {}

        ################
        class_counts = {} # number of reads in each class e.g. "ref1_UNMODIFIED" -> 50

        alleles_list = [] #will be turned into df with rows with information for each variant (allele)

        #for each reference, the following are computed individually
        all_insertion_count_vectors = {} #all insertions (including quantification window bases)
        all_insertion_left_count_vectors = {} #all insertions (including quantification window bases)
        all_deletion_count_vectors = {}
        all_substitution_count_vectors = {}
        all_indelsub_count_vectors = {}
        all_substitution_base_vectors = {}
        all_base_count_vectors = {} #number of times each base is seen

        insertion_count_vectors = {} #insertions that are in the quantification window
        deletion_count_vectors = {}
        substitution_count_vectors = {}
        indelsub_count_vectors = {}
        substitution_base_vectors = {} #these two are taken as a subset of all_substitution_base_vectors afterward to save time
        base_count_vectors = {} #calculated after as well

        insertion_count_vectors_noncoding = {} #insertions that are in the noncoding region
        deletion_count_vectors_noncoding = {}
        substitution_count_vectors_noncoding = {}


        insertion_length_vectors = {}
        deletion_length_vectors = {}

        inserted_n_lists = {} # list of number of insertions for all reads
        deleted_n_lists = {}
        substituted_n_lists = {}
        effective_len_lists = {} # list of effective lengths for all reads

        hists_inframe = {}
        hists_frameshift = {}

        #initialize data structures for each ref
        for ref_name in ref_names:
            this_len_amplicon = refs[ref_name]['sequence_length']

            counts_total                         [ref_name] = 0
            counts_modified                      [ref_name] = 0
            counts_unmodified                    [ref_name] = 0
            counts_discarded                     [ref_name] = 0
            counts_modified_frameshift           [ref_name] = 0
            counts_modified_non_frameshift       [ref_name] = 0
            counts_non_modified_non_frameshift   [ref_name] = 0
            counts_splicing_sites_modified       [ref_name] = 0

            counts_insertion                     [ref_name] = 0
            counts_deletion                      [ref_name] = 0
            counts_substitution                  [ref_name] = 0

            counts_only_insertion                [ref_name]  = 0
            counts_only_deletion                 [ref_name] = 0
            counts_only_substitution             [ref_name] = 0
            counts_insertion_and_deletion        [ref_name] = 0
            counts_insertion_and_substitution    [ref_name] = 0
            counts_deletion_and_substitution     [ref_name] = 0
            counts_insertion_and_deletion_and_substitution [ref_name] = 0

            all_insertion_count_vectors         [ref_name] = np.zeros(this_len_amplicon)
            all_insertion_left_count_vectors    [ref_name] = np.zeros(this_len_amplicon)
            all_deletion_count_vectors          [ref_name] = np.zeros(this_len_amplicon)
            all_substitution_count_vectors      [ref_name] = np.zeros(this_len_amplicon)
            all_indelsub_count_vectors          [ref_name] = np.zeros(this_len_amplicon)

            insertion_count_vectors             [ref_name] = np.zeros(this_len_amplicon)
            deletion_count_vectors              [ref_name] = np.zeros(this_len_amplicon)
            substitution_count_vectors          [ref_name] = np.zeros(this_len_amplicon)
            indelsub_count_vectors              [ref_name] = np.zeros(this_len_amplicon)

            insertion_count_vectors_noncoding   [ref_name] = np.zeros(this_len_amplicon)
            deletion_count_vectors_noncoding    [ref_name] = np.zeros(this_len_amplicon)
            substitution_count_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)

            #count times substitutions occur
            all_substitution_base_vectors        [ref_name+"_A" ] = np.zeros(this_len_amplicon)
            all_substitution_base_vectors        [ref_name+"_C" ] = np.zeros(this_len_amplicon)
            all_substitution_base_vectors        [ref_name+"_G" ] = np.zeros(this_len_amplicon)
            all_substitution_base_vectors        [ref_name+"_T" ] = np.zeros(this_len_amplicon)
            all_substitution_base_vectors        [ref_name+"_N" ] = np.zeros(this_len_amplicon)

            all_base_count_vectors               [ref_name+"_A" ] = np.zeros(this_len_amplicon)
            all_base_count_vectors               [ref_name+"_C" ] = np.zeros(this_len_amplicon)
            all_base_count_vectors               [ref_name+"_G" ] = np.zeros(this_len_amplicon)
            all_base_count_vectors               [ref_name+"_T" ] = np.zeros(this_len_amplicon)
            all_base_count_vectors               [ref_name+"_N" ] = np.zeros(this_len_amplicon)
            all_base_count_vectors               [ref_name+"_-" ] = np.zeros(this_len_amplicon)


            #length by position in amplicon
            insertion_length_vectors             [ref_name] = np.zeros(this_len_amplicon)
            deletion_length_vectors              [ref_name] = np.zeros(this_len_amplicon)


            inserted_n_lists                    [ref_name] = []
            deleted_n_lists                     [ref_name] = []
            substituted_n_lists                 [ref_name] = []
            effective_len_lists                 [ref_name] = []

            hists_inframe                       [ref_name] = defaultdict(lambda :0)
            hists_frameshift                    [ref_name] = defaultdict(lambda :0)
        #end initialize data structures for each ref

        #take care of empty seqs
        cache_fastq_seq = ''
        variantCache[cache_fastq_seq]['count'] = 0

        ###iterate through variants
        for variant in variantCache:
            #skip variant if there were none observed
            variantCount = variantCache[variant]['count']
            if (variantCount == 0):
                continue
            N_TOTAL += variantCount

            aln_ref_names = variantCache[variant]['aln_ref_names'] #list of references this seq aligned to
            aln_ref_scores = variantCache[variant]['aln_scores']
            class_name = variantCache[variant]['class_name'] #for classifying read e.g. 'HDR_MODIFIED' for pie chart

            if not args.expand_ambiguous_alignments and len(aln_ref_names) > 1: #got 'Ambiguous' -- don't count toward totals
                variantPayload = variantCache[variant]["variant_"+aln_ref_names[0]]
                alleleRow = {'#Reads':variantCount,
                           'Aligned_Sequence':variantPayload['aln_seq'],
                           'Reference_Sequence':variantPayload['aln_ref'],
                           'n_inserted':variantPayload['insertion_n'],
                           'n_deleted':variantPayload['deletion_n'],
                           'n_mutated':variantPayload['substitution_n'],
                           'Reference_Name':'AMBIGUOUS_'+aln_ref_names[0],
                           'Read_Status':variantPayload['classification'],
                           'Aligned_Reference_Names':"&".join(aln_ref_names),
                           'Aligned_Reference_Scores':"&".join([str(x) for x in aln_ref_scores]),
                           'ref_positions':variantPayload['ref_positions']
                }
                alleles_list.append(alleleRow)

                class_name = 'AMBIGUOUS'
                if class_name not in class_counts:
                        class_counts[class_name] = 0
                class_counts[class_name]+=variantCount

            else:
                if class_name not in class_counts:
                        class_counts[class_name] = 0
                class_counts[class_name]+=variantCount

                #iterate through payloads -- if a read aligned equally-well to two references, it could have more than one payload
                for ref_name in aln_ref_names:
                    variantPayload = variantCache[variant]["variant_"+ref_name]
                    if args.discard_indel_reads and (variantPayload['deletion_n'] > 0 or variantPayload['insertion_n'] > 0):
                        counts_discarded[ref_name] += variantCount
                        continue


                    counts_total[ref_name] += variantCount
                    if variantPayload['classification'] == 'MODIFIED':
                        counts_modified[ref_name] += variantCount
                    else:
                        counts_unmodified[ref_name] += variantCount

                    alleleRow = {'#Reads':variantCount,
                               'Aligned_Sequence':variantPayload['aln_seq'],
                               'Reference_Sequence':variantPayload['aln_ref'],
                               'n_inserted':variantPayload['insertion_n'],
                               'n_deleted':variantPayload['deletion_n'],
                               'n_mutated':variantPayload['substitution_n'],
                               'Reference_Name':variantPayload['ref_name'],
                               'Read_Status':variantPayload['classification'],
                               'Aligned_Reference_Names':"&".join(aln_ref_names),
                               'Aligned_Reference_Scores':"&".join([str(x) for x in aln_ref_scores]),
                               'ref_positions':variantPayload['ref_positions']
                    }
                    alleles_list.append(alleleRow)

                    this_effective_len= refs[ref_name]['sequence_length'] #how long is this alignment (insertions increase length, deletions decrease length)

                    this_has_insertions = False
                    all_insertion_count_vectors[ref_name][variantPayload['all_insertion_positions']]+=variantCount
                    all_insertion_left_count_vectors[ref_name][variantPayload['all_insertion_left_positions']]+=variantCount
                    all_indelsub_count_vectors[ref_name][variantPayload['all_insertion_positions']]+=variantCount

                    if not args.ignore_insertions:
                        inserted_n_lists[ref_name].extend([variantPayload['insertion_n']]*variantCount)
                        insertion_count_vectors[ref_name][variantPayload['insertion_positions']]+=variantCount
                        indelsub_count_vectors[ref_name][variantPayload['insertion_positions']]+=variantCount
                        this_effective_len = this_effective_len + variantPayload['insertion_n']
                        if variantPayload['insertion_n'] > 0:
                             counts_insertion[ref_name] += variantCount
                             this_has_insertions = True


                    this_has_deletions = False
                    all_deletion_count_vectors[ref_name][variantPayload['all_deletion_positions']]+=variantCount
                    all_indelsub_count_vectors[ref_name][variantPayload['all_deletion_positions']]+=variantCount
                    if not args.ignore_deletions:
                        deleted_n_lists[ref_name].extend([variantPayload['deletion_n']]*variantCount)
                        deletion_count_vectors[ref_name][variantPayload['deletion_positions']]+=variantCount
                        indelsub_count_vectors[ref_name][variantPayload['deletion_positions']]+=variantCount
                        this_effective_len = this_effective_len - variantPayload['deletion_n']
                        if variantPayload['deletion_n'] > 0:
                             counts_deletion[ref_name] += variantCount
                             this_has_deletions = True

                    effective_len_lists[ref_name].extend([this_effective_len]*variantCount)

                    this_has_substitutions = False
                    all_substitution_count_vectors[ref_name][variantPayload['all_substitution_positions']] += variantCount
                    all_indelsub_count_vectors[ref_name][variantPayload['all_substitution_positions']] += variantCount

                    if not args.ignore_substitutions:
                        substituted_n_lists[ref_name].extend([variantPayload['substitution_n']] * variantCount)
                        substitution_count_vectors[ref_name][variantPayload['substitution_positions']] += variantCount
                        indelsub_count_vectors[ref_name][variantPayload['substitution_positions']] += variantCount
                        if variantPayload['substitution_n'] > 0:
                             counts_substitution[ref_name] += variantCount
                             this_has_substitutions = True

                        nucs = ['A','T','C','G','N']
                        for nuc in nucs:
                            isNuc = variantPayload['all_substitution_values'] == ord(nuc)
                            if(np.sum(isNuc) > 0):
                                locs = np.array(variantPayload['all_substitution_positions'])[isNuc]
                                all_substitution_base_vectors[ref_name + "_" + nuc ][locs] += variantCount


                    if this_has_deletions:
                        if this_has_insertions:
                            if this_has_substitutions:
                                counts_insertion_and_deletion_and_substitution[ref_name] += variantCount
                            else:
                                counts_insertion_and_deletion[ref_name] += variantCount
                        else:
                            if this_has_substitutions:
                                counts_deletion_and_substitution[ref_name] += variantCount
                            else:
                                counts_only_deletion[ref_name] += variantCount
                    else: #no deletions
                        if this_has_insertions:
                            if this_has_substitutions:
                                counts_insertion_and_substitution[ref_name] += variantCount
                            else:
                                counts_only_insertion[ref_name] += variantCount
                        else:
                            if this_has_substitutions:
                                counts_only_substitution[ref_name] += variantCount

                    #set all_base_count_vectors
                    aln_seq = variantPayload['aln_seq']
                    ref_pos = variantPayload['ref_positions']
                    for i in range(len(aln_seq)):
                        if ref_pos[i] < 0:
                            continue
                        nuc = aln_seq[i]
                        all_base_count_vectors[ref_name + "_" + nuc][ref_pos[i]] += variantCount

                    exon_positions = refs[ref_name]['exon_positions']
                    exon_len_mods = refs[ref_name]['exon_len_mods'] #for each exon, how much length did this reference modify it?
                    tot_exon_len_mod = sum(exon_len_mods) #for all exons, how much length was modified?
                    splicing_positions = refs[ref_name]['splicing_positions']
                    insertion_coordinates = variantPayload['insertion_coordinates']
                    insertion_sizes = variantPayload['insertion_sizes']
                    all_insertion_positions = variantPayload['all_insertion_positions']
                    all_insertion_left_positions = variantPayload['all_insertion_left_positions']
                    insertion_positions = variantPayload['insertion_positions']
                    deletion_coordinates = variantPayload['deletion_coordinates']
                    deletion_sizes = variantPayload['deletion_sizes']
                    all_deletion_positions = variantPayload['all_deletion_positions']
                    deletion_positions = variantPayload['deletion_positions']
                    all_substitution_positions = variantPayload['all_substitution_positions']
                    substitution_positions = variantPayload['substitution_positions']


                    if this_has_insertions or this_has_deletions or this_has_substitutions or tot_exon_len_mod != 0: #only count modified reads
                        length_modified_positions_exons=[]
                        current_read_exons_modified = False
                        current_read_spliced_modified = False

                        for idx_ins,(ins_start,ins_end) in enumerate(insertion_coordinates):
                            insertion_length_vectors[ref_name][ins_start]+=(insertion_sizes[idx_ins]*variantCount)
                            insertion_length_vectors[ref_name][ins_end]+=(insertion_sizes[idx_ins]*variantCount)

                            if refs[ref_name]['contains_coding_seq']:
                                if set(exon_positions).intersection((ins_start, ins_end)): # check that we are inserting in one exon
                                    current_read_exons_modified = True
                                    set1 = set(exon_positions).intersection((ins_start, ins_end))
                                    length_modified_positions_exons.append((insertion_sizes[idx_ins]))

                        for idx_del, (del_start,del_end) in enumerate(deletion_coordinates):
                            deletion_length_vectors[ref_name][list(range(del_start,del_end))] += (deletion_sizes[idx_del]*variantCount)

                        if refs[ref_name]['contains_coding_seq']:
                            del_positions_to_append = sorted(set(exon_positions).intersection(set(deletion_positions)))
                            if del_positions_to_append:
                                current_read_exons_modified = True
                                length_modified_positions_exons.append(-len(del_positions_to_append))

                            if set(exon_positions).intersection(substitution_positions):
                                current_read_exons_modified = True

                            #splicing modifications
                            if set(splicing_positions).intersection(deletion_positions):
                                current_read_spliced_modified = True

                            if set(splicing_positions).intersection(insertion_positions):
                                current_read_spliced_modified = True

                            if set(splicing_positions).intersection(substitution_positions):
                                current_read_spliced_modified = True

                            if current_read_spliced_modified:
                                counts_splicing_sites_modified[ref_name] += variantCount

                            if tot_exon_len_mod != 0:
                                effective_length = sum(length_modified_positions_exons) + tot_exon_len_mod
                                if (effective_length % 3) == 0:
                                    # indels have restored exon frame
                                    counts_modified_non_frameshift[ref_name] += variantCount
                                    hists_inframe[ref_name][effective_length] += variantCount
                                else:
                                    counts_modified_frameshift[ref_name] += variantCount
                                    hists_frameshift[ref_name][effective_length] += variantCount
                            # if modified check if frameshift
                            elif current_read_exons_modified:

                                if not length_modified_positions_exons:
                                    # there are no indels
                                    counts_modified_non_frameshift[ref_name] += variantCount
                                    hists_inframe[ref_name][0] += variantCount
                                else:
                                    effective_length = sum(length_modified_positions_exons)

                                    if (effective_length % 3) == 0:
                                        counts_modified_non_frameshift[ref_name] += variantCount
                                        hists_inframe[ref_name][effective_length] += variantCount
                                    else:
                                        counts_modified_frameshift[ref_name] += variantCount
                                        hists_frameshift[ref_name][effective_length] += variantCount

                            # the indels and subtitutions are outside the exon/s  so we don't care!
                            else:
                                counts_non_modified_non_frameshift[ref_name] += variantCount
                                insertion_count_vectors_noncoding[ref_name][insertion_positions] += variantCount
                                deletion_count_vectors_noncoding[ref_name][deletion_positions] += variantCount
                                substitution_count_vectors_noncoding[ref_name][substitution_positions] += variantCount
                                hists_inframe[ref_name][0] += variantCount
                    #if unmodified but the tot_exon_len_mod is != 0
                    elif tot_exon_len_mod != 0:
                        if (tot_exon_len_mod % 3) == 0:
                            # indels have restored exon frame
                            counts_modified_non_frameshift[ref_name] += variantCount
                            hists_inframe[ref_name][effective_length] += variantCount
                        else:
                            counts_modified_frameshift[ref_name] += variantCount
                            hists_frameshift[ref_name][effective_length] += variantCount



        for ref_name in ref_names:
            this_include_idx = refs[ref_name]['include_idxs']
            substitution_base_vectors      [ref_name+"_A" ] = [all_substitution_base_vectors[ref_name+"_A"][x] for x in this_include_idx]
            substitution_base_vectors      [ref_name+"_C" ] = [all_substitution_base_vectors[ref_name+"_C"][x] for x in this_include_idx]
            substitution_base_vectors      [ref_name+"_G" ] = [all_substitution_base_vectors[ref_name+"_G"][x] for x in this_include_idx]
            substitution_base_vectors      [ref_name+"_T" ] = [all_substitution_base_vectors[ref_name+"_T"][x] for x in this_include_idx]
            substitution_base_vectors      [ref_name+"_N" ] = [all_substitution_base_vectors[ref_name+"_N"][x] for x in this_include_idx]

            base_count_vectors             [ref_name+"_A" ] = [all_base_count_vectors[ref_name+"_A"][x] for x in this_include_idx]
            base_count_vectors             [ref_name+"_C" ] = [all_base_count_vectors[ref_name+"_C"][x] for x in this_include_idx]
            base_count_vectors             [ref_name+"_G" ] = [all_base_count_vectors[ref_name+"_G"][x] for x in this_include_idx]
            base_count_vectors             [ref_name+"_T" ] = [all_base_count_vectors[ref_name+"_T"][x] for x in this_include_idx]
            base_count_vectors             [ref_name+"_N" ] = [all_base_count_vectors[ref_name+"_N"][x] for x in this_include_idx]
            base_count_vectors             [ref_name+"_-" ] = [all_base_count_vectors[ref_name+"_-"][x] for x in this_include_idx]

        # For HDR work, create a few more arrays, where all reads are aligned to ref1 (the first reference) and the indels are computed with regard to ref1
        if args.expected_hdr_amplicon_seq != "":
            ref1_name = ref_names[0]
            ref1_len = refs[ref1_name]['sequence_length']

            ref1_all_insertion_count_vectors = {} #all insertions (including quantification window bases) with respect to ref1
            ref1_all_deletion_count_vectors = {}
            ref1_all_substitution_count_vectors = {}
            for ref_name in ref_names:
                ref1_all_insertion_count_vectors[ref_name] = np.zeros(ref1_len)
                ref1_all_deletion_count_vectors[ref_name] = np.zeros(ref1_len)
                ref1_all_substitution_count_vectors[ref_name] = np.zeros(ref1_len)

            #for ref1 we will add all other indels to indels that have already been found..
            ref1_all_insertion_count_vectors[ref_names[0]] = all_insertion_count_vectors[ref_names[0]].copy()
            ref1_all_deletion_count_vectors[ref_names[0]] = all_deletion_count_vectors[ref_names[0]].copy()
            ref1_all_substitution_count_vectors[ref_names[0]] = all_substitution_count_vectors[ref_names[0]].copy()

            #then go through and align and add other reads
            for variant in variantCache:
                #skip variant if there were none observed
                variantCount = variantCache[variant]['count']
                if (variantCount == 0):
                    continue

                aln_ref_names = variantCache[variant]['aln_ref_names'] #list of references this seq aligned to
                if len(aln_ref_names) == 1 and ref_names[0] == aln_ref_names[0]: #if this read was only aligned to ref1, skip it because we already included the indels when we initialized the ref1_all_deletion_count_vectors array
                    continue

                #align this variant to ref1 sequence
                fws1,fws2,fwscore=CRISPResso2Align.global_align(variant, refs[ref1_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref1_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                rvs1,rvs2,rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(variant), refs[ref1_name]['sequence'],matrix=aln_matrix,gap_incentive=refs[ref1_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                s1 = fws1
                s2 = fws2
                score = fwscore
                if (rvscore > fwscore):
                    s1 = rvs1
                    s2 = rvs2
                    score = rvscore
                payload=CRISPRessoCOREResources.find_indels_substitutions(s1,s2,refs[ref1_name]['include_idxs'])

                #indels in this alignment against ref1 should be recorded for each ref it was originally assigned to, as well as for ref1
                #for example, if this read aligned to ref3, align this read to ref1, and add the resulting indels to ref1_all_insertion_count_vectors[ref3] as well as ref1_all_insertion_count_vectors[ref1]
                #   Thus, ref1_all_insertion_count_vectors[ref3] will show the position of indels of reads that aligned to ref3, but mapped onto ref1
                #   And ref1_alle_insertion_count_vectors[ref1] will show the position of indels of all reads, mapped onto ref1
                for ref_name in aln_ref_names:
                    if ref_name == ref_names[0]:
                        continue
                    ref1_all_insertion_count_vectors[ref_name][payload['all_insertion_positions']]+=variantCount
                    ref1_all_deletion_count_vectors[ref_name][payload['all_deletion_positions']]+=variantCount
                    ref1_all_substitution_count_vectors[ref_name][payload['all_substitution_positions']]+=variantCount

                #add these indel counts to the ref1
                ref1_all_insertion_count_vectors[ref_names[0]][payload['all_insertion_positions']]+=variantCount
                ref1_all_deletion_count_vectors[ref_names[0]][payload['all_deletion_positions']]+=variantCount
                ref1_all_substitution_count_vectors[ref_names[0]][payload['all_substitution_positions']]+=variantCount

        info('Done!')

        #order class_counts
        decorated_class_counts = []
        for class_count_name in class_counts:
            thisRefInd = 100
            thisIsMod = 1
            for idx,ref_name in enumerate(ref_names):
                if class_count_name.startswith(ref_name):
                    thisRefInd = idx
                    break
            if "UNMODIFIED" in class_count_name:
                thisIsMod = 0
            decorated_class_counts.append((thisRefInd,thisIsMod,class_count_name))
        decorated_class_counts.sort()
        class_counts_order = [class_count_name for thisRefInd,thisIsMod,class_count_name in decorated_class_counts]

        if N_TOTAL == 0:
            raise CRISPRessoShared.NoReadsAlignedException('Error: No alignments were found')

        #create alleles table
        info('Calculating allele frequencies...')

        #set up allele table
        df_alleles = pd.DataFrame(alleles_list)
        #df_alleles['%Reads']=df_alleles['#Reads']/df_alleles['#Reads'].sum()*100 # sum of #reads will be >= N_TOTAL because an allele appears once for each reference it aligns to
        df_alleles['%Reads']=df_alleles['#Reads']/N_TOTAL*100
        df_alleles[['n_deleted','n_inserted','n_mutated']] = df_alleles[['n_deleted','n_inserted','n_mutated']].astype(int)

        if np.sum(np.array(map(int,pd.__version__.split('.')))*(100,10,1))< 170:
           df_alleles.sort('#Reads',ascending=False,inplace=True)
        else:
           df_alleles.sort_values(by='#Reads',ascending=False,inplace=True)

        def calculate_range(l):
            try:
                r=max(15,int(np.round(np.percentile(l[np.nonzero(l)],99))))
            except:
                r=15
            return r


        all_insertion_pct_vectors = {} #all insertions/tot (including quantification window bases)
        all_deletion_pct_vectors = {}
        all_substitution_pct_vectors = {}
        all_indelsub_pct_vectors = {}

        insertion_pct_vectors = {} #insertions that are in the quantification window
        deletion_pct_vectors = {}
        substitution_pct_vectors = {}
        indelsub_pct_vectors = {}

        insertion_pct_vectors_noncoding = {} #insertions that are in the noncoding region
        deletion_pct_vectors_noncoding = {}
        substitution_pct_vectors_noncoding = {}

        for ref_name in ref_names:
            #save these values in the ref object-- we need to print them later

            ref_len = refs[ref_name]['sequence_length']
            if refs[ref_name]['contains_guide']:
                min_cut=min(refs[ref_name]['sgRNA_cut_points'])
                max_cut=max(refs[ref_name]['sgRNA_cut_points'])
                xmin,xmax=-min_cut,ref_len-max_cut
            else:
                min_cut=ref_len/2
                max_cut=ref_len/2
                xmin,xmax=-min_cut,+max_cut

            refs[ref_name]['xmin'] = xmin
            refs[ref_name]['xmax'] = xmax
            refs[ref_name]['min_cut'] = min_cut
            refs[ref_name]['max_cut'] = max_cut

            range_mut=calculate_range(substituted_n_lists[ref_name])
            range_ins=calculate_range(inserted_n_lists[ref_name])
            range_del=calculate_range(deleted_n_lists[ref_name])

            y_values_mut,x_bins_mut=plt.histogram(substituted_n_lists[ref_name],bins=range(0,range_mut))
            y_values_ins,x_bins_ins=plt.histogram(inserted_n_lists[ref_name],bins=range(0,range_ins))
            y_values_del,x_bins_del=plt.histogram(deleted_n_lists[ref_name],bins=range(0,range_del))

            refs[ref_name]['y_values_mut'] = y_values_mut
            refs[ref_name]['x_bins_mut'] = x_bins_mut
            refs[ref_name]['y_values_ins'] = y_values_ins
            refs[ref_name]['x_bins_ins'] = x_bins_ins
            refs[ref_name]['y_values_del'] = y_values_del
            refs[ref_name]['x_bins_del'] = x_bins_del

            vals = np.array(effective_len_lists[ref_name]) - ref_len
            hdensity,hlengths=np.histogram(vals,np.arange(xmin,xmax))
            hlengths=hlengths[:-1]
            center_index=np.nonzero(hlengths==0)[0][0]

            refs[ref_name]['hdensity'] = hdensity
            refs[ref_name]['hlengths'] = hlengths
            refs[ref_name]['center_index'] = center_index

            if not dict(hists_inframe[ref_name]):
                hists_inframe[ref_name]={0:0}
            else:
                hists_inframe[ref_name]=dict(hists_inframe[ref_name])

            if not dict(hists_frameshift[ref_name]):
                hists_frameshift[ref_name]={0:0}
            else:
                hists_frameshift[ref_name]=dict(hists_frameshift[ref_name])

            count_tot = counts_total[ref_name]
            if count_tot > 0:
                #normalize effect vectors
                all_insertion_pct_vectors[ref_name] = 100 * all_insertion_count_vectors[ref_name] / count_tot
                all_deletion_pct_vectors[ref_name] = 100 * all_deletion_count_vectors[ref_name] / count_tot
                all_substitution_pct_vectors[ref_name] = 100 * all_substitution_count_vectors[ref_name] / count_tot
                all_indelsub_pct_vectors[ref_name] = 100 * all_indelsub_count_vectors[ref_name] / count_tot

                insertion_pct_vectors[ref_name] = 100 * insertion_count_vectors[ref_name] / count_tot
                deletion_pct_vectors[ref_name] = 100 * deletion_count_vectors[ref_name] / count_tot
                substitution_pct_vectors[ref_name] = 100 * substitution_count_vectors[ref_name] / count_tot
                indelsub_pct_vectors[ref_name] = 100 * indelsub_count_vectors[ref_name] / count_tot

                insertion_pct_vectors_noncoding[ref_name] = 100 * insertion_count_vectors_noncoding[ref_name] / count_tot
                deletion_pct_vectors_noncoding[ref_name] = 100 * deletion_count_vectors_noncoding[ref_name] / count_tot
                substitution_pct_vectors_noncoding[ref_name] = 100 * substitution_count_vectors_noncoding[ref_name] / count_tot

                deletion_length_vectors_sum = deletion_length_vectors[ref_name]
                deletion_length_vectors[ref_name] = np.zeros(ref_len)
                delMask = deletion_count_vectors[ref_name] > 0
                deletion_length_vectors[ref_name][delMask] = deletion_length_vectors_sum[delMask]/deletion_count_vectors[ref_name][delMask]

                insertion_length_vectors_sum = insertion_length_vectors[ref_name]
                insertion_length_vectors[ref_name] = np.zeros(ref_len)
                insMask = insertion_count_vectors[ref_name] > 0
                insertion_length_vectors[ref_name][insMask] = insertion_length_vectors_sum[insMask]/insertion_count_vectors[ref_name][insMask]
            else: #no reads for this ref
                this_len_amplicon = refs[ref_name]['sequence_length']
                all_insertion_pct_vectors[ref_name] = np.zeros(this_len_amplicon)
                all_deletion_pct_vectors[ref_name]  = np.zeros(this_len_amplicon)
                all_substitution_pct_vectors[ref_name] = np.zeros(this_len_amplicon)
                all_indelsub_pct_vectors[ref_name] = np.zeros(this_len_amplicon)

                insertion_pct_vectors[ref_name] = np.zeros(this_len_amplicon)
                deletion_pct_vectors[ref_name] = np.zeros(this_len_amplicon)
                substitution_pct_vectors[ref_name] = np.zeros(this_len_amplicon)
                indelsub_pct_vectors[ref_name] = np.zeros(this_len_amplicon)

                insertion_pct_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)
                deletion_pct_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)
                substitution_pct_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)

                deletion_length_vectors[ref_name] = np.zeros(ref_len)

                insertion_length_vectors[ref_name] = np.zeros(ref_len)

        info('Done!')

        if args.dump:
            ref_info_file_name = _jp('CRISPResso_reference_info.txt')
            ref_info_file = open(ref_info_file_name,'w')
            refString = ( 'name' + "\t" +
                'sequence' + "\t" +
                'sequence_length' + "\t" +
                'min_aln_score' + "\t" +
                'gap_incentive' + "\t" +
                'sgRNA_cut_points' + "\t" +
                'sgRNA_intervals' + "\t" +
                'sgRNA_plot_idxs' + "\t" +
                'sgRNA_sequences' + "\t" +
                'contains_guide' + "\t" +
                'contains_coding_seq' + "\t" +
                'exon_positions' + "\t" +
                'exon_intervals' + "\t" +
                'exon_len_mods' + "\t" +
                'splicing_positions' + "\t" +
                'include_idxs' + "\t" +
                'exclude_idxs' + "\t" +
                'idx_cloned_from' + "\n")
            ref_info_file.write(refString)
            np.set_printoptions(linewidth=1000**1000) #no line breaks
            for ref_name in ref_names:
                refString = ( refs[ref_name]['name'] + "\t" +
                    str(refs[ref_name]['sequence']) + "\t" +
                    str(refs[ref_name]['sequence_length']) + "\t" +
                    str(refs[ref_name]['min_aln_score']) + "\t" +
                    str(refs[ref_name]['gap_incentive']) + "\t" +
                    str(refs[ref_name]['sgRNA_cut_points']) + "\t" +
                    str(refs[ref_name]['sgRNA_intervals']) + "\t" +
                    str(refs[ref_name]['sgRNA_sequences']) + "\t" +
                    str(refs[ref_name]['sgRNA_plot_idxs']) + "\t" +
                    str(refs[ref_name]['contains_guide']) + "\t" +
                    str(refs[ref_name]['contains_coding_seq']) + "\t" +
                    str(refs[ref_name]['exon_positions']) + "\t" +
                    str(refs[ref_name]['exon_intervals']) + "\t" +
                    str(refs[ref_name]['exon_len_mods']) + "\t" +
                    str(refs[ref_name]['splicing_positions']) + "\t" +
                    str(refs[ref_name]['include_idxs']) + "\t" +
                    str(refs[ref_name]['exclude_idxs']) + "\t" +
                    str(refs[ref_name]['idx_cloned_from']) + "\n")
                ref_info_file.write(refString)
            ref_info_file.close()

        crispresso2_info['ref_names'] = ref_names
        crispresso2_info['refs'] = refs

        info('Saving processed data...')

        #write alleles table
        #crispresso1Cols = ["Aligned_Sequence","Reference_Sequence","NHEJ","UNMODIFIED","HDR","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        #df_alleles.ix[:,crispresso1Cols].to_csv(_jp('Alleles_frequency_table.txt'),sep='\t',header=True,index=None)
        #crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
#        crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads","Aligned_Reference_Names","Aligned_Reference_Scores"]
#        crispresso2Cols = ["Read_Sequence","Amplicon_Sequence","Amplicon_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]

        allele_frequency_table_filename = 'Alleles_frequency_table.txt'
        allele_frequency_table_fileLoc = _jp(allele_frequency_table_filename)

        allele_frequency_table_zip_filename = _jp('Alleles_frequency_table.zip')

        df_alleles.ix[:,crispresso2Cols].to_csv(allele_frequency_table_fileLoc,sep='\t',header=True,index=None)
        with zipfile.ZipFile(allele_frequency_table_zip_filename,'w',zipfile.ZIP_DEFLATED) as myzip:
            myzip.write(allele_frequency_table_fileLoc,allele_frequency_table_filename)
        os.remove(allele_frequency_table_fileLoc)
        crispresso2_info['allele_frequency_table_filename'] = os.path.basename(allele_frequency_table_filename) #filename is the name of the file in the zip
        crispresso2_info['allele_frequency_table_zip_filename'] = os.path.basename(allele_frequency_table_zip_filename)

        if args.crispresso1_mode:
            with open(_jp('Quantification_of_editing_frequency.txt'),'w+') as outfile:
                outfile.write("Quantification of editing frequency:\n")
                for ref_name in ref_names:
                    n_unmod = counts_unmodified[ref_name]
                    n_mod = counts_modified[ref_name]
                    n_discarded = counts_discarded[ref_name]

                    n_insertion = counts_insertion[ref_name]
                    n_deletion = counts_deletion[ref_name]
                    n_substitution = counts_substitution[ref_name]

                    outfile.write("%s: Unmodified: %d Modified: %d Discarded: %d\n" % (ref_name,n_unmod,n_mod,n_discarded))
                    outfile.write("(%d reads with insertions, %d reads with deletions, %d reads with substitutions)\n" % (n_insertion,n_deletion,n_substitution))

                outfile.write('Total Aligned:%d reads ' % N_TOTAL)

        quant_of_editing_freq_filename =_jp('CRISPResso_quantification_of_editing_frequency.txt')
        with open(quant_of_editing_freq_filename,'w+') as outfile:
            outfile.write('Amplicon\tUnmodified%\tModified%\tReads_aligned\tReads_total\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions\n')
            for ref_name in ref_names:
                n_aligned = counts_total[ref_name]
                n_unmod = counts_unmodified[ref_name]
                n_mod = counts_modified[ref_name]
                n_discarded = counts_discarded[ref_name]

                n_insertion = counts_insertion[ref_name]
                n_deletion = counts_deletion[ref_name]
                n_substitution = counts_substitution[ref_name]
                n_only_insertion = counts_only_insertion[ref_name]
                n_only_deletion = counts_only_deletion[ref_name]
                n_only_substitution = counts_only_substitution[ref_name]
                n_insertion_and_deletion = counts_insertion_and_deletion[ref_name]
                n_insertion_and_substitution = counts_insertion_and_substitution[ref_name]
                n_deletion_and_substitution = counts_deletion_and_substitution[ref_name]
                n_insertion_and_deletion_and_substitution = counts_insertion_and_deletion_and_substitution[ref_name]

                unmod_pct = "NA"
                mod_pct = "NA"
                if n_aligned > 0:
                    unmod_pct = 100*n_unmod/float(n_aligned)
                    mod_pct = 100*n_mod/float(n_aligned)

                vals = [ref_name]
                vals.extend([str(x) for x in [round(unmod_pct,8),round(mod_pct,8),n_aligned,N_TOTAL,n_unmod,n_mod,n_discarded,n_insertion,n_deletion,n_substitution,n_only_insertion,n_only_deletion,n_only_substitution,n_insertion_and_deletion,n_insertion_and_substitution,n_deletion_and_substitution,n_insertion_and_deletion_and_substitution]])
                outfile.write("\t".join(vals) + "\n")

        crispresso2_info['quant_of_editing_freq_filename'] = os.path.basename(quant_of_editing_freq_filename)


        #write statistics
        if args.crispresso1_mode:
            with open(_jp('Mapping_statistics.txt'),'w+') as outfile:
                outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d\n' % (N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_TOTAL))

        mapping_stats_filename = _jp('CRISPResso_mapping_statistics.txt')
        with open(mapping_stats_filename,'w+') as outfile:
            outfile.write('READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n')
            outfile.write("\t".join([str(x) for x in[N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_TOTAL,aln_stats['N_COMPUTED_ALN'],aln_stats['N_CACHED_ALN'],aln_stats['N_COMPUTED_NOTALN'],aln_stats['N_CACHED_NOTALN']]]) + "\n")
        crispresso2_info['aln_stats'] = aln_stats
        crispresso2_info['mapping_stats_filename'] = os.path.basename(mapping_stats_filename)

        def save_vector_to_file(vector,filename):
            #np.savetxt(_jp('%s.txt' %name), np.vstack([(np.arange(len(vector))+1),vector]).T, fmt=['%d','%.18e'],delimiter='\t', newline='\n', header='amplicon position\teffect',footer='', comments='# ')
            np.savetxt(filename, np.vstack([(np.arange(len(vector))+1),vector]).T, fmt=['%d','%.18e'],delimiter='\t', newline='\n', header='amplicon position\teffect',footer='', comments='# ')

        def save_count_vectors_to_file(vectors,vectorNames,refSeq,filename):
            outfile = open(filename,"w")
            outfile.write("Sequence\t"+"\t".join(list(refSeq))+"\n") #first row: reference sequence
            for vector,vectorName in zip(vectors,vectorNames):
                outfile.write(vectorName +"\t" + "\t".join([str(x) for x in vector]) + "\n") #next, vectors are printed
            outfile.close()

        crispresso2_info['insertion_pct_vectors'] = insertion_pct_vectors
        crispresso2_info['deletion_pct_vectors'] = insertion_pct_vectors
        crispresso2_info['substitution_pct_vectors'] = insertion_pct_vectors
        crispresso2_info['indelsub_pct_vectors'] = insertion_pct_vectors

        for ref_name in ref_names:
            #only show reference name in filenames if more than one reference
            ref_plot_name = ref_name+"."
            if len(ref_names) == 1 and ref_names[0] == "Reference":
                ref_plot_name = ""

            #n_this_category = counts_total[ref_name]
            #if n_this_category < 1:
            #    continue

            if not args.suppress_plots:
                ins_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_insertion.txt')
                save_vector_to_file(insertion_pct_vectors[ref_name],ins_pct_vector_filename)
                crispresso2_info['refs'][ref_name]['insertion_pct_vector_filename'] = os.path.basename(ins_pct_vector_filename)

                del_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_deletion.txt')
                save_vector_to_file(deletion_pct_vectors[ref_name],del_pct_vector_filename)
                crispresso2_info['refs'][ref_name]['deletion_pct_vector_filename'] = os.path.basename(del_pct_vector_filename)

                sub_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_substitution.txt')
                save_vector_to_file(substitution_pct_vectors[ref_name],sub_pct_vector_filename)
                crispresso2_info['refs'][ref_name]['substitution_pct_vector_filename'] = os.path.basename(sub_pct_vector_filename)

                indelsub_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_combined.txt')
                save_vector_to_file(indelsub_pct_vectors[ref_name],indelsub_pct_vector_filename)
                crispresso2_info['refs'][ref_name]['combined_pct_vector_filename'] = os.path.basename(indelsub_pct_vector_filename)

            #save mods in quantification window
            quant_window_mod_count_filename = _jp(ref_plot_name+'Quantification_window_modification_count_vectors.txt')
            save_count_vectors_to_file([insertion_count_vectors[ref_name],
                        deletion_count_vectors[ref_name],
                        substitution_count_vectors[ref_name],
                        indelsub_count_vectors[ref_name],
                        [counts_total[ref_name]]*refs[ref_name]['sequence_length']],
                        ['Insertions','Deletions','Substitutions','All_modifications','Total'],
                            refs[ref_name]['sequence'],quant_window_mod_count_filename)
            crispresso2_info['refs'][ref_name]['quant_window_mod_count_filename'] = os.path.basename(quant_window_mod_count_filename)

            #save all mods
            mod_count_filename = _jp(ref_plot_name+'Modification_count_vectors.txt')
            save_count_vectors_to_file([all_insertion_count_vectors[ref_name],
                        all_insertion_left_count_vectors[ref_name],
                        all_deletion_count_vectors[ref_name],
                        all_substitution_count_vectors[ref_name],
                        all_indelsub_count_vectors[ref_name],
                        [counts_total[ref_name]]*refs[ref_name]['sequence_length']],
                        ['Insertions','Insertions_Left','Deletions','Substitutions','All_modifications','Total'],
                            refs[ref_name]['sequence'],mod_count_filename)
            crispresso2_info['refs'][ref_name]['mod_count_filename'] = os.path.basename(mod_count_filename)
            crispresso2_info['refs'][ref_name]['mod_count_filename_caption'] = "A tab-separated file showing the number of modifications for each position in the amplicon. " \
                "The first row shows the amplicon sequence, and successive rows show the number of reads with insertions (row 2), insertions_left (row 3), deletions (row 4), substitutions (row 5) and the sum of all modifications (row 6)." \
                "Additionally, the last row shows the number of reads aligned. If an insertion occurs between bases 5 and 6, the insertions vector will be incremented at bases 5 and 6. " \
                "However, the insertions_left vector will only be incremented at base 5 so the sum of the insertions_left row represents an accurate count of the number of insertions, " \
                "whereas the sum of the insertions row will yield twice the number of insertions."

            if (refs[ref_name]['contains_coding_seq']): #PERFORM FRAMESHIFT ANALYSIS
                MODIFIED_FRAMESHIFT = counts_modified_frameshift[ref_name]
                MODIFIED_NON_FRAMESHIFT = counts_modified_non_frameshift[ref_name]
                NON_MODIFIED_NON_FRAMESHIFT = counts_non_modified_non_frameshift[ref_name]
                SPLICING_SITES_MODIFIED = counts_splicing_sites_modified[ref_name]
                frameshift_analysis_filename = _jp(ref_plot_name+'Frameshift_analysis.txt')
                with open(frameshift_analysis_filename,'w+') as outfile:
                        outfile.write('Frameshift analysis:\n\tNoncoding mutation:%d reads\n\tIn-frame mutation:%d reads\n\tFrameshift mutation:%d reads\n' %(NON_MODIFIED_NON_FRAMESHIFT, MODIFIED_NON_FRAMESHIFT ,MODIFIED_FRAMESHIFT))
                crispresso2_info['refs'][ref_name]['frameshift_analysis_filename'] = os.path.basename(frameshift_analysis_filename)
                crispresso2_info['refs'][ref_name]['frameshift_analysis_filename_caption'] = "A text file describing the number of noncoding, in-frame, and frameshift mutations. This report file is produced when the amplicon contains a coding sequence."

                splice_sites_analysis_filename = _jp(ref_plot_name+'Splice_sites_analysis.txt')
                with open(splice_sites_analysis_filename,'w+') as outfile:
                        outfile.write('Splice sites analysis:\n\tUnmodified:%d reads\n\tPotential splice sites modified:%d reads\n' %(counts_total[ref_name]- SPLICING_SITES_MODIFIED, SPLICING_SITES_MODIFIED))
                crispresso2_info['refs'][ref_name]['splice_sites_analysis_filename'] = os.path.basename(splice_sites_analysis_filename)
                crispresso2_info['refs'][ref_name]['splice_sites_analysis_filename_caption'] = "A text file describing the number of splicing sites that are unmodified and modified. This file report is produced when the amplicon contains a coding sequence."

                ins_pct_vector_noncoding_filename = _jp(ref_plot_name+'Effect_vector_insertion_noncoding.txt')
                save_vector_to_file(insertion_pct_vectors_noncoding[ref_name],ins_pct_vector_noncoding_filename)
                crispresso2_info['refs'][ref_name]['insertion_pct_vector_noncoding_filename'] = os.path.basename(ins_pct_vector_noncoding_filename)
                crispresso2_info['refs'][ref_name]['insertion_pct_vector_noncoding_filename_caption'] = "A tab-separated text file with a one-row header that shows the percentage of reads with a noncoding insertion at each base in the " + ref_name + " sequence. " \
                    "The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding insertion at that location. This report file is produced when the amplicon contains a coding sequence."

                del_pct_vector_noncoding_filename = _jp(ref_plot_name+'Effect_vector_deletion_noncoding.txt')
                save_vector_to_file(deletion_pct_vectors_noncoding[ref_name],del_pct_vector_noncoding_filename)
                crispresso2_info['refs'][ref_name]['deletion_pct_vector_noncoding_filename'] = os.path.basename(del_pct_vector_noncoding_filename)
                crispresso2_info['refs'][ref_name]['deletion_pct_vector_noncoding_filename_caption'] = "A tab-separated text file with a one-row header that shows the percentage of reads with a noncoding deletion at each base in the " + ref_name + " sequence. " \
                    "The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding deletion at that location. This report file is produced when the amplicon contains a coding sequence."

                sub_pct_vector_noncoding_filename = _jp(ref_plot_name+'Effect_vector_substitution_noncoding.txt')
                save_vector_to_file(substitution_pct_vectors_noncoding[ref_name],sub_pct_vector_noncoding_filename)
                crispresso2_info['refs'][ref_name]['substitution_pct_vector_noncoding_filename'] = os.path.basename(sub_pct_vector_noncoding_filename)
                crispresso2_info['refs'][ref_name]['substitution_pct_vector_noncoding_filename_caption'] = "A tab-separated text file with a one-row header that shows the percentage of reads with a noncoding substitution at each base in the " + ref_name + " sequence. " \
                    "The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a nondcoding substitution at that location. This report file is produced when the amplicon contains a coding sequence."

            if args.dump:
                if refs[ref_name]['sgRNA_cut_points']:
                    cp.dump(refs[ref_name]['sgRNA_cut_points'], open( _jp(ref_plot_name+'Cut_points.pickle'), 'wb' ) )

                if refs[ref_name]['sgRNA_intervals']:
                    cp.dump(refs[ref_name]['sgRNA_intervals'], open( _jp(ref_plot_name+'sgRNA_intervals.pickle'), 'wb' ) )

            hdensity = refs[ref_name]['hdensity']
            hlengths = refs[ref_name]['hlengths']
            center_index = refs[ref_name]['center_index']

            y_values_mut = refs[ref_name]['y_values_mut']
            x_bins_mut = refs[ref_name]['x_bins_mut']
            y_values_ins = refs[ref_name]['y_values_ins']
            x_bins_ins = refs[ref_name]['x_bins_ins']
            y_values_del = refs[ref_name]['y_values_del']
            x_bins_del = refs[ref_name]['x_bins_del']

            if not args.suppress_plots:

                indel_histogram_file = _jp(ref_plot_name+'Indel_histogram.txt')
                pd.DataFrame(np.vstack([hlengths,hdensity]).T,columns=['indel_size','fq']).to_csv(indel_histogram_file,index=None,sep='\t')
                crispresso2_info['refs'][ref_name]['indel_histogram_filename'] = os.path.basename(indel_histogram_file)
                crispresso2_info['refs'][ref_name]['indel_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the length of indels (both insertions and deletions) in the " + ref_name +" sequence in the quantification window. " \
                            "Indels outside of the quantification window are not included. The indel_size column shows the number of substitutions, and the fq column shows the number of reads having an indel of that length."



                insertion_histogram_file = _jp(ref_plot_name+'Insertion_histogram.txt')
                pd.DataFrame(np.vstack([x_bins_ins[:-1],y_values_ins]).T,columns=['ins_size','fq']).to_csv(insertion_histogram_file,index=None,sep='\t')
                crispresso2_info['refs'][ref_name]['insertion_histogram_filename'] = os.path.basename(insertion_histogram_file)
                crispresso2_info['refs'][ref_name]['insertion_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of insertions in the " + ref_name +" sequence in the quantification window. " \
                    "Insertions outside of the quantification window are not included. The ins_size column shows the number of insertions, and the fq column shows the number of reads having that number of insertions."


                deletion_histogram_file = _jp(ref_plot_name+'Deletion_histogram.txt')
                pd.DataFrame(np.vstack([-x_bins_del[:-1],y_values_del]).T,columns=['del_size','fq']).to_csv(deletion_histogram_file,index=None,sep='\t')
                crispresso2_info['refs'][ref_name]['deletion_histogram_filename'] = os.path.basename(deletion_histogram_file)
                crispresso2_info['refs'][ref_name]['deletion_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of deletions in the " + ref_name +" sequence in the quantification window. " \
                "Deletions outside of the quantification window are not included. The del_size column shows the number of deletions, and the fq column shows the number of reads having that number of deletions."


                substitution_histogram_file = _jp(ref_plot_name+'Substitution_histogram.txt')
                pd.DataFrame(np.vstack([x_bins_mut[:-1],y_values_mut]).T,columns=['sub_count','fq']).to_csv(substitution_histogram_file,index=None,sep='\t')
                crispresso2_info['refs'][ref_name]['substitution_histogram_filename'] = os.path.basename(substitution_histogram_file)
                crispresso2_info['refs'][ref_name]['substitution_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of substitutions in the " + ref_name +" sequence in the quantification window. " \
                "Substitutions outside of the quantification window are not included. The sub_size column shows the number of substitutions, and the fq column shows the number of reads having that number of substitutions."



            if args.dump:
                np.savez(_jp(ref_plot_name+'Effect_vector_insertion'),insertion_pct_vectors[ref_name])
                np.savez(_jp(ref_plot_name+'Effect_vector_deletion'),deletion_pct_vectors[ref_name])
                np.savez(_jp(ref_plot_name+'Effect_vector_substitution'),substitution_pct_vectors[ref_name])

                np.savez(_jp(ref_plot_name+'Effect_vector_combined'),indelsub_pct_vectors[ref_name])

                np.savez(_jp(ref_plot_name+'Position_dependent_vector_avg_insertion_size'),insertion_length_vectors[ref_name])
                np.savez(_jp(ref_plot_name+'Position_dependent_vector_avg_deletion_size'),deletion_length_vectors[ref_name])

        if args.dump:
            info('Dumping all the processed data...')

        if not args.suppress_plots:
            info('Making Plots...')
        ###############################################################################################################################################
        save_png = True
        if args.suppress_report:
            save_png = False

        n_refs = len(ref_names)
        #helper function .. if there is only one reference, don't print the name on the top of every plot
        def get_plot_title_with_ref_name(plotTitle,ref_name):
            if n_refs > 1:
                return (plotTitle + ": " + ref_name)
            return plotTitle

        def count_alternate_alleles(sub_base_vectors,ref_name,ref_sequence,ref_total_aln_reads):
            #create vectors with all allele frequencies -- not just the substitution (the reference allele will not be 0)
            alph = ['A','C','G','T','N']

            #count the total number of times each substitution occurs
            count_sub_base_vectors = {}
            alt_nuc_counts = {}
            for a in alph:
                alt_nuc_counts[a] = {}
                count_sub_base_vectors[a] = list(sub_base_vectors[ref_name+"_"+a])
                for b in alph:
                    alt_nuc_counts[a][b] = 0

            for idx,c in enumerate(ref_sequence):
                tot_sub_at_idx = 0
                for a in alph:
                    sub = sub_base_vectors[ref_name+"_" + a][idx]
                    alt_nuc_counts[c][a] += sub
                    tot_sub_at_idx += sub

            #df_subs = pd.DataFrame([count_sub_base_vectors["A"],count_sub_base_vectors["C"],count_sub_base_vectors["G"],count_sub_base_vectors["T"],count_sub_base_vectors["N"]])
            df_subs = pd.DataFrame([count_sub_base_vectors[a] for a in alph])
            df_subs.index = alph
            df_subs.columns = list(ref_sequence)
            return (df_subs,alt_nuc_counts)

            ############

        ###############################################################################################################################################
        ### FIGURE 1: Alignment
        #fig=plt.figure(figsize=(12*1.5,12*1.5))
        if not args.suppress_plots:
            fig=plt.figure(figsize=(12,12))
            ax = plt.subplot(111)
            labels=['READS\nIN INPUTS\n(%d)'%N_READS_INPUT,'READS AFTER\nPREPROCESSING\n(%d)'%N_READS_AFTER_PREPROCESSING,'READS\nALIGNED\n(%d)'%N_TOTAL]
            sizes = [N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_TOTAL]
            rects = ax.bar(np.arange(len(sizes)),sizes,color='silver')
            #label each bar
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()/2, height + 0.05,height, ha='center', va='bottom')

            ax.set_xticks(np.arange(len(sizes)))
            ax.set_xticklabels(labels)
            ax.set_ylabel('Sequences % (no.)')
            y_label_values= np.round(np.linspace(0, max(N_READS_INPUT,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
            ax.set_yticks(y_label_values)
            ax.set_yticklabels(['%.1f%% (%.0f)' % (100*cnt/N_READS_INPUT,cnt) for cnt in y_label_values])
            #if too many barplots, flip the labels
            plt.ylim(0,max(sizes)*1.1)

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            plt.tight_layout()

            plot_root = _jp("1a.Read_barplot")
            plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
            if save_png:
                plt.savefig(plot_root+'.png',bbox_inches='tight')
            plt.close()

            crispresso2_info['plot_1a_root'] = os.path.basename(plot_root)
            crispresso2_info['plot_1a_caption'] = "Figure 1a: The number of reads in input fastqs, after preprocessing, and after alignment to amplicons."
            crispresso2_info['plot_1a_data'] = [('Mapping statistics',os.path.basename(mapping_stats_filename))]


            #(1b) a piechart of classes
            labels = []
            sizes = []
            for class_name in class_counts_order:
                if args.expected_hdr_amplicon_seq != "" and class_name == ref_names[0]+"_MODIFIED":
                    labels.append("NHEJ" + "\n(" + str(class_counts[class_name]) + " reads)")
                elif args.expected_hdr_amplicon_seq != "" and class_name == "HDR_MODIFIED":
                    labels.append("Imperfect HDR" + "\n(" + str(class_counts[class_name]) + " reads)")
                elif args.expected_hdr_amplicon_seq != "" and class_name == "HDR_UNMODIFIED":
                    labels.append("HDR" + "\n(" + str(class_counts[class_name]) + " reads)")
                else:
                    display_class_name = class_name
                    if len(ref_names) == 1:
                        display_class_name = display_class_name.replace('Reference_','')

                    labels.append(display_class_name + "\n(" + str(class_counts[class_name]) + " reads)")

                sizes.append(100*class_counts[class_name]/float(N_TOTAL))


            #fig=plt.figure(figsize=(8.3,8))
            fig=plt.figure(figsize=(12,12))
            ax = plt.subplot(111)
            patches, texts, autotexts =ax.pie(sizes, labels=labels,\
                            autopct='%1.2f%%')

            plt.axis('off')

            proptease = fm.FontProperties()
    #        proptease.set_size('x-large')
    #        plt.setp(autotexts, fontproperties=proptease)
    #        plt.setp(texts, fontproperties=proptease)
            plt.axis("equal")

            plot_root = _jp("1b.Alignment_pie_chart")
            plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
            if save_png:
                plt.savefig(plot_root+'.png',bbox_inches='tight')
            plt.close()

            crispresso2_info['plot_1b_root'] = os.path.basename(plot_root)
            crispresso2_info['plot_1b_caption'] = "Figure 1b: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles."
            crispresso2_info['plot_1b_data'] = [('Quantification of editing',os.path.basename(quant_of_editing_freq_filename))]

            #(1c) a barchart of classes
            fig=plt.figure(figsize=(12,12))
            ax = plt.subplot(111)

            rects = ax.bar(np.arange(len(sizes)),sizes,color='silver')
            #label each bar
            for rect in rects:
                height = rect.get_height()
                if len(sizes) > 4:
                    #ax.text(rect.get_x() + rect.get_width()/2, height + 0.05,"%.2f%%"%height, ha='center', va='bottom',fontsize=12)
                    ax.text(rect.get_x() + rect.get_width()/2, height + 0.05,"%.2f%%"%height, ha='center', va='bottom')
                else:
                    ax.text(rect.get_x() + rect.get_width()/2, height + 0.05,"%.2f%%"%height, ha='center', va='bottom')


            ax.set_xticks(np.arange(len(sizes)))

            ax.set_xticklabels(labels)
            ax.set_xticklabels([X.replace("&"," & ").replace("_UNMODIFIED"," UNMODIFIED").replace("_MODIFIED"," MODIFIED") for X in labels])
            ax.set_ylabel('Sequences % (no.)')
            y_label_values= np.round(np.linspace(0, min(100,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
            ax.set_yticks(y_label_values)
            ax.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*N_TOTAL) for pct in y_label_values])
            #if too many barplots, flip the labels
            if len(sizes) > 4:
                #plt.setp(ax.get_xticklabels(), fontsize=12, rotation='vertical',multialignment='right')
                plt.setp(ax.get_xticklabels(),rotation='vertical',multialignment='right')
            plt.ylim(0,max(sizes)*1.1)

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            plt.tight_layout()

            plot_root = _jp('1c.Alignment_barplot')
            plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
            if save_png:
                plt.savefig(plot_root+'.png',bbox_inches='tight')
            plt.close()

            crispresso2_info['plot_1c_root'] = os.path.basename(plot_root)
            crispresso2_info['plot_1c_caption'] = "Figure 1c: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles."
            crispresso2_info['plot_1c_data'] = [('Quantification of editing',os.path.basename(quant_of_editing_freq_filename))]
        ###############################################################################################################################################


        for ref_name in ref_names:
            ref_len = refs[ref_name]['sequence_length']
            ref_seq = refs[ref_name]['sequence']
            xmin = refs[ref_name]['xmin']
            xmax = refs[ref_name]['xmax']
            min_cut = refs[ref_name]['min_cut']
            max_cut = refs[ref_name]['max_cut']
            hdensity = refs[ref_name]['hdensity']
            hlengths = refs[ref_name]['hlengths']
            center_index = refs[ref_name]['center_index']
            include_idxs_list = refs[ref_name]['include_idxs']
            quantification_window_ref_seq = [list(ref_seq)[x] for x in include_idxs_list]
            sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
#            print('debug 2120 here: '+str(sgRNA_sequences))
            cut_points = refs[ref_name]['sgRNA_cut_points']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
            tot_aln_reads = counts_total[ref_name]
            n_this_category = counts_total[ref_name]

            #only show reference name in filenames if more than one reference
            ref_plot_name = ref_name+"."
            if len(ref_names) == 1 and ref_names[0] == "Reference":
                ref_plot_name = ""

            if n_this_category < 1:
                continue

            #plot quilt for this amplicon  (if not crispresso1 mode)
            if not args.crispresso1_mode:
                ##nucleotide counts
                df_nuc_freq = pd.DataFrame([base_count_vectors[ref_name+"_A"],base_count_vectors[ref_name+"_C"],base_count_vectors[ref_name+"_G"],base_count_vectors[ref_name+"_T"],base_count_vectors[ref_name+"_N"],base_count_vectors[ref_name+'_-']])
                df_nuc_freq.index = ['A','C','G','T','N','-']
                df_nuc_freq.columns = quantification_window_ref_seq
                #print table showing nuc frequencies (sum to total alleles) (in quantification window)
                quant_window_nuc_freq_filename = _jp(ref_plot_name + 'Quantification_window_nucleotide_frequency_table.txt')
                df_nuc_freq.to_csv(quant_window_nuc_freq_filename,sep='\t',header=True,index=True)
                crispresso2_info['refs'][ref_name]['quant_window_nuc_freq_filename'] = os.path.basename(quant_window_nuc_freq_filename)

                df_nuc_pct = df_nuc_freq.divide(tot_aln_reads)
                quant_window_nuc_pct_filename = _jp(ref_plot_name + 'Quantification_window_nucleotide_percentage_table.txt')
                df_nuc_pct.to_csv(quant_window_nuc_pct_filename,sep='\t',header=True,index=True)
                crispresso2_info['refs'][ref_name]['quant_window_nuc_pct_filename'] = os.path.basename(quant_window_nuc_pct_filename)

                df_nuc_freq_all = pd.DataFrame([all_base_count_vectors[ref_name+"_A"],all_base_count_vectors[ref_name+"_C"],all_base_count_vectors[ref_name+"_G"],all_base_count_vectors[ref_name+"_T"],all_base_count_vectors[ref_name+"_N"],all_base_count_vectors[ref_name+'_-']])
                df_nuc_freq_all.index = ['A','C','G','T','N','-']
                df_nuc_freq_all.columns = list(ref_seq)
                #print table showing nuc frequencies (sum to total alleles) (in entire region)
                nuc_freq_filename = _jp(ref_plot_name + 'Nucleotide_frequency_table.txt')
                df_nuc_freq_all.to_csv(nuc_freq_filename,sep='\t',header=True,index=True)
                crispresso2_info['refs'][ref_name]['nuc_freq_filename'] = os.path.basename(nuc_freq_filename)

                df_nuc_pct_all = df_nuc_freq_all.divide(tot_aln_reads)
                nuc_pct_filename = _jp(ref_plot_name + 'Nucleotide_percentage_table.txt')
                df_nuc_pct_all.to_csv(nuc_pct_filename,sep='\t',header=True,index=True)
                crispresso2_info['refs'][ref_name]['nuc_pct_filename'] = os.path.basename(nuc_pct_filename)

                #substitution frequencies
                df_sub_freq,alt_nuc_counts = count_alternate_alleles(
                    sub_base_vectors = substitution_base_vectors,
                    ref_name = ref_name,
                    ref_sequence = quantification_window_ref_seq,
                    ref_total_aln_reads = tot_aln_reads
                    )

                #print table showing sub frequencies
                quant_window_sub_freq_filename =_jp(ref_plot_name + 'Quantification_window_substitution_frequency_table.txt')
                df_sub_freq.to_csv(quant_window_sub_freq_filename,sep='\t',header=True,index=True)
                crispresso2_info['refs'][ref_name]['quant_window_sub_freq_filename'] = os.path.basename(quant_window_sub_freq_filename)


                df_sub_freq_all,alt_nuc_counts_all = count_alternate_alleles(
                    sub_base_vectors = all_substitution_base_vectors,
                    ref_name = ref_name,
                    ref_sequence = ref_seq,
                    ref_total_aln_reads = tot_aln_reads
                    )

                sub_freq_table_filename = _jp(ref_plot_name + 'Substitution_frequency_table.txt')
                df_sub_freq_all.to_csv(sub_freq_table_filename,sep='\t',header=True,index=True)
                crispresso2_info['refs'][ref_name]['sub_freq_table_filename'] = os.path.basename(sub_freq_table_filename)

                if not args.suppress_plots:


                    mod_pcts = []
                    tot = float(counts_total[ref_name])
                    mod_pcts.append(np.concatenate((['Insertions'], np.array(all_insertion_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Insertions_Left'], np.array(all_insertion_left_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Deletions'], np.array(all_deletion_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Substitutions'], np.array(all_substitution_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['All_modifications'], np.array(all_indelsub_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Total'],[counts_total[ref_name]]*refs[ref_name]['sequence_length'])))
                    colnames = ['Modification']+list(ref_seq)
                    modification_percentage_summary_df = pd.DataFrame(mod_pcts,columns=colnames)

                    nuc_df_for_plot = df_nuc_pct_all.reset_index().rename(columns={'index':'Nucleotide'})
                    nuc_df_for_plot.insert(0,'Batch',ref_name) #this function was designed for plottin batch... so just add a column in there to make it happy
                    mod_df_for_plot = modification_percentage_summary_df.copy()
                    mod_df_for_plot.insert(0,'Batch',ref_name)

                    plot_root = _jp('2a.'+ref_plot_name + 'Nucleotide_percentage_quilt')
                    CRISPRessoPlot.plot_nucleotide_quilt(nuc_df_for_plot,mod_df_for_plot,plot_root,save_png,sgRNA_intervals=sgRNA_intervals,quantification_window_idxs=include_idxs_list)
                    crispresso2_info['refs'][ref_name]['plot_2a_root'] = os.path.basename(plot_root)
                    crispresso2_info['refs'][ref_name]['plot_2a_caption'] = "Figure 2a: Nucleotide distribution across amplicon. At each base in the reference amplicon, the percentage of each base as observed in sequencing reads is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the percentage of reads for which that base was deleted. Brown bars between bases show the percentage of reads having an insertion at that position."
                    crispresso2_info['refs'][ref_name]['plot_2a_data'] = [('Nucleotide frequency table',os.path.basename(nuc_freq_filename))]

                    crispresso2_info['refs'][ref_name]['plot_2b_roots'] = []
                    crispresso2_info['refs'][ref_name]['plot_2b_captions'] = []
                    crispresso2_info['refs'][ref_name]['plot_2b_datas'] = []
                    for sgRNA,cut_point in zip(sgRNA_sequences,cut_points):
                        #get nucleotide columns to print for this sgRNA
                        sel_cols = [0,1]
                        plot_half_window = max(1,args.plot_window_size)
                        new_sel_cols_start = max(2,cut_point-plot_half_window+1)
                        new_sel_cols_end = min(ref_len,cut_point+plot_half_window+1)
                        sel_cols.extend(list(range(new_sel_cols_start+2,new_sel_cols_end+2)))
                        #get new intervals
                        new_sgRNA_intervals = []
                        #add annotations for each sgRNA (to be plotted on this sgRNA's plot)
                        for (int_start,int_end) in refs[ref_name]['sgRNA_intervals']:
                            new_sgRNA_intervals += [(int_start - new_sel_cols_start,int_end - new_sel_cols_start)]
                        new_include_idx = []
                        for x in include_idxs_list:
                            new_include_idx += [x - new_sel_cols_start]
                        plot_root = _jp('2b.'+ref_plot_name + 'Nucleotide_percentage_quilt_around_sgRNA_' + sgRNA)
                        CRISPRessoPlot.plot_nucleotide_quilt(
                                nuc_df_for_plot.iloc[:,sel_cols],
                                mod_df_for_plot.iloc[:,sel_cols],
                                plot_root,
                                save_png,sgRNA_intervals=new_sgRNA_intervals,quantification_window_idxs=new_include_idx)
                        crispresso2_info['refs'][ref_name]['plot_2b_roots'].append(os.path.basename(plot_root))
                        crispresso2_info['refs'][ref_name]['plot_2b_captions'].append('Figure 2b: Nucleotide distribution around sgRNA ' + sgRNA + '.')
                        crispresso2_info['refs'][ref_name]['plot_2b_datas'].append([('Nucleotide frequency in quantification window',os.path.basename(quant_window_nuc_freq_filename))])


            ###############################################################################################################################################
            #(3)visualize effective lengths of reads aligning to this amplicon

            if not args.suppress_plots:

                n_this_category = counts_total[ref_name]
                if n_this_category < 1:
                    continue

                plt.figure(figsize=(10,10))
                ax = plt.subplot(111)
                densityPct = 0.0
                densityPcts = [0.0]*len(hdensity)
                if hdensity.sum() > 0:
                    densityPct = hdensity[center_index]/(float(hdensity.sum()))*100.0
                    densityPcts = hdensity/(float(hdensity.sum()))*100.0
                ax.bar(0,densityPct,color='red',linewidth=0)
                #plt.hold(True)
                barlist=plt.bar(hlengths,densityPcts,align='center',linewidth=0)
                barlist[center_index].set_color('r')
                ax.set_xlim([xmin,xmax])
                ax.set_title(get_plot_title_with_ref_name('Indel size distribution',ref_name))
                ax.set_ylabel('Sequences % (no.)')
                y_label_values= np.round(np.linspace(0, min(100,max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                ax.set_yticks(y_label_values)
                ax.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*n_this_category) for pct in y_label_values])
                ax.set_xlabel('Indel size (bp)')
                #lgd=plt.legend(['No indel','Indel'])
                lgd=ax.legend(['No indel','Indel'],loc='center', bbox_to_anchor=(0.5, -0.18),ncol=1, fancybox=True, shadow=True)
                lgd.legendHandles[0].set_height(3)
                lgd.legendHandles[1].set_height(3)

                plot_root = _jp('3a.'+ref_plot_name+'Indel_size_distribution')
                plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_inches='tight')
                plt.close()

                crispresso2_info['refs'][ref_name]['plot_3a_root'] = os.path.basename(plot_root)
                crispresso2_info['refs'][ref_name]['plot_3a_caption'] = "Figure 3a: Frequency distribution of alleles with indels (blue) and without indels (red)."
                crispresso2_info['refs'][ref_name]['plot_3a_data'] = [('Indel histogram',os.path.basename(crispresso2_info['refs'][ref_name]['indel_histogram_filename']))]
                ###############################################################################################################################################

                ###############################################################################################################################################
                #(3b) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);

                y_values_mut = refs[ref_name]['y_values_mut']
                x_bins_mut = refs[ref_name]['x_bins_mut']
                y_values_ins = refs[ref_name]['y_values_ins']
                x_bins_ins = refs[ref_name]['x_bins_ins']
                y_values_del = refs[ref_name]['y_values_del']
                x_bins_del = refs[ref_name]['x_bins_del']


                fig=plt.figure(figsize=(26,6.5))

                ax=fig.add_subplot(1,3,1)
                ax.bar(x_bins_ins[:-1],y_values_ins,align='center',linewidth=0,color=(0,0,1))
                barlist=ax.bar(x_bins_ins[:-1],y_values_ins,align='center',linewidth=0,color=(0,0,1))
                barlist[0].set_color('r')
                plt.title(get_plot_title_with_ref_name('Insertions',ref_name))
                plt.xlabel('Size (bp)')
                plt.ylabel('Sequences % (no.)')
                lgd=plt.legend(['Non-insertion','Insertion'][::-1], bbox_to_anchor=(.82, -0.22),ncol=1, fancybox=True, shadow=True)
                lgd.legendHandles[0].set_height(6)
                lgd.legendHandles[1].set_height(6)
                plt.xlim(xmin=-1)
                y_label_values= np.round(np.linspace(0, min(counts_total[ref_name],max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/counts_total[ref_name]*100,n_reads) for n_reads in y_label_values])

                ax=fig.add_subplot(1,3,2)
                ax.bar(-x_bins_del[:-1],y_values_del,align='center',linewidth=0,color=(0,0,1))
                barlist=ax.bar(-x_bins_del[:-1],y_values_del,align='center',linewidth=0,color=(0,0,1))
                barlist[0].set_color('r')
                plt.title(get_plot_title_with_ref_name('Deletions',ref_name))
                plt.xlabel('Size (bp)')
                plt.ylabel('Sequences % (no.)')
                lgd=plt.legend(['Non-deletion','Deletion'][::-1], bbox_to_anchor=(.82, -0.22),ncol=1, fancybox=True, shadow=True)
                lgd.legendHandles[0].set_height(6)
                lgd.legendHandles[1].set_height(6)
                plt.xlim(xmax=1)
                y_label_values= np.round(np.linspace(0, min(counts_total[ref_name],max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/counts_total[ref_name]*100,n_reads) for n_reads in y_label_values])



                ax=fig.add_subplot(1,3,3)
                ax.bar(x_bins_mut[:-1],y_values_mut,align='center',linewidth=0,color=(0,0,1))
                barlist=ax.bar(x_bins_mut[:-1],y_values_mut,align='center',linewidth=0,color=(0,0,1))
                barlist[0].set_color('r')
                plt.title(get_plot_title_with_ref_name('Substitutions',ref_name))
                plt.xlabel('Positions substituted (number)')
                plt.ylabel('Sequences % (no.)')
                lgd=plt.legend(['Non-substitution','Substitution'][::-1] ,bbox_to_anchor=(.82, -0.22),ncol=1, fancybox=True, shadow=True)
                lgd.legendHandles[0].set_height(6)
                lgd.legendHandles[1].set_height(6)
                plt.xlim(xmin=-1)
                y_label_values= np.round(np.linspace(0, min(counts_total[ref_name],max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/counts_total[ref_name]*100,n_reads) for n_reads in y_label_values])


                plt.tight_layout()

                plot_root = _jp('3b.'+ref_plot_name+'Insertion_deletion_substitutions_size_hist')
                plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_inches='tight')
                plt.close()

                crispresso2_info['refs'][ref_name]['plot_3b_root'] = os.path.basename(plot_root)
                crispresso2_info['refs'][ref_name]['plot_3b_caption'] = "Figure 3b: Left panel, frequency distribution of sequence modifications that increase read length with respect to the reference amplicon, classified as insertions (positive indel size). Middle panel, frequency distribution of sequence modifications that reduce read length with respect to the reference amplicon, classified as deletions (negative indel size). Right panel, frequency distribution of sequence modifications that do not alter read length with respect to the reference amplicon, which are classified as substitutions (number of substituted positions shown)."
                crispresso2_info['refs'][ref_name]['plot_3b_data'] = []


                #(4) another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);
                #Indels location Plots

                n_this_category_modified = 0
                modifiedName = ref_name + "_MODIFIED"
                if modifiedName in class_counts:
                    n_this_category_modified = class_counts[modifiedName]

                plt.figure(figsize=(10,10))

                y_max=max(all_indelsub_count_vectors[ref_name])*1.1

                #shade quantification window
                if len(include_idxs_list) > 1:
                    lastStart = include_idxs_list[0]
                    lastIdx = include_idxs_list[0]
                    for idx in range(1,len(include_idxs_list)):
                        if include_idxs_list[idx] == lastIdx + 1:
                            lastIdx = include_idxs_list[idx]
                        else:
                            p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                            plt.gca().add_patch(p) #gca = get current axis
                            lastStart = include_idxs_list[idx]
                            lastIdx = include_idxs_list[idx]
                    p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
                    plt.gca().add_patch(p)

                plt.plot(all_indelsub_count_vectors[ref_name],'r',lw=3,label=get_plot_title_with_ref_name('Combined Insertions/Deletions/Substitutions',ref_name))
                 #plt.hold(True)

                if cut_points:
                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                            plt.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                        else:
                            plt.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',lw=2,label='_nolegend_')

                    for idx,sgRNA_int in enumerate(sgRNA_intervals):
                        if idx==0:
                            plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                        else:
                            plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


                lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
                ylabel_values = np.arange(0,1,1.0/6.0)
                if y_max > 0:
                    y_label_values=np.arange(0,y_max,y_max/6.0)
                if len(ref_names) == 1:
                    plt.ylabel('Sequences: % Total ( no. )')
                    plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])
                else:
                    plt.ylabel('Sequences: % Total ( % '+ref_name+', no. )')
                    plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])
                plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )

                plt.title(get_plot_title_with_ref_name('Mutation position distribution',ref_name))
                plt.xlabel('Reference amplicon position (bp)')
                plt.ylim(0,max(1,y_max))
                plt.xlim(0,ref_len-1)

                plot_root = _jp('4a.'+ref_plot_name+'Combined_insertion_deletion_substitution_locations')
                plt.savefig(plot_root+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
                plt.close()

                crispresso2_info['refs'][ref_name]['plot_4a_root'] = os.path.basename(plot_root)
                crispresso2_info['refs'][ref_name]['plot_4a_caption'] = "Figure 4a: Combined frequency of any modification across the amplicon. Modifications outside of the quantification window are also shown."
                crispresso2_info['refs'][ref_name]['plot_4a_data'] = []

    #            print("subs: " + ref_name + ":"+ str(all_substitution_count_vectors[ref_name]))
                plt.figure(figsize=(10,10))

                #shade quantification window
                if len(include_idxs_list) > 1:
                    lastStart = include_idxs_list[0]
                    lastIdx = include_idxs_list[0]
                    for idx in range(1,len(include_idxs_list)):
                        if include_idxs_list[idx] == lastIdx + 1:
                            lastIdx = include_idxs_list[idx]
                        else:
                            p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                            plt.gca().add_patch(p) #gca = get current axis
                            lastStart = include_idxs_list[idx]
                            lastIdx = include_idxs_list[idx]
                    p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
                    plt.gca().add_patch(p)

                plt.plot(all_insertion_count_vectors[ref_name],'r',lw=3,label='Insertions')
                #plt.hold(True)
                plt.plot(all_deletion_count_vectors[ref_name],'m',lw=3,label='Deletions')
                plt.plot(all_substitution_count_vectors[ref_name],'g',lw=3,label='Substitutions')

                y_max=max(max(all_insertion_count_vectors[ref_name]),max(all_deletion_count_vectors[ref_name]),max(all_substitution_count_vectors[ref_name]))*1.1


                if cut_points:
                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                            plt.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                        else:
                            plt.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',lw=2,label='_nolegend_')

                    for idx,sgRNA_int in enumerate(sgRNA_intervals):
                        if idx==0:
                            plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                        else:
                            plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.31),ncol=1, fancybox=True, shadow=True)
                y_label_values = np.arange(0,1,1.0/6.0)
                if y_max > 0:
                    y_label_values=np.arange(0,y_max,y_max/6.0)
                plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )

                plt.xlabel('Reference amplicon position (bp)')
                if len(ref_names) == 1:
                    plt.ylabel('Sequences: % Total ( no. )')
                    #plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(N_MODIFIED)*100, n_reads) for n_reads in y_label_values])
                    plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])
                else:
                    plt.ylabel('Sequences: % Total ( % '+ref_name+', no. )')
                    plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])

                plt.ylim(0,max(1,y_max))
                plt.xlim(0,ref_len-1)

                plt.title(get_plot_title_with_ref_name('Mutation position distribution',ref_name))

                plot_root = _jp('4b.'+ref_plot_name+'Insertion_deletion_substitution_locations')
                plt.savefig(plot_root+'.pdf',bbox_extra_artists=(lgd,),pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
                plt.close()
                crispresso2_info['refs'][ref_name]['plot_4b_root'] = os.path.basename(plot_root)
                crispresso2_info['refs'][ref_name]['plot_4b_caption'] = "Figure 4b: Frequency of insertions (red), deletions (purple), and substitutions (green) across the entire amplicon, including modifications outside of the quantification window."
                crispresso2_info['refs'][ref_name]['plot_4b_data'] = [('Modification frequency',os.path.basename(mod_count_filename))]


                plt.figure(figsize=(10,10))

                y_max=max(max(insertion_count_vectors[ref_name]),max(deletion_count_vectors[ref_name]),max(substitution_count_vectors[ref_name]),1)*1.1

                #shade quantification window
                if len(include_idxs_list) > 1:
                    lastStart = include_idxs_list[0]
                    lastIdx = include_idxs_list[0]
                    for idx in range(1,len(include_idxs_list)):
                        if include_idxs_list[idx] == lastIdx + 1:
                            lastIdx = include_idxs_list[idx]
                        else:
                            p = matplotlib.patches.Rectangle((lastStart -0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                            plt.gca().add_patch(p) #gca = get current axis
                            lastStart = include_idxs_list[idx]
                            lastIdx = include_idxs_list[idx]
                    p = matplotlib.patches.Rectangle((lastStart -0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
                    plt.gca().add_patch(p)



                plt.plot(insertion_count_vectors[ref_name],'r',linewidth=3,label='Insertions')
                plt.plot(deletion_count_vectors[ref_name],'m',linewidth=3,label='Deletions')
                plt.plot(substitution_count_vectors[ref_name],'g',linewidth=3,label='Substitutions')


                if cut_points:
                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                            plt.plot([cut_point+.5,cut_point+.5],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                        else:
                            plt.plot([cut_point+.5,cut_point+.5],[0,y_max],'--k',linewidth=2,label='_nolegend_')


                    for idx,sgRNA_int in enumerate(sgRNA_intervals):
                        if idx==0:
                            plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                        else:
                            plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


                lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.31),ncol=1, fancybox=True, shadow=True)
                y_label_values = np.arange(0,1,1.0/6.0)
                if y_max > 0:
                    y_label_values=np.arange(0,y_max,y_max/6.0).astype(int)
                plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])
                plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )

                plt.xlabel('Reference amplicon position (bp)')
                if len(ref_names) == 1:
                    plt.ylabel('Sequences: % Total ( no. )')
                    plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])
                else:
                    plt.ylabel('Sequences: % Total ( % '+ref_name+', no. )')
                    plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])


                plt.ylim(0,max(1,y_max))
                plt.xlim(0,ref_len-1)
                plt.title(get_plot_title_with_ref_name('Mutation position distribution',ref_name))
                plot_root = _jp('4c.'+ref_plot_name+'Quantification_window_insertion_deletion_substitution_locations')
                plt.savefig(plot_root+'.pdf',bbox_extra_artists=(lgd,),pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
                plt.close()

                crispresso2_info['refs'][ref_name]['plot_4c_root'] = os.path.basename(plot_root)
                crispresso2_info['refs'][ref_name]['plot_4c_caption'] = "Figure 4c: Frequency of insertions (red), deletions (purple), and substitutions (green) across the entire amplicon, considering only modifications that overlap with the quantification window."
                crispresso2_info['refs'][ref_name]['plot_4c_data'] = [('Modification frequency in quantification window',os.path.basename(quant_window_mod_count_filename))]

                #Position dependent indels plot
                fig=plt.figure(figsize=(24,10))
                ax1=fig.add_subplot(1,2,1)
                markerline, stemlines, baseline=ax1.stem(insertion_length_vectors[ref_name])
                plt.setp(markerline, 'markerfacecolor', 'r', 'markersize', 8)
                plt.setp(baseline, 'linewidth', 0)
                plt.setp(stemlines, 'color', 'r','linewidth',3)
                #plt.hold(True)
                y_max=max(insertion_length_vectors[ref_name])*1.1
                if cut_points:

                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                            ax1.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                        else:
                            ax1.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='_nolegend_')

                plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )
                plt.xlabel('Reference amplicon position (bp)')
                plt.ylabel('Average insertion length')
                plt.ylim(0,max(1,y_max))
                plt.xlim(0,ref_len-1)
                ax1.set_title(get_plot_title_with_ref_name('Position dependent insertion size',ref_name))
                plt.tight_layout()

                ax2=fig.add_subplot(1,2,2)
                markerline, stemlines, baseline=ax2.stem(deletion_length_vectors[ref_name])
                plt.setp(markerline, 'markerfacecolor', 'm', 'markersize', 8)
                plt.setp(baseline, 'linewidth', 0)
                plt.setp(stemlines, 'color', 'm','linewidth',3)
                #plt.hold(True)
                y_max=max(deletion_length_vectors[ref_name])*1.1
                if cut_points:

                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                            ax2.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                        else:
                            ax2.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='_nolegend_')

                plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )
                plt.xlabel('Reference amplicon position (bp)')
                plt.ylabel('Average deletion length')

                ymin, ymax = ax2.yaxis.get_view_interval()
                plt.ylim(ymin=0,ymax=max(1,y_max))
                plt.xlim(0,ref_len-1)
                ax2.set_title(get_plot_title_with_ref_name('Position dependent deletion size', ref_name))

                plt.tight_layout()


                plot_root = _jp('4d.'+ref_plot_name+'Position_dependent_average_indel_size')
                plt.savefig(plot_root+'.pdf',bbox_extra_artists=(lgd,),pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
                plt.close()

                crispresso2_info['refs'][ref_name]['plot_4d_root'] = os.path.basename(plot_root)
                crispresso2_info['refs'][ref_name]['plot_4d_caption'] = "Figure 4d: Position dependent insertion size(left) and deletion size (right), including only modifications that overlap with the quantification window."
                crispresso2_info['refs'][ref_name]['plot_4d_data'] = []

                ###############################################################################################################################################
                #4e : for HDR, global modifications with respect to reference 1
                ###############################################################################################################################################

                if args.expected_hdr_amplicon_seq != "" and (ref_name == ref_names[0] or ref_name == "HDR"):
                    plt.figure(figsize=(10,10))
                    ref1_all_insertion_positions = ref1_all_insertion_count_vectors[ref_name]
                    ref1_all_deletion_positions = ref1_all_deletion_count_vectors[ref_name]
                    ref1_all_substitution_positions = ref1_all_substitution_count_vectors[ref_name]

                    y_max=max(max(ref1_all_insertion_positions),max(ref1_all_deletion_positions),max(ref1_all_substitution_positions),1)*1.1

                    plt.plot(ref1_all_insertion_positions,'r',linewidth=3,label='Insertions')
                    plt.plot(ref1_all_deletion_positions,'m',linewidth=3,label='Deletions')
                    plt.plot(ref1_all_substitution_positions,'g',linewidth=3,label='Substitutions')


                    ref1_cut_points = refs[ref_names[0]]['sgRNA_cut_points']
                    ref1_sgRNA_intervals = refs[ref_names[0]]['sgRNA_intervals']
                    ref1_include_idxs_list = sorted(list(refs[ref_names[0]]['include_idxs']))
                    #shade quantification window
                    if len(include_idxs_list) > 1:
                        lastStart = include_idxs_list[0]
                        lastIdx = include_idxs_list[0]
                        for idx in range(1,len(include_idxs_list)):
                            if include_idxs_list[idx] == lastIdx + 1:
                                lastIdx = include_idxs_list[idx]
                            else:
                                p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                                plt.gca().add_patch(p) #gca = get current axis
                                lastStart = include_idxs_list[idx]
                                lastIdx = include_idxs_list[idx]
                        p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
                        plt.gca().add_patch(p)

                    if ref1_cut_points:
                        for idx,ref1_cut_point in enumerate(ref1_cut_points):
                            if idx==0:
                                plt.plot([ref1_cut_point+0.5,ref1_cut_point+0.5],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                            else:
                                plt.plot([ref1_cut_point+0.5,ref1_cut_point+0.5],[0,y_max],'--k',linewidth=2,label='_nolegend_')


                        for idx,sgRNA_int in enumerate(ref1_sgRNA_intervals):
                            if idx==0:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                            else:
                                plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


                    lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.31),ncol=1, fancybox=True, shadow=True)
                    y_label_values = np.arange(0,1,1.0/6.0)
                    if y_max > 0:
                        y_label_values=np.arange(0,y_max,y_max/6.0).astype(int)

                    plt.ylabel('Sequences: % Total ( no. )')
                    plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])

                    plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )
                    plt.xlabel(ref_names[0] + ' position (bp)')

                    plt.ylim(0,max(1,y_max))
                    plt.xlim(0,ref_len-1)
                    if ref_name == ref_names[0]:
                        plt.title('Mutation position distribution in all reads with reference to %s'%(ref_names[0]))
                        plot_root = _jp('4e.' + ref_names[0] + '.Global_mutations_in_all_reads')
                        crispresso2_info['refs'][ref_names[0]]['plot_4e_root'] = os.path.basename(plot_root)
                        crispresso2_info['refs'][ref_names[0]]['plot_4e_caption'] = "Figure 4e: Modifications in all reads when aligned to the reference sequence. Insertions: red, deletions: purple, substitutions: green. All modifications (including those outside the quantification window) are shown."
                        crispresso2_info['refs'][ref_names[0]]['plot_4e_data'] = []
                    elif ref_name == "HDR":
                        plt.title('Mutation position distribution in %s reads with reference to %s'%(ref_name,ref_names[0]))
                        plot_root = _jp('4f.' + ref_names[0] + '.Global_mutations_in_HDR_reads_with_reference_to_'+ref_names[0])
                        crispresso2_info['refs'][ref_names[0]]['plot_4f_root'] = os.path.basename(plot_root)
                        crispresso2_info['refs'][ref_names[0]]['plot_4f_caption'] = "Figure 4f: Modifications in HDR reads with respect to the reference sequence. Insertions: red, deletions: purple, substitutions: green. All modifications (including those outside the quantification window) are shown."
                        crispresso2_info['refs'][ref_names[0]]['plot_4f_data'] = []

                    plt.savefig(plot_root+'.pdf',bbox_extra_artists=(lgd,),pad_inches=1,bbox_inches='tight')
                    if save_png:
                        plt.savefig(plot_root+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
                    plt.close()


            ###############################################################################################################################################
            #(5, 6) frameshift analyses plots
            if (refs[ref_name]['contains_coding_seq']): #PERFORM FRAMESHIFT ANALYSIS
                #make frameshift plots
                ref_len = refs[ref_name]['sequence_length']
                cut_points = refs[ref_name]['sgRNA_cut_points']
                sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
                MODIFIED_FRAMESHIFT = counts_modified_frameshift[ref_name]
                MODIFIED_NON_FRAMESHIFT = counts_modified_non_frameshift[ref_name]
                NON_MODIFIED_NON_FRAMESHIFT = counts_non_modified_non_frameshift[ref_name]
                SPLICING_SITES_MODIFIED = counts_splicing_sites_modified[ref_name]
                exon_intervals = refs[ref_name]['exon_intervals']
                exon_len_mods = refs[ref_name]['exon_len_mods']
                count_total = counts_total[ref_name]
                count_modified = counts_modified[ref_name]
                count_unmodified = counts_unmodified[ref_name]

                if not args.suppress_plots:
                    fig=plt.figure(figsize=(12*1.5,14.5*1.5))
                    ax1 = plt.subplot2grid((6,3), (0, 0), colspan=3, rowspan=5)

                    patches, texts, autotexts =ax1.pie([MODIFIED_FRAMESHIFT,\
                                                        MODIFIED_NON_FRAMESHIFT,\
                                                        NON_MODIFIED_NON_FRAMESHIFT],\
                                                       labels=['Frameshift mutation\n(%d reads)' %MODIFIED_FRAMESHIFT,\
                                                              'In-frame mutation\n(%d reads)' % MODIFIED_NON_FRAMESHIFT,\
                                                              'Noncoding mutation\n(%d reads)' %NON_MODIFIED_NON_FRAMESHIFT],\
                                                       explode=(0.0,0.0,0.0),\
                                                       colors=[(0.89019608,  0.29019608,  0.2, 0.8),(0.99215686,  0.73333333,  0.51764706,0.8),(0.99607843,  0.90980392,  0.78431373,0.8)],\
                                                       autopct='%1.1f%%')

                    ax2 = plt.subplot2grid((6,3), (5, 0), colspan=3, rowspan=1)
                    ax2.plot([0,ref_len],[0,0],'-k',linewidth=2,label=ref_name+' sequence')
                    #plt.hold(True)

                    for idx,exon_interval in enumerate(exon_intervals):
                        if idx==0:
                            ax2.plot(exon_interval,[0,0],'-',linewidth=10,c=(0,0,1,0.5),label='Coding sequence/s',solid_capstyle='butt')
                        else:
                            ax2.plot(exon_interval,[0,0],'-',linewidth=10,c=(0,0,1,0.5),label='_nolegend_',solid_capstyle='butt')

                    if cut_points:
#                        print('cut poitns are ' + str(cut_points))
#                        ax2.plot(cut_points+1,np.zeros(len(cut_points)),'vr', ms=25,label='Predicted cleavage position')

                        for idx,cut_point in enumerate(cut_points):
                            if not args.base_editor_output:
                                if idx==0:
                                        ax2.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                                else:
                                        ax2.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='_nolegend_')

                            for idx,sgRNA_int in enumerate(sgRNA_intervals):
                                if idx==0:
                                   ax2.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                                else:
                                   ax2.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                    plt.legend(bbox_to_anchor=(0, 0, 1., 0),  ncol=1, mode="expand", borderaxespad=0.,numpoints=1)
                    plt.xlim(0,ref_len)
                    plt.axis('off')

                    proptease = fm.FontProperties()
                    proptease.set_size('xx-large')
                    plt.setp(autotexts, fontproperties=proptease)
                    plt.setp(texts, fontproperties=proptease)

                    plot_root = _jp('5.'+ref_plot_name+'Frameshift_in-frame_mutations_pie_chart')
                    plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                    if save_png:
                        plt.savefig(plot_root+'.png',bbox_inches='tight')
                    plt.close()
                    crispresso2_info['refs'][ref_name]['plot_5_root'] = os.path.basename(plot_root)
                    crispresso2_info['refs'][ref_name]['plot_5_caption'] = "Figure 5: Frameshift analysis of coding sequence reads affected by modifications (unmodified reads are excluded from this analysis)."
                    crispresso2_info['refs'][ref_name]['plot_5_data'] = []



                     #profiles-----------------------------------------------------------------------------------
                    fig=plt.figure(figsize=(22,10))
                    ax1=fig.add_subplot(2,1,1)
                    x,y=map(np.array,zip(*[a for a in hists_frameshift[ref_name].iteritems()]))
                    if sum(hists_frameshift[ref_name].values()) != 0:
                        y=y/float(sum(hists_frameshift[ref_name].values()))*100
                    ax1.bar(x-0.1,y)
                    ax1.set_xlim(-30.5,30.5)
                    ax1.set_frame_on(False)
                    ax1.set_xticks([idx for idx in range(-30,31) if idx % 3])
                    ax1.tick_params(which='both',      # both major and minor ticks are affected
                       bottom=False,      # ticks along the bottom edge are off
                       top=False,         # ticks along the top edge are off
                       labelbottom=True) # labels along the bottom edge are off)
                    ax1.yaxis.tick_left()
                    xmin, xmax = ax1.get_xaxis().get_view_interval()
                    ymin, ymax = ax1.get_yaxis().get_view_interval()
                    ax1.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if idx % 3]],rotation='vertical')
                    plt.title(get_plot_title_with_ref_name('Frameshift profile',ref_name))
                    ax1.tick_params(axis='both', which='major', labelsize=24)
                    ax1.tick_params(axis='both', which='minor', labelsize=24)
                    plt.tight_layout()
                    ax1.set_ylabel('Sequences % (no.)')
                    y_label_values= np.round(np.linspace(0, min(100,max(ax1.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                    ax1.set_yticks(y_label_values)
                    ax1.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*sum(hists_frameshift[ref_name].values())) for pct in y_label_values])

                    ax2=fig.add_subplot(2,1,2)
                    x,y=map(np.array,zip(*[a for a in hists_inframe[ref_name].iteritems()]))
                    if sum(hists_inframe[ref_name].values()) > 0:
                        y=y/float(sum(hists_inframe[ref_name].values()))*100
                    #ax2.bar(x-0.5,y,color=(0,1,1,0.2))
                    ax2.bar(x-0.1,y,color=(0,1,1,0.2))
                    ax2.set_xlim(-30.5,30.5)
                    ax2.set_frame_on(False)
                    ax2.set_xticks([idx for idx in range(-30,31) if (idx % 3 ==0) ])
                    ax2.tick_params(which='both',      # both major and minor ticks are affected
                       bottom=False,      # ticks along the bottom edge are off
                       top=False,         # ticks along the top edge are off
                       labelbottom=True) # labels along the bottom edge are off)
                    ax2.yaxis.tick_left()
                    xmin, xmax = ax2.xaxis.get_view_interval()
                    ymin, ymax = ax2.yaxis.get_view_interval()
                    ax2.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if (idx % 3==0)]],rotation='vertical',horizontalalignment="center")
                    plt.title(get_plot_title_with_ref_name('In-frame profile',ref_name))
                    ax2.set_ylabel('Sequences % (no.)')
                    y_label_values= np.round(np.linspace(0, min(100,max(ax2.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                    ax2.set_yticks(y_label_values)
                    ax2.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*sum(hists_inframe[ref_name].values())) for pct in y_label_values])

                    ax2.tick_params(axis='both', which='major', labelsize=24)
                    ax2.tick_params(axis='both', which='minor', labelsize=24)
                    plt.tight_layout()

                    plot_root = _jp('6.'+ref_plot_name+'Frameshift_in-frame_mutation_profiles')
                    plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                    if save_png:
                        plt.savefig(plot_root+'.png',bbox_inches='tight')
                    plt.close()
                    crispresso2_info['refs'][ref_name]['plot_6_root'] = os.path.basename(plot_root)
                    crispresso2_info['refs'][ref_name]['plot_6_caption'] = "Figure 6: Frameshift and in-frame mutagenesis profiles indicating position affected by modification."
                    crispresso2_info['refs'][ref_name]['plot_6_data'] = []


                     #-----------------------------------------------------------------------------------------------------------
                    fig=plt.figure(figsize=(12*1.5,12*1.5))
                    ax=fig.add_subplot(1,1,1)
                    patches, texts, autotexts =ax.pie([SPLICING_SITES_MODIFIED,\
                                                      (count_total - SPLICING_SITES_MODIFIED)],\
                                                      labels=['Potential splice sites modified\n(%d reads)' %SPLICING_SITES_MODIFIED,\
                                                              'Unmodified\n(%d reads)' % (count_total - SPLICING_SITES_MODIFIED)],\
                                                      explode=(0.0,0),\
                                                      colors=[(0.89019608,  0.29019608,  0.2, 0.8),(0.99607843,  0.90980392,  0.78431373,0.8)],\
                                                      autopct='%1.1f%%')
                    proptease = fm.FontProperties()
                    proptease.set_size('xx-large')
                    plt.setp(autotexts, fontproperties=proptease)
                    plt.setp(texts, fontproperties=proptease)
                    plot_root = _jp('8.'+ref_plot_name+'Potential_splice_sites_pie_chart')
                    plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                    if save_png:
                        plt.savefig(plot_root+'.png',bbox_inches='tight')
                    plt.close()
                    crispresso2_info['refs'][ref_name]['plot_8_root'] = os.path.basename(plot_root)
                    crispresso2_info['refs'][ref_name]['plot_8_caption'] = "Figure 8: Predicted impact on splice sites. Potential splice sites modified refers to reads in which the either of the two intronic positions adjacent to exon junctions are disrupted."
                    crispresso2_info['refs'][ref_name]['plot_8_data'] = []


                    #non coding
                    plt.figure(figsize=(10,10))
                    plt.plot(insertion_count_vectors_noncoding[ref_name],'r',linewidth=3,label='Insertions')
                    #plt.hold(True)
                    plt.plot(deletion_count_vectors_noncoding[ref_name],'m',linewidth=3,label='Deletions')
                    plt.plot(substitution_count_vectors_noncoding[ref_name],'g',linewidth=3,label='Substitutions')

                    y_max=max(max(insertion_count_vectors_noncoding[ref_name]),max(deletion_count_vectors_noncoding[ref_name]),max(substitution_count_vectors_noncoding[ref_name]))*1.1

                    #shade quantification window
                    if len(include_idxs_list) > 1:
                        lastStart = include_idxs_list[0]
                        lastIdx = include_idxs_list[0]
                        for idx in range(1,len(include_idxs_list)):
                            if include_idxs_list[idx] == lastIdx + 1:
                                lastIdx = include_idxs_list[idx]
                            else:
                                p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                                plt.gca().add_patch(p) #gca = get current axis
                                lastStart = include_idxs_list[idx]
                                lastIdx = include_idxs_list[idx]
                        p = matplotlib.patches.Rectangle((lastStart-0.5, 0), 1+(lastIdx-lastStart)-0.5, y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
                        plt.gca().add_patch(p)

                    if cut_points:
                        for idx,cut_point in enumerate(cut_points):
                            if not args.base_editor_output:
                                if idx==0:
                                        plt.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                                else:
                                        plt.plot([cut_point+0.5,cut_point+0.5],[0,y_max],'--k',linewidth=2,label='_nolegend_')

                            for idx,sgRNA_int in enumerate(sgRNA_intervals):
                                if idx==0:
                                   plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                                else:
                                   plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                    lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.31),ncol=1, fancybox=True, shadow=True)
                    plt.xticks(np.arange(0,ref_len,max(3,(ref_len/6) - (ref_len/6)%5)).astype(int) )

                    plt.xlabel('Reference amplicon position (bp)')
                    plt.ylabel('Sequences (no.)')
                    plt.ylim(0,max(1,y_max))
                    plt.xlim(0,ref_len-1)
                    plt.title(get_plot_title_with_ref_name('Noncoding mutation position distribution',ref_name))

                    plot_root = _jp('7.'+ref_plot_name+'Insertion_deletion_substitution_locations_noncoding')
                    plt.savefig(plot_root+'.pdf',bbox_extra_artists=(lgd,),pad_inches=1,bbox_inches='tight')
                    if save_png:
                        plt.savefig(plot_root+'.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
                    plt.close()
                    crispresso2_info['refs'][ref_name]['plot_7_root'] = os.path.basename(plot_root)
                    crispresso2_info['refs'][ref_name]['plot_7_caption'] = "Figure 7: Reads with insertions (red), deletions (purple), and substitutions (green) mapped to reference amplicon position exclusively in noncoding region/s (that is, without mutations affecting coding sequences). The predicted cleavage site is indicated by a vertical dashed line. Only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are plotted."
                    crispresso2_info['refs'][ref_name]['plot_7_data'] = []
            #end contains coding seq


            ######PLOT
            if not args.crispresso1_mode and args.base_editor_output:
                if not args.suppress_plots:

                    fig_filename_root= _jp('10a.'+ref_plot_name+'Substitution_frequencies_at_each_bp')
                    CRISPRessoPlot.plot_subs_across_ref(
                        ref_len = ref_len,
                        ref_seq = ref_seq,
                        ref_name = ref_name,
                        ref_count = tot_aln_reads,
                        all_substitution_base_vectors = all_substitution_base_vectors,

                        plot_title = get_plot_title_with_ref_name('Substitution frequency', ref_name),
                        fig_filename_root= fig_filename_root,
                        save_also_png = save_png,
                        quantification_window_idxs=include_idxs_list
                        )
                    crispresso2_info['refs'][ref_name]['plot_10a_root'] = os.path.basename(fig_filename_root)
                    crispresso2_info['refs'][ref_name]['plot_10a_caption'] = "Figure 10a: Substitution frequencies across the amplicon."
                    crispresso2_info['refs'][ref_name]['plot_10a_data'] = [('Nucleotide frequencies',os.path.basename(nuc_freq_filename))]


                    #plot all substitution rates in entire region
                    fig_filename_root = _jp('10b.'+ref_plot_name+'Substitution_frequency_barplot')
                    CRISPRessoPlot.plot_sub_freqs(
                        alt_nuc_counts = alt_nuc_counts_all,
                        plot_title = get_plot_title_with_ref_name('Substitution frequency\nin entire amplicon', ref_name),
                        fig_filename_root = fig_filename_root,
                        save_also_png = save_png
                        )
                    crispresso2_info['refs'][ref_name]['plot_10b_root'] = os.path.basename(fig_filename_root)
                    crispresso2_info['refs'][ref_name]['plot_10b_caption'] = "Figure 10b: Substitution frequencies across the amplicon."
                    crispresso2_info['refs'][ref_name]['plot_10b_data'] = [('Nucleotide frequencies',os.path.basename(nuc_freq_filename))]

                    #plot all substitution rates in quantification_window
                    fig_filename_root = _jp('10c.'+ref_plot_name+'Substitution_frequency_barplot_around_sgRNA_'+sgRNA)
                    CRISPRessoPlot.plot_sub_freqs(
                        alt_nuc_counts = alt_nuc_counts,
                        plot_title = get_plot_title_with_ref_name('Substitution frequency\nin quantification window', ref_name),
                        fig_filename_root = fig_filename_root,
                        save_also_png = save_png
                        )
                    crispresso2_info['refs'][ref_name]['plot_10c_root'] = os.path.basename(fig_filename_root)
                    crispresso2_info['refs'][ref_name]['plot_10c_caption'] = "Figure 10c: Substitution frequencies in the quantification window"
                    crispresso2_info['refs'][ref_name]['plot_10c_data'] = [('Nucleotide frequencies in quantification window',os.path.basename(quant_window_sub_freq_filename))]



            ##new plots alleles around cut_sites
            sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
            sgRNA_cut_points = refs[ref_name]['sgRNA_cut_points']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
            sgRNA_plot_idxs = refs[ref_name]['sgRNA_plot_idxs']

            crispresso2_info['refs'][ref_name]['plot_9_roots'] = []
            crispresso2_info['refs'][ref_name]['plot_9_captions'] = []
            crispresso2_info['refs'][ref_name]['plot_9_datas'] = []
            crispresso2_info['refs'][ref_name]['allele_frequency_files'] = []

            crispresso2_info['refs'][ref_name]['plot_10c_roots'] = []
            crispresso2_info['refs'][ref_name]['plot_10c_captions'] = []
            crispresso2_info['refs'][ref_name]['plot_10c_datas'] = []

            crispresso2_info['refs'][ref_name]['plot_10d_roots'] = []
            crispresso2_info['refs'][ref_name]['plot_10d_captions'] = []
            crispresso2_info['refs'][ref_name]['plot_10d_datas'] = []

            crispresso2_info['refs'][ref_name]['plot_10e_roots'] = []
            crispresso2_info['refs'][ref_name]['plot_10e_captions'] = []
            crispresso2_info['refs'][ref_name]['plot_10e_datas'] = []

            crispresso2_info['refs'][ref_name]['plot_10f_roots'] = []
            crispresso2_info['refs'][ref_name]['plot_10f_captions'] = []
            crispresso2_info['refs'][ref_name]['plot_10f_datas'] = []

            crispresso2_info['refs'][ref_name]['plot_10g_roots'] = []
            crispresso2_info['refs'][ref_name]['plot_10g_captions'] = []
            crispresso2_info['refs'][ref_name]['plot_10g_datas'] = []

            for ind,sgRNA in enumerate(sgRNA_sequences):
                cut_point = sgRNA_cut_points[ind]
                plot_idxs = sgRNA_plot_idxs[ind]
                df_allele_around_cut=CRISPRessoShared.get_dataframe_around_cut(df_alleles.loc[df_alleles['Reference_Name'] == ref_name],cut_point,plot_half_window)

                #write alleles table to file
                allele_filename = _jp(ref_plot_name+'Alleles_frequency_table_around_sgRNA_'+sgRNA+'.txt')
                df_allele_around_cut.to_csv(allele_filename,sep='\t',header=True)
                crispresso2_info['refs'][ref_name]['allele_frequency_files'].append(os.path.basename(allele_filename))

                ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_half_window+1:cut_point+plot_half_window+1]
                fig_filename_root = _jp('9.'+ref_plot_name+'Alleles_frequency_table_around_sgRNA_'+sgRNA)
                n_good = df_allele_around_cut.ix[df_allele_around_cut['%Reads']>=args.min_frequency_alleles_around_cut_to_plot].shape[0]
                if not args.suppress_plots and n_good > 0:
                    CRISPRessoPlot.plot_alleles_table(ref_seq_around_cut,df_alleles=df_allele_around_cut,fig_filename_root=fig_filename_root,
                        MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot,MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot,SAVE_ALSO_PNG=save_png,base_editor_output=args.base_editor_output,sgRNA_intervals=sgRNA_intervals)
                    crispresso2_info['refs'][ref_name]['plot_9_roots'].append(os.path.basename(fig_filename_root))
                    crispresso2_info['refs'][ref_name]['plot_9_captions'].append("Figure 9: Visualization of the distribution of identified alleles around each the cleavage site for the guide " + sgRNA + ". Nucleotides are indicated by unique colors (A = green; C = red; G = yellow; T = purple). Substitutions are shown in bold font. Red rectangles highlight inserted sequences. Horizontal dashed lines indicate deleted sequences. The vertical dashed line indicates the predicted cleavage site.")
                    crispresso2_info['refs'][ref_name]['plot_9_datas'].append([('Allele frequency table',os.path.basename(allele_filename))])



                if not args.crispresso1_mode and args.base_editor_output:
                    ###guide-specific base editor plots
                    plot_ref_seq = ref_seq_around_cut
                    plot_nuc_pcts = df_nuc_pct_all.iloc[:,plot_idxs]
                    plot_nuc_freqs = df_nuc_freq_all.iloc[:,plot_idxs]

                    #get computation window in plotted region
                    is_window = np.zeros(ref_len)
                    for ind in include_idxs_list:
                        is_window[ind] = 1
                    plot_is_window = np.zeros(len(plot_idxs)) #binary array whether base should be plotted
                    plot_quant_window_idxs = []
                    for ind, loc in enumerate(plot_idxs):
                        plot_is_window[ind] = is_window[loc]
                        if is_window[loc]:
                            plot_quant_window_idxs.append(ind-2)

                    from_nuc_indices = [pos for pos, char in enumerate(list(plot_nuc_pcts.columns.values)) if char == args.conversion_nuc_from]
                    just_sel_nuc_pcts = plot_nuc_pcts.iloc[:,from_nuc_indices].copy() #only nucleotides targeted by base editing
                    #just_sel_nuc_pcts.columns = [char + str(pos+1) for pos,char in enumerate(list(just_sel_nuc_pcts.columns.values))]
                    just_sel_nuc_pcts.columns = [args.conversion_nuc_from + str(pos+1) for pos in from_nuc_indices]
                    just_sel_nuc_freqs = plot_nuc_freqs.iloc[:,from_nuc_indices].copy()
                    just_sel_nuc_freqs.columns = [args.conversion_nuc_from + str(pos+1) for pos in from_nuc_indices]

                    quant_window_sel_nuc_pct_filename = _jp(ref_plot_name + 'Selected_nucleotide_percentage_table_around_sgRNA_'+sgRNA+'.txt')
                    just_sel_nuc_pcts.to_csv(quant_window_sel_nuc_pct_filename,sep='\t',header=True,index=True)
#                   not storing the name because it is unique to this sgRNA
#                    crispresso2_info['quant_window_sel_nuc_pct_filename'] = os.path.basename(quant_window_sel_nuc_pct_filename)

                    quant_window_sel_nuc_freq_filename = _jp(ref_plot_name + 'Selected_nucleotide_frequency_table_around_sgRNA_'+sgRNA+'.txt')
                    just_sel_nuc_freqs.to_csv(quant_window_sel_nuc_freq_filename,sep='\t',header=True,index=True)
#                   not storing the name because it is unique to this sgRNA
#                    crispresso2_info['quant_window_sel_nuc_freq_filename'] = os.path.basename(quant_window_sel_nuc_freq_filename)

                    #print table showing all nuc frequencies (sum to total alleles) (in entire region)

        #                CRISPRessoPlot.plot_nuc_freqs(
        #                    df_nuc_freq = df_nuc_freq,
        #                    tot_aln_reads = tot_aln_reads,
        #                    plot_title = get_plot_title_with_ref_name('Nucleotide Frequencies',ref_name),
        #                    fig_filename_root = _jp('14a.'+ref_name+'.nucleotide_frequency'),
        #                    save_also_png = save_png
        #                    )

                    if not args.suppress_plots:

                        fig_filename_root = _jp('10d.'+ref_plot_name+'Log2_nucleotide_frequency_around_sgRNA_'+sgRNA)
                        CRISPRessoPlot.plot_log_nuc_freqs(
                            df_nuc_freq = plot_nuc_freqs,
                            tot_aln_reads = tot_aln_reads,
                            plot_title = get_plot_title_with_ref_name('Log2 Nucleotide Frequencies Around sgRNA ' + sgRNA,ref_name),
                            fig_filename_root = fig_filename_root,
                            save_also_png = save_png,
                            quantification_window_idxs = plot_quant_window_idxs
                            )
                        crispresso2_info['refs'][ref_name]['plot_10d_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['refs'][ref_name]['plot_10d_captions'].append("Figure 10d: Log2 nucleotide frequencies for each position in the plotting window around the sgRNA " + sgRNA + ". The quantification window is outlined by the dotted box.")
                        crispresso2_info['refs'][ref_name]['plot_10d_datas'].append([])


                        fig_filename_root = _jp('10e.'+ref_plot_name+'Selected_conversion_at_'+args.conversion_nuc_from+'s_around_sgRNA_'+sgRNA)
                        CRISPRessoPlot.plot_conversion_at_sel_nucs(
                            df_subs = plot_nuc_pcts,
                            ref_name = ref_name,
                            ref_sequence = plot_ref_seq,
                            plot_title = get_plot_title_with_ref_name('Substitution Frequencies at '+args.conversion_nuc_from+'s around sgRNA ' + sgRNA,ref_name),
                            conversion_nuc_from = args.conversion_nuc_from,
                            fig_filename_root = fig_filename_root,
                            save_also_png = save_png
                            )
                        crispresso2_info['refs'][ref_name]['plot_10e_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['refs'][ref_name]['plot_10e_captions'].append("Figure 10e: Proportion of each base at each nucleotide targeted by base editors in the plotting window around the sgRNA " + sgRNA + ". The number of each target base is annotated on the reference sequence at the bottom of the plot.")
                        crispresso2_info['refs'][ref_name]['plot_10e_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from + 's',os.path.basename(quant_window_sel_nuc_freq_filename))])

                        fig_filename_root = _jp('10f.'+ref_plot_name+'Selected_conversion_no_ref_at_'+args.conversion_nuc_from+'s_around_sgRNA_'+sgRNA)
                        CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref(
                            df_subs = plot_nuc_pcts,
                            ref_name = ref_name,
                            ref_sequence = plot_ref_seq,
                            plot_title = get_plot_title_with_ref_name('Substitution Frequencies at '+args.conversion_nuc_from+'s around sgRNA ' + sgRNA,ref_name),
                            conversion_nuc_from = args.conversion_nuc_from,
                            fig_filename_root = fig_filename_root,
                            save_also_png = save_png
                            )
                        crispresso2_info['refs'][ref_name]['plot_10f_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['refs'][ref_name]['plot_10f_captions'].append("Figure 10f: Non-reference base proportions. For target nucleotides in the plotting window, this plot shows the proportion of non-reference (non-"+args.conversion_nuc_from + ") bases as a percentage of all non-reference sequences. The number of each target base is annotated on the reference sequence at the bottom of the plot.")
                        crispresso2_info['refs'][ref_name]['plot_10f_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from + 's',os.path.basename(quant_window_sel_nuc_freq_filename))])

                        fig_filename_root = _jp('10g.'+ref_plot_name+'Selected_conversion_no_ref_scaled_at_'+args.conversion_nuc_from+'s_around_sgRNA_'+sgRNA)
                        CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref_scaled(
                            df_subs = plot_nuc_pcts,
                            ref_name = ref_name,
                            ref_sequence = plot_ref_seq,
                            plot_title = get_plot_title_with_ref_name('Substitution Frequencies at '+args.conversion_nuc_from+'s around sgRNA ' + sgRNA,ref_name),
                            conversion_nuc_from = args.conversion_nuc_from,
                            fig_filename_root = fig_filename_root,
                            save_also_png = save_png
                            )
                        crispresso2_info['refs'][ref_name]['plot_10g_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['refs'][ref_name]['plot_10g_captions'].append("Figure 10g: Non-reference base counts. For target nucleotides in the plotting window, this plot shows the number of non-reference (non-" + args.conversion_nuc_from + ") bases. The number of each target base is annotated on the reference sequence at the bottom of the plot.")
                        crispresso2_info['refs'][ref_name]['plot_10g_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from +'s',os.path.basename(quant_window_sel_nuc_freq_filename))])

            #END GUIDE SPECIFIC PLOTS

        #(5, 6) GLOBAL frameshift analyses plots
        if args.coding_seq:
            global_MODIFIED_FRAMESHIFT = 0
            global_MODIFIED_NON_FRAMESHIFT = 0
            global_NON_MODIFIED_NON_FRAMESHIFT = 0
            global_SPLICING_SITES_MODIFIED = 0
            global_hists_frameshift = defaultdict(lambda :0)
            global_hists_inframe = defaultdict(lambda :0)

            global_count_total = 0
            global_count_modified = 0
            global_count_unmodified = 0
            global_exon_len_mods = []

            for ref_name in ref_names:
                if refs[ref_name]['contains_coding_seq']: #PERFORM FRAMESHIFT ANALYSIS
                    if ref_name == "HDR":
                        global_MODIFIED_FRAMESHIFT += counts_modified_frameshift[ref_name]
                        global_MODIFIED_NON_FRAMESHIFT += counts_modified_non_frameshift[ref_name]
                        global_NON_MODIFIED_NON_FRAMESHIFT += counts_non_modified_non_frameshift[ref_name]
                        global_SPLICING_SITES_MODIFIED += counts_splicing_sites_modified[ref_name]

                        #for HDR, add all unmodified reads to those that have modifications not in exons
                        global_NON_MODIFIED_NON_FRAMESHIFT += counts_unmodified[ref_name]
                    else:
                        global_MODIFIED_FRAMESHIFT += counts_modified_frameshift[ref_name]
                        global_MODIFIED_NON_FRAMESHIFT += counts_modified_non_frameshift[ref_name]
                        global_NON_MODIFIED_NON_FRAMESHIFT += counts_non_modified_non_frameshift[ref_name]
                        global_SPLICING_SITES_MODIFIED += counts_splicing_sites_modified[ref_name]

                    for (exon_len,count) in hists_frameshift[ref_name].iteritems():
                        global_hists_frameshift[exon_len] += count
                    for (exon_len,count) in hists_inframe[ref_name].iteritems():
                        global_hists_inframe[exon_len] += count


                    global_count_total += counts_total[ref_name]
                    global_count_modified += counts_modified[ref_name]
                    global_count_unmodified += counts_unmodified[ref_name]

            if not args.suppress_plots:
                fig=plt.figure(figsize=(12*1.5,14.5*1.5))
                ax1 = plt.subplot2grid((6,3), (0, 0), colspan=3, rowspan=5)

                patches, texts, autotexts =ax1.pie([global_MODIFIED_FRAMESHIFT,\
                                                    global_MODIFIED_NON_FRAMESHIFT,\
                                                    global_NON_MODIFIED_NON_FRAMESHIFT],\
                                                   labels=['Frameshift mutation\n(%d reads)' %global_MODIFIED_FRAMESHIFT,\
                                                          'In-frame mutation\n(%d reads)' % global_MODIFIED_NON_FRAMESHIFT,\
                                                          'Noncoding mutation\n(%d reads)' %global_NON_MODIFIED_NON_FRAMESHIFT],\
                                                   explode=(0.0,0.0,0.0),\
                                                   colors=[(0.89019608,  0.29019608,  0.2, 0.8),(0.99215686,  0.73333333,  0.51764706,0.8),(0.99607843,  0.90980392,  0.78431373,0.8)],\
                                                   autopct='%1.1f%%')

                proptease = fm.FontProperties()
                proptease.set_size('xx-large')
                plt.setp(autotexts, fontproperties=proptease)
                plt.setp(texts, fontproperties=proptease)
                plot_root = _jp('5a.Global_frameshift_in-frame_mutations_pie_chart')
                plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_inches='tight')
                plt.close()
                crispresso2_info['plot_5a_root'] = os.path.basename(plot_root)
                crispresso2_info['plot_5a_caption'] = "Figure 5a: Frameshift analysis of coding sequence reads affected by modifications for all reads. Unmodified reference reads are excluded from this plot, and all HDR reads are included in this plot."
                crispresso2_info['plot_5a_data'] = []

                 #profiles-----------------------------------------------------------------------------------
                fig=plt.figure(figsize=(22,10))
                ax1=fig.add_subplot(2,1,1)
                x,y=map(np.array,zip(*[a for a in global_hists_frameshift.iteritems()]))
                if sum(global_hists_frameshift.values()) != 0:
                    y=y/float(sum(global_hists_frameshift.values()))*100
                ax1.bar(x-0.1,y)
                ax1.set_xlim(-30.5,30.5)
                ax1.set_frame_on(False)
                ax1.set_xticks([idx for idx in range(-30,31) if idx % 3])
                ax1.tick_params(which='both',      # both major and minor ticks are affected
                   bottom=False,      # ticks along the bottom edge are off
                   top=False,         # ticks along the top edge are off
                   labelbottom=True) # labels along the bottom edge are off)
                ax1.yaxis.tick_left()
                xmin, xmax = ax1.get_xaxis().get_view_interval()
                ymin, ymax = ax1.get_yaxis().get_view_interval()
                ax1.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if idx % 3]],rotation='vertical')
                plt.title('Global Frameshift profile')
                ax1.tick_params(axis='both', which='major', labelsize=24)
                ax1.tick_params(axis='both', which='minor', labelsize=24)
                plt.tight_layout()
                ax1.set_ylabel('Sequences % (no.)')
                y_label_values= np.round(np.linspace(0, min(100,max(ax1.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                ax1.set_yticks(y_label_values)
                ax1.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*sum(global_hists_frameshift.values())) for pct in y_label_values])

                ax2=fig.add_subplot(2,1,2)
                x,y=map(np.array,zip(*[a for a in global_hists_inframe.iteritems()]))
                if sum(global_hists_inframe.values()) > 0:
                    y=y/float(sum(global_hists_inframe.values()))*100
                #ax2.bar(x-0.5,y,color=(0,1,1,0.2))
                ax2.bar(x-0.1,y,color=(0,1,1,0.2))
                ax2.set_xlim(-30.5,30.5)
                ax2.set_frame_on(False)
                ax2.set_xticks([idx for idx in range(-30,31) if (idx % 3 ==0) ])
                ax2.tick_params(which='both',      # both major and minor ticks are affected
                   bottom=False,      # ticks along the bottom edge are off
                   top=False,         # ticks along the top edge are off
                   labelbottom=True) # labels along the bottom edge are off)
                ax2.yaxis.tick_left()
                xmin, xmax = ax2.xaxis.get_view_interval()
                ymin, ymax = ax2.yaxis.get_view_interval()
                ax2.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if (idx % 3==0)]],rotation='vertical',horizontalalignment="center")
                plt.title('Global In-frame profile')
                ax2.set_ylabel('Sequences % (no.)')
                y_label_values= np.round(np.linspace(0, min(100,max(ax2.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
                ax2.set_yticks(y_label_values)
                ax2.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*sum(global_hists_inframe.values())) for pct in y_label_values])

                ax2.tick_params(axis='both', which='major', labelsize=24)
                ax2.tick_params(axis='both', which='minor', labelsize=24)
                plt.tight_layout()

                plot_root = _jp('6a.Global_frameshift_in-frame_mutation_profiles')
                plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_inches='tight')
                plt.close()
                crispresso2_info['plot_6a_root'] = os.path.basename(plot_root)
                crispresso2_info['plot_6a_caption'] = "Figure 6a: Frameshift and in-frame mutagenesis profiles for all reads indicating position affected by modification."
                crispresso2_info['plot_6a_data'] = []


                 #-----------------------------------------------------------------------------------------------------------
                fig=plt.figure(figsize=(12*1.5,12*1.5))
                ax=fig.add_subplot(1,1,1)
                patches, texts, autotexts =ax.pie([global_SPLICING_SITES_MODIFIED,\
                                                  (global_count_total - global_SPLICING_SITES_MODIFIED)],\
                                                  labels=['Potential splice sites modified\n(%d reads)' %global_SPLICING_SITES_MODIFIED,\
                                                          'Unmodified\n(%d reads)' % (global_count_total - global_SPLICING_SITES_MODIFIED)],\
                                                  explode=(0.0,0),\
                                                  colors=[(0.89019608,  0.29019608,  0.2, 0.8),(0.99607843,  0.90980392,  0.78431373,0.8)],\
                                                  autopct='%1.1f%%')
                proptease = fm.FontProperties()
                proptease.set_size('xx-large')
                plt.setp(autotexts, fontproperties=proptease)
                plt.setp(texts, fontproperties=proptease)
                plot_root = _jp('8a.Global_potential_splice_sites_pie_chart')
                plt.savefig(plot_root+'.pdf',pad_inches=1,bbox_inches='tight')
                if save_png:
                    plt.savefig(plot_root+'.png',bbox_inches='tight')
                plt.close()
                crispresso2_info['plot_8a_root'] = os.path.basename(plot_root)
                crispresso2_info['plot_8a_caption'] = "Figure 8a: Predicted impact on splice sites for all reads. Potential splice sites modified refers to reads in which the either of the two intronic positions adjacent to exon junctions are disrupted."
                crispresso2_info['plot_8a_data'] = []

            #end global coding seq plots


        info('Done!')

        if not args.keep_intermediate:
            info('Removing Intermediate files...')

            if args.fastq_r2!='':
                files_to_remove=[processed_output_filename,flash_hist_filename,flash_histogram_filename,\
                        flash_not_combined_1_filename,flash_not_combined_2_filename]
            else:
                files_to_remove=[processed_output_filename]

            if args.trim_sequences and args.fastq_r2!='':
                files_to_remove+=[output_forward_paired_filename,output_reverse_paired_filename,\
                        output_forward_unpaired_filename,output_reverse_unpaired_filename]

            if args.split_paired_end:
                files_to_remove+=splitted_files_to_remove

            if args.min_average_read_quality>0 or args.min_single_bp_quality>0 or args.min_bp_quality_or_N>0:
                if args.fastq_r2!='':
                    files_to_remove+=[args.fastq_r1,args.fastq_r2]
                else:
                    files_to_remove+=[args.fastq_r1]

            for file_to_remove in files_to_remove:
                try:
                    if os.path.islink(file_to_remove):
                        os.unlink(file_to_remove)
                    else:
                        os.remove(file_to_remove)
                except Exception as e:
                    warn('Skipping removal of: %s: %s' %(file_to_remove,e))

        crispresso2_info['counts_total'] = counts_total
        crispresso2_info['counts_modified'] = counts_modified
        crispresso2_info['counts_unmodified'] = counts_unmodified
        crispresso2_info['counts_discarded'] = counts_discarded

        crispresso2_info['counts_insertion'] = counts_insertion
        crispresso2_info['counts_deletion'] = counts_deletion
        crispresso2_info['counts_substitution'] = counts_substitution

        crispresso2_info['counts_only_insertion'] = counts_only_insertion
        crispresso2_info['counts_only_deletion'] = counts_only_deletion
        crispresso2_info['counts_only_substitution'] = counts_only_substitution
        crispresso2_info['counts_insertion_and_deletion'] = counts_insertion_and_deletion
        crispresso2_info['counts_insertion_and_substitution'] = counts_insertion_and_substitution
        crispresso2_info['counts_deletion_and_substitution'] = counts_deletion_and_substitution
        crispresso2_info['counts_insertion_and_deletion_and_substitution'] = counts_insertion_and_deletion_and_substitution
        crispresso2_info['counts_modified_frameshift'] = counts_modified_frameshift
        crispresso2_info['counts_modified_non_frameshift'] = counts_modified_non_frameshift
        crispresso2_info['counts_non_modified_non_frameshift'] = counts_non_modified_non_frameshift
        crispresso2_info['counts_splicing_sites_modified'] = counts_splicing_sites_modified
        crispresso2_info['class_counts'] = class_counts

        end_time =  datetime.now()
        end_time_string =  end_time.strftime('%Y-%m-%d %H:%M:%S')
        running_time = end_time - start_time
        running_time_string =  str(running_time)

        crispresso2_info['end_time'] = end_time
        crispresso2_info['end_time_string'] = end_time_string
        crispresso2_info['running_time'] = running_time
        crispresso2_info['running_time_string'] = running_time_string


        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = os.path.join(OUTPUT_DIRECTORY,"CRISPResso2_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_report(crispresso2_info,report_name,OUTPUT_DIRECTORY,_ROOT)
            crispresso2_info['report_location'] = report_name
            crispresso2_info['report_filename'] = os.path.basename(report_name)

        cp.dump(crispresso2_info, open(crispresso2_info_file, 'wb' ) )

        info('Analysis Complete!')
        print(CRISPRessoShared.get_crispresso_footer())

        sys.exit(0)

    except CRISPRessoShared.NTException as e:
        print_stacktrace_if_debug()
        error('Alphabet error, please check your input.\n\nERROR: %s' % e)
        sys.exit(1)
    except CRISPRessoShared.SgRNASequenceException as e:
        print_stacktrace_if_debug()
        error('sgRNA error, please check your input.\n\nERROR: %s' % e)
        sys.exit(2)
    except CRISPRessoShared.TrimmomaticException as e:
        print_stacktrace_if_debug()
        error('Trimming error, please check your input.\n\nERROR: %s' % e)
        sys.exit(4)
    except CRISPRessoShared.FlashException as e:
        print_stacktrace_if_debug()
        error('Merging error, please check your input.\n\nERROR: %s' % e)
        sys.exit(5)
    except CRISPRessoShared.BadParameterException as e:
        print_stacktrace_if_debug()
        error('Parameter error, please check your input.\n\nERROR: %s' % e)
        sys.exit(6)
    except CRISPRessoShared.NoReadsAlignedException as e:
        print_stacktrace_if_debug()
        error('Alignment error, please check your input.\n\nERROR: %s' % e)
        sys.exit(7)
    except CRISPRessoShared.AutoException as e:
        print_stacktrace_if_debug()
        error('Autorun error. This sample cannot be run in auto mode.\n\nERROR: %s' % e)
        sys.exit(8)
    except CRISPRessoShared.AlignmentException as e:
        print_stacktrace_if_debug()
        error('Alignment error, please check your input.\n\nERROR: %s' % e)
        sys.exit(8)
    except CRISPRessoShared.ExonSequenceException as e:
        print_stacktrace_if_debug()
        error('Coding sequence error, please check your input.\n\nERROR: %s' % e)
        sys.exit(11)
    except CRISPRessoShared.DuplicateSequenceIdException as e:
        print_stacktrace_if_debug()
        error('Fastq file error, please check your input.\n\nERROR: %s' % e)
        sys.exit(12)
    except CRISPRessoShared.NoReadsAfterQualityFilteringException as e:
        print_stacktrace_if_debug()
        error('Filtering error, please check your input.\n\nERROR: %s' % e)
        sys.exit(13)
    except Exception as e:
        print_stacktrace_if_debug()
        error('Unexpected error, please check your input.\n\nERROR: %s' % e)
        sys.exit(-1)
