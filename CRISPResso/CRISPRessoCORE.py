#!/usr/bin/env python
# -*- coding: utf8 -*-

'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2017 The General Hospital Corporation. All Rights Reserved.
'''



import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import argparse
from collections import defaultdict
import errno
import gzip
import os
import re
import signal
import shutil
import subprocess as sb
import traceback
import unicodedata

if running_python3:
    import pickle as cp #python 3
else:
    import cPickle as cp #python 2.7

from CRISPResso import CRISPRessoCOREResources
from CRISPResso import CRISPRessoShared
from CRISPResso import CRISPRessoPlot
from CRISPResso import cnwalign

from datetime import datetime
present = datetime.now()
d1 = datetime.strptime('07/07/2018','%d/%m/%Y')
if present > d1:
    print('\nYour version of CRISPResso2 is out of date. Please download a new version.\n')
    sys.exit(1)

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

def get_data(path):
    return os.path.join(_ROOT, 'data', path)

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


def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        files_in_dir = os.listdir()
        raise BadParameterException("The specified file '"+filename + "' cannot be opened. Available files: " + str(files_in_dir))

def force_symlink(src, dst):

    if os.path.exists(dst) and os.path.samefile(src,dst):
        return

    try:
        os.symlink(src, dst)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)
        elif exc.errno == errno.EPROTO:
            #in docker on windows 7, symlinks don't work so well, so we'll just copy the file.
            shutil.copyfile(src, dst)

def get_avg_read_length_fastq(fastq_filename):
     cmd=('z' if fastq_filename.endswith('.gz') else '' ) +('cat < %s' % fastq_filename)+\
                  r''' | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' '''
     p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
     return int(p.communicate()[0].strip())

def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     return int(float(p.communicate()[0])/4.0)

import time
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

check_program('java')
check_program('flash')

#start = time.time()
sns=check_library('seaborn')
#end = time.time()
sns.set_context('poster')
sns.set(font_scale=2.2)
sns.set_style('white')

#########################################



###EXCEPTIONS############################
class FlashException(Exception):
    pass

class TrimmomaticException(Exception):
    pass

class NoReadsAlignedException(Exception):
    pass

class SgRNASequenceException(Exception):
    pass

class NTException(Exception):
    pass

class ExonSequenceException(Exception):
    pass

class DuplicateSequenceIdException(Exception):
    pass

class NoReadsAfterQualityFiltering(Exception):
    pass

class BadParameterException(Exception):
    pass

class AutoException(Exception):
    pass

#########################################


def process_fastq(fastq_filename,variantCache,ref_names,refs,args):
    """process_fastq processes each of the reads contained in a fastq file, given a cache of pre-computed variants
        fastqIn: name of fastq (e.g. output of FLASh)
            This file can be gzipped or plain text

        variantCache: dictionary of sequence>payload
            The payload is an array. The first element is the count of sequences producing that payload,
                and the successive items in the array are payload objects
            # payload object:
                ### from CRISPRessoCOREResources.find_indels_substitutions
                # 'all_insertion_positions' #arr with 1's where there are insertions (including those outside of include_idxs mask)
                # 'all_insertion_left_positions' #arr with 1's to the left of where the insertion occurs
                # 'insertion_positions' # arr with 1's where there are insertions (1bp before and 1bp after insertion) that overlap with include_idxs mask
                # 'insertion_coordinates' # one entry per insertion, tuple of (start,end)
                # 'insertion_sizes'
                # 'insertion_n'
                # 'all_deletion_positions' #arr with 1's where there are insertions
                # 'deletion_positions' #arr with 1's where there are insertions that overlap the include_idxs mask
                # 'deletion_coordinates' # one entry per deletion
                # 'deletion_sizes' # correspond to entries in 'deletion_coordinates'
                # 'deletion_n'
                # 'all_substitution_positions'
                # 'substitution_positions'
                # 'substitution_n'
                # 'substitution_values'
                # 'ref_positions'
                ### added in this function
                # 'closest_aln_name' # name of sequence that it most closely aligns to
                # 'closest_aln_length' # sequence length
                # 'classification' # MODIFIED or UNMODIFIED
                # 'aln_scores' # scores of alignment to each other reference sequence
                # 'aln_seq' # NW-aligned sequence
                # 'aln_ref' # NW-aligned sequence of corresponding reference (closest_aln_name)

        refNameList: list of reference names
        refs: dictionary of sequences name>ref object
            ##ref object:
                # 'name'
                # 'sequence'
                # 'sequence_length'
                # 'min_aln_score' #sequence must align with at least this score
                # 'cut_points'
                # 'gap_incentive' #incentive for gaps at each position of the reference - to force gaps at the cut points, the indices of these cut points are set to 1  i.e. gap_incentive[4] = 1 would incentivise alignments with insertions to the right of the 4th character in the reference, or deletions of the 4th character in the reference.
                # 'offset_plots'
                # 'sgRNA_intervals'
                # 'sgRNA_sequences'
                # 'contains_guide'
                # 'contains_coding_seq'
                # 'exon_positions'
                # 'exon_intervals'
                # 'splicing_positions'
                # 'include_idxs'
                # 'exclude_idxs'
                # 'plot_idxs'
                # 'idx_cloned_from' #if this reference didn't contain a guide, it was aligned to 'idx_cloned_from' reference, and cut_points, gap_incentive, sgRNA_intervals, and inculde_idx were cloned from it (at the appropriate indices)
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

    alnMatrixLoc = os.path.join(_ROOT,"EDNAFULL")
    check_file(alnMatrixLoc)
    alnMatrix = cnwalign.read_matrix(alnMatrixLoc)

    if (args.needleman_wunsch_gap_open > 0):
        raise BadParameterException("Needleman Wunsch gap open penalty must be <= 0")
    if (args.needleman_wunsch_gap_extend > 0):
        raise BadParameterException("Needleman Wunsch gap extend penalty must be <= 0")


    not_aln = {} #cache for reads that don't align

    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename)
    else:
        fastq_handle=open(fastq_filename)

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
            variantCache[fastq_seq][0] += 1

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
                fws1,fws2,fwscore=cnwalign.global_align(fastq_seq, refs[ref_name]['sequence'],matrix=alnMatrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                rvs1,rvs2,rvscore=cnwalign.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'],matrix=alnMatrix,gap_incentive=refs[ref_name]['gap_incentive'],gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
#                print "for " + ref_name + " got fws1: " + str(fws1) + " and fws2: " + str(fws2) + " score: " +str(fwscore)

                s1 = fws1
                s2 = fws2
                score = fwscore
                if (rvscore > fwscore):
                    s1 = rvs1
                    s2 = rvs2
                    score = rvscore
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
                variantCache[fastq_seq] = []
                variantCache[fastq_seq].append(1) #first element in variant cache is the number of sequences

                if not args.expand_ambiguous_alignments and len(best_match_names) > 1:
                    idx = 0
                    payload=CRISPRessoCOREResources.find_indels_substitutions(best_match_s1s[idx],best_match_s2s[idx],refs[best_match_names[idx]]['include_idxs'])
                    payload['closest_aln_name'] = 'AMBIGUOUS'
                    payload['aln_scores'] = aln_scores
                    payload['classification'] = ''


                    payload['aln_seq'] = best_match_s1s[idx]
                    payload['aln_ref'] = best_match_s2s[idx]
                    payload['closest_aln_length'] = refs[best_match_name]['sequence_length']

                    variantCache[fastq_seq].append(payload) #successive elements are payloads corresponding to mappings to references

                else:
                    for idx, best_match_name in enumerate(best_match_names):
                        payload=CRISPRessoCOREResources.find_indels_substitutions(best_match_s1s[idx],best_match_s2s[idx],refs[best_match_name]['include_idxs'])
                        payload['closest_aln_name'] = best_match_name
                        payload['aln_scores'] = aln_scores


                        # If there is an insertion/deletion/substitution in the target window, the read is modified.
                        if not args.ignore_deletions and payload['deletion_n'] > 0:
                            payload['classification'] = 'MODIFIED'
                        elif not args.ignore_insertions and payload['insertion_n'] > 0:
                            payload['classification'] = 'MODIFIED'
                        elif not args.ignore_substitutions and payload['substitution_n'] > 0:
                            payload['classification'] = 'MODIFIED'
                        else:
                            payload['classification'] = 'UNMODIFIED'


                        payload['aln_seq'] = best_match_s1s[idx]
                        payload['aln_ref'] = best_match_s2s[idx]
                        payload['closest_aln_length'] = refs[best_match_name]['sequence_length']

                        variantCache[fastq_seq].append(payload) #successive elements are payloads corresponding to mappings to references

            else:
                N_COMPUTED_NOTALN+=1
                not_aln[fastq_seq] = 1


    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS,N_COMPUTED_ALN,N_CACHED_ALN,N_COMPUTED_NOTALN,N_CACHED_NOTALN))
    alnStats = {"N_TOT_READS" : N_TOT_READS,
               "N_CACHED_ALN" : N_CACHED_ALN,
               "N_CACHED_NOTALN" : N_CACHED_NOTALN,
               "N_COMPUTED_ALN" : N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN" : N_COMPUTED_NOTALN}
    return(alnStats)


def add_hist(hist_to_add,hist_global):
    for key,value in hist_to_add.iteritems():
        hist_global[key]+=value
    return hist_global


def slugify(value): #adapted from the Django project

    value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '_', value).strip())
    value = unicode(re.sub('[-\s]+', '-', value))

    return str(value)

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
        raise Exception('Error in splitting read pairs from a single file')

    return output_filename_r1,output_filename_r2


def main():

    def print_stacktrace_if_debug():
        debug_flag = False
        if 'args' in vars() and 'debug' in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)

    try:
        description = ['~~~CRISPResso 2~~~','-Analysis of CRISPR/Cas9 outcomes from deep sequencing data-']
        header = CRISPRessoShared.get_crispresso_header(description=description,header_str=None)
        print(header)

        args = CRISPRessoShared.getCRISPRessoArgParser(_ROOT,requiredParams={'fastq_r1':True}).parse_args()

        alnMatrixLoc = os.path.join(_ROOT,"EDNAFULL")
        check_file(alnMatrixLoc)
        alnMatrix = cnwalign.read_matrix(alnMatrixLoc)

        #check files
        check_file(args.fastq_r1)
        if args.fastq_r2:
            check_file(args.fastq_r2)

        #normalize name and remove not allowed characters
        if args.name:
            clean_name=slugify(args.name)
            if args.name!= clean_name:
                warn('The specified name %s contained invalid characters and was changed to: %s' % (args.name,clean_name))
                args.name=clean_name

        #### ASSERT GUIDE(S)
        guides = []
        if args.guide_seq:
            for current_guide_seq in args.guide_seq.split(','):
                wrong_nt=CRISPRessoShared.find_wrong_nt(current_guide_seq)
                if wrong_nt:
                    raise NTException('The sgRNA sequence contains bad characters:%s'  % ' '.join(wrong_nt))
                guides.append(current_guide_seq)

        ###FRAMESHIFT SUPPORT###
        coding_seqs = []
        if args.coding_seq:
            for exon_seq in args.coding_seq.strip().upper().split(','):
                #check for wrong NT
                wrong_nt=CRISPRessoShared.find_wrong_nt(exon_seq)
                if wrong_nt:
                    raise NTException('The coding sequence contains bad characters:%s' % ' '.join(wrong_nt))

                coding_seqs.append(exon_seq)

        ####SET REFERENCES TO COMPARE###
        ref_names = [] #ordered list of names
        refs = {} #dict of ref_name > ref object

#        #if we should automatically infer amplicon sequence, pull out the most frequent read and assign it to be the amplicon
        if args.auto:
            #paste <(zcat R1.fastq.gz) <(zcat R2.fastq.gz) | head -n 500 | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | flash - --interleaved --to-stdout 2>/dev/null | awk '((NR-2)%4==0){print $1}' | sort | uniq -c | sort -nr | head
            number_of_reads_to_consider = 1000 * 4 #1000 fastq sequences (4 lines each)
            view_cmd_1 = 'cat'
            if args.fastq_r1.endswith('.gz'):
                view_cmd_1 = 'zcat'
            file_generation_command = "%s %s | head -n %d "%(view_cmd_1,args.fastq_r1,number_of_reads_to_consider)

            if args.fastq_r2:
                view_cmd_2 = 'cat'
                if args.fastq_r2.endswith('.gz'):
                    view_cmd_2 = 'zcat'
                file_generation_command = "paste <(%s %s) <(%s %s) | head -n %d | paste - - - - | awk -v OFS=\"\\n\" -v FS=\"\\t\" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | flash - --interleaved-input --min-overlap %d --to-stdout 2>/dev/null " %(view_cmd_1,args.fastq_r1,view_cmd_2,args.fastq_r2,number_of_reads_to_consider,args.min_paired_end_reads_overlap)
            count_frequent_cmd = file_generation_command + " | awk \"((NR-2)%4==0){print $1}\" | sort | uniq -c | sort -nr "
            def default_sigpipe():
                signal.signal(signal.SIGPIPE, signal.SIG_DFL)

            p = sb.Popen(count_frequent_cmd, shell=True,stdout=sb.PIPE,preexec_fn=default_sigpipe)
            top_unaligned = p.communicate()[0]
            if p.poll() != 0:
                raise AutoException('Cannot retrieve most frequent amplicon sequences. Got nonzero return code.')
            seq_lines = top_unaligned.strip().split("\n")
            if len(seq_lines) == 0:
                raise AutoException('Cannot parse any frequent amplicons sequences.')

            curr_amplicon_id = 1

            amplicon_seq_arr = []
            amplicon_name_arr = []

            #add most frequent amplicon to the list
            count,seq = seq_lines[0].strip().split()
            amplicon_seq_arr.append(seq)
            amplicon_name_arr.append('Amplicon')
            curr_amplicon_id += 1

            #for the remainder of the amplicons, test them before adding
            for i in range(1,len(seq_lines)):
                count,seq = seq_lines[i].strip().split()
                last_count,last_seq = seq_lines[i-1].strip().split()
                #if this allele is present in at least 20% of the samples
                if float(last_count)/float(number_of_reads_to_consider) > 0.01:
                    for amp_seq in amplicon_seq_arr:
                        ref_incentive = np.zeros(len(amp_seq)+1,dtype=np.int)
                        fws1,fws2,fwscore=cnwalign.global_align(seq,amp_seq,matrix=alnMatrix,gap_incentive=ref_incentive,gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                        rvs1,rvs2,rvscore=cnwalign.global_align(CRISPRessoShared.reverse_complement(seq),amp_seq,matrix=alnMatrix,gap_incentive=ref_incentive,gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
                        #if the sequence is similar to a previously-seen read, don't add it
                        if fwscore > 0.95 or rvscore > 0.95:
                            continue
                        else:
                            amplicon_seq_arr.append(seq)
                            amplicon_name_arr.append('Amplicon_%d'%curr_amplicon_id)
                            curr_amplicon_id += 1
                            continue
                else:
                    break


            amplicon_min_alignment_score_arr = []
            plural_string = ""
            if len(amplicon_seq_arr) > 1:
                plural_string = "s"
            info("Auto-detected %d reference amplicon%s"%(len(amplicon_seq_arr),plural_string))
        else: #not auto
            amplicon_seq_arr = args.amplicon_seq.split(",")
            amplicon_name_arr = args.amplicon_name.split(",")
            amplicon_min_alignment_score_arr = args.amplicon_min_alignment_score.split(",")

        found_guide_seq = [False]*len(guides)
        found_coding_seq = [False]*len(coding_seqs)

        max_amplicon_len = 0 #for flash

        for idx,seq in enumerate(amplicon_seq_arr):
            this_seq = seq.strip().upper()
            this_seq_length = len(this_seq)
            if this_seq_length > max_amplicon_len:
                max_amplicon_len = this_seq_length

            this_name = 'Reference'+str(idx)
            if idx < len(amplicon_name_arr):
                this_name = amplicon_name_arr[idx]

            wrong_nt=CRISPRessoShared.find_wrong_nt(this_seq)
            if wrong_nt:
                raise NTException('Amplicon sequence %d (%s) contains invalid characters:%s' % idx,this_name, ' '.join(wrong_nt))

            this_min_aln_score = args.default_min_aln_score
            if idx < len(amplicon_min_alignment_score_arr):
                this_min_aln_score = float(amplicon_min_alignment_score_arr[idx])

            # Calculate cut sites for this reference
            this_contains_guide = False
            this_cut_points = []
            this_sgRNA_intervals = []
            this_sgRNA_sequences = []
            this_offset_plots = []

            for guide_idx, current_guide_seq in enumerate(guides):
                offset_fw=args.cleavage_offset+len(current_guide_seq)-1
                offset_rc=(-args.cleavage_offset)-1
                this_cut_points+=[m.start() + offset_fw for m in re.finditer(current_guide_seq, this_seq)]+\
                                 [m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse_complement(current_guide_seq), this_seq)]+\
                                 [m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse(current_guide_seq), this_seq)]
                this_sgRNA_intervals+=[(m.start(),m.start()+len(current_guide_seq)-1) for m in re.finditer(current_guide_seq, this_seq)]+\
                                      [(m.start(),m.start()+len(current_guide_seq)-1) for m in re.finditer(CRISPRessoShared.reverse_complement(current_guide_seq), this_seq)]+\
                                      [(m.start(),m.start()+len(current_guide_seq)-1) for m in re.finditer(CRISPRessoShared.reverse(current_guide_seq), this_seq)]
                this_sgRNA_sequences.append(current_guide_seq)

                if this_cut_points:
                    found_guide_seq[guide_idx] = True
                    this_contains_guide=True
                    this_offset_plots.append(1)
                else:
                    this_offset_plots.append(0)

            # Calculate coding sequence for this reference
            this_exon_positions = set()
            this_exon_intervals = []
            this_splicing_positions = []
            this_contains_coding_seq = False
            for exon_idx, exon_seq in enumerate(coding_seqs):
                st_exon = this_seq.find(exon_seq)
                if st_exon >= 0:
                    found_coding_seq[exon_idx] = True
                    this_contains_coding_seq = True
                    en_exon = st_exon + len(exon_seq)  # this do not include the upper bound as usual in python
                    this_exon_intervals.append((st_exon, en_exon))
                    this_exon_positions = this_exon_positions.union(set(range(st_exon, en_exon)))

                    # consider 2 base pairs before and after each exon
                    this_splicing_positions += [max(0, st_exon - 2), max(0, st_exon - 1), min(this_seq_length - 1, en_exon), min(this_seq_length - 1, en_exon + 1)]

            this_exon_positions = sorted(this_exon_positions)

            # protect from the wrong splitting of exons by the users to avoid false splicing sites
            this_splicing_positions = set(this_splicing_positions).difference(this_exon_positions)


            #create mask of positions in which to include/exclude indels for the analysis window
            this_include_idxs=[]
            #first, if base editor mode is set, set the guide as the analysis swindow
            if args.base_editor_mode and len(this_sgRNA_intervals) > 0:
                for sgRNA_int in this_sgRNA_intervals:
                    this_include_idxs.extend(range(sgRNA_int[0],sgRNA_int[1]))
            #otherwise, if exact coordinates have been given, set those
            elif args.analysis_window_coordinates is not None and len(args.analysis_window_coordinates.split(",")) > idx :
                theseCoords = args.analysis_window_coordinates.split(",")[idx].split("_")
                for coord in theseCoords:
                    coordRE = re.match(r'^(\d+)-(\d+)$',coord)
                    if coordRE:
                        start = int(coordRE.group(1))
                        end = int(coordRE.group(2)) + 1
                        if end > this_seq_length:
                            raise NTException("End coordinate " + str(end) + " for '" + str(coord) + "' in '" + str(theseCoords) + "' is longer than the sequence length ("+str(this_seq_length)+")")
                        this_include_idxs.extend(range(start,end))
                    else:
                        raise NTException("Cannot parse analysis window coordinate '" + str(coord) + "' in '" + str(theseCoords) + "'. Coordinates must be given in the form start-end e.g. 5-10 . Please check the --analysis_window_coordinate parameter.")
            elif this_cut_points and args.window_around_sgrna>0:
                if args.crispresso1_mode:
                    half_window=max(1,args.window_around_sgrna/2)
                    for cut_p in this_cut_points:
                        st=max(0,cut_p-half_window+1)
                        en=min(len(seq)-1,cut_p+half_window+1)
                        this_include_idxs.extend(range(st,en))
                else:
                    for cut_p in this_cut_points:
                        st=max(0,cut_p-args.window_around_sgrna+1)
                        en=min(len(seq)-1,cut_p+args.window_around_sgrna+1)
                        this_include_idxs.extend(range(st,en))
            else:
               this_include_idxs=range(len(seq))

            this_exclude_idxs=[]

            if args.exclude_bp_from_left:
               this_exclude_idxs+=range(args.exclude_bp_from_left)

            if args.exclude_bp_from_right:
               this_exclude_idxs+=range(this_seq_length)[-args.exclude_bp_from_right:]

            #flatten the arrays to avoid errors with old numpy library
            this_include_idxs=np.ravel(this_include_idxs)
            this_exclude_idxs=np.ravel(this_exclude_idxs)

            this_include_idxs=set(np.setdiff1d(this_include_idxs,this_exclude_idxs))

            this_plot_idxs=[]
            if this_cut_points and args.offset_around_cut_to_plot>0:
                window_around_cut = args.offset_around_cut_to_plot
                if args.crispresso1_mode:
                    window_around_cut=max(1,args.offset_around_cut_to_plot/2)
                for cut_p in this_cut_points:
                    if cut_p - window_around_cut + 1 < 0:
                        raise BadParameterException('Offset around cut would extend to the left of the amplicon. Please decrease offset_around_cut_to_plot parameter')
                    if cut_p - window_around_cut > len(seq)-1:
                        raise BadParameterException('Offset around cut would be greater than sequence length . Please decrease offset_around_cut_to_plot parameter')
                    st=max(0,cut_p-window_around_cut+1)
                    en=min(len(seq)-1,cut_p+window_around_cut+1)
                    this_plot_idxs.append(range(st,en))
            else:
               this_plot_idxs=range(len(seq))

            this_plot_idxs = np.ravel(this_plot_idxs)

            this_gap_incentive = np.zeros(len(seq)+1,dtype=np.int)
            for cut_point in this_cut_points:
                this_gap_incentive[cut_point+1] = args.needleman_wunsch_gap_incentive

            refObj = {'name':this_name,
                   'sequence':seq,
                   'sequence_length':this_seq_length,
                   'min_aln_score':this_min_aln_score,
                   'cut_points':this_cut_points,
                   'gap_incentive':this_gap_incentive,
                   'offset_plots':np.array(this_offset_plots),
                   'sgRNA_intervals':this_sgRNA_intervals,
                   'sgRNA_sequences':this_sgRNA_sequences,
                   'contains_guide':this_contains_guide,
                   'contains_coding_seq':this_contains_coding_seq,
                   'exon_positions':this_exon_positions,
                   'exon_intervals':this_exon_intervals,
                   'splicing_positions':this_splicing_positions,
                   'include_idxs':this_include_idxs,
                   'exclude_idxs':this_exclude_idxs,
                   'plot_idxs':this_plot_idxs,
                   'idx_cloned_from':None,
                   }
            ref_names.append(this_name)
            refs[this_name] = refObj


        #throw error if guides, or coding seqs not found in any reference
        if args.guide_seq:
            for idx, presence_bool in enumerate(found_guide_seq):
                if not presence_bool:
                    raise SgRNASequenceException('The guide sequence %d (%s) provided is not present in the amplicon sequences!\n\nPlease check your input!' % (idx, guides[idx]))

        if args.coding_seq:
            for idx, presence_bool in enumerate(found_coding_seq):
                if not presence_bool:
                    raise ExonSequenceException('The coding subsequence %d (%s) provided is not contained in any amplicon sequence!\n\nPlease check your input!' % (idx,coding_seqs[idx]))


        #clone cut points and include idx from first reference where those are set
        clone_ref_name = None
        for ref_name in ref_names:
            cut_points = refs[ref_name]['cut_points']
            if cut_points:
                if len(ref_names) > 1:
                    info("Using cut points from %s as template for other references"%ref_name)
                clone_ref_name = ref_name
                break

        if clone_ref_name is not None:
            for ref_name in ref_names:
                cut_points = refs[ref_name]['cut_points']
                sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
                if cut_points:
                    if len(ref_names) > 1:
                        info("Reference '%s' has cut points defined: %s. Not cloning."%(ref_name,cut_points))
                    continue
                if sgRNA_intervals:
                    if len(ref_names) > 1:
                        info("Reference '%s' has sgRNA_intervals defined: %s. Not cloning."%(ref_name,sgRNA_intervals))
                    continue

                fws1,fws2,fwscore=cnwalign.global_align(refs[ref_name]['sequence'], refs[clone_ref_name]['sequence'],matrix=alnMatrix,gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,gap_incentive=refs[clone_ref_name]['gap_incentive'])
                if fwscore < 60:
                    continue
                info("Reference '%s' has NO cut points or sgRNA intervals idxs defined. Cloning from %s."%(ref_name,clone_ref_name))
                print('fwscore is '+str(fwscore))
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

                this_cut_points = [s1inds[X] for X in refs[clone_ref_name]['cut_points']]
                this_gap_incentive = np.zeros(len(seq)+1,dtype=np.int)
                for cut_point in this_cut_points:
                    this_gap_incentive[cut_point + 1] = args.needleman_wunsch_gap_incentive

                this_sgRNA_intervals = []
                for (sgRNA_interval_start,sgRNA_interval_end) in refs[clone_ref_name]['sgRNA_intervals']:
                    this_sgRNA_intervals.append((s1inds[sgRNA_interval_start],s1inds[sgRNA_interval_end]))

                this_include_idxs = [s1inds[X] for X in refs[clone_ref_name]['include_idxs']]
                #subtract any indices in 'exclude_idxs' -- e.g. in case some of the cloned include_idxs were near the read ends (exlcuded)
                this_exclude_idxs = set(refs[ref_name]['exclude_idxs'])
                this_include_idxs = set(np.setdiff1d(this_include_idxs,this_exclude_idxs))

                refs[ref_name]['cut_points'] = this_cut_points
                refs[ref_name]['gap_incentive'] = this_gap_incentive
                refs[ref_name]['sgRNA_intervals'] = this_sgRNA_intervals
                refs[ref_name]['include_idxs'] = this_include_idxs
                refs[ref_name]['idx_cloned_from'] = clone_ref_name

        #create output directory
        get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')

        if not args.name:
            if args.fastq_r2!='':
                database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1),get_name_from_fasta(args.fastq_r2))
            else:
                database_id='%s' % get_name_from_fasta(args.fastq_r1)

        else:
            database_id=args.name


        OUTPUT_DIRECTORY='CRISPResso_on_%s' % database_id

        if args.output_folder:
            OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPResso_RUNNING_LOG.txt')


        try:
            os.makedirs(OUTPUT_DIRECTORY)
            info('Creating Folder %s' % OUTPUT_DIRECTORY)
#            info('Done!') #crispresso2 doesn't announce that the folder is created... save some electricity here
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        finally:
            logging.getLogger().addHandler(logging.FileHandler(log_filename))

            with open(log_filename,'w+') as outfile:
                outfile.write('CRISPResso version %s\n[Command used]:\n%s\n\n[Execution log]:\n' %(CRISPRessoShared.__version__,' '.join(sys.argv)))



        if args.split_paired_end:
            if args.fastq_r2!='':
                raise BadParameterException('The option --split_paired_end is available only when a single fastq file is specified!')
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
                force_symlink(os.path.abspath(args.fastq_r1),symlink_filename)
                output_forward_filename=symlink_filename
            else:
                output_forward_filename=_jp('reads.trimmed.fq.gz')
                #Trimming with trimmomatic
                cmd='java -jar %s SE -phred33 %s  %s %s >>%s 2>&1'\
                % (get_data('trimmomatic-0.33.jar'),args.fastq_r1,
                   output_forward_filename,
                   args.trimmomatic_options_string.replace('NexteraPE-PE.fa','TruSeq3-SE.fa'),
                   log_filename)
                #print cmd
                TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)

                if TRIMMOMATIC_STATUS:
                        raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')


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
                cmd='java -jar %s PE -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
                    % (get_data('trimmomatic-0.33.jar'),
                        args.fastq_r1,args.fastq_r2,output_forward_paired_filename,
                        output_forward_unpaired_filename,output_reverse_paired_filename,
                        output_reverse_unpaired_filename,args.trimmomatic_options_string,log_filename)
                #print cmd
                TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)
                if TRIMMOMATIC_STATUS:
                    raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

                info('Done!')


            info('Estimating average read length...')
            if get_n_reads_fastq(output_forward_paired_filename):
                avg_read_length=get_avg_read_length_fastq(output_forward_paired_filename)
            else:
               raise NoReadsAfterQualityFiltering('No reads survived the average or single bp quality filtering.')

            #Merging with Flash
            info('Merging paired sequences with Flash...')
            cmd='flash %s %s --min-overlap %d --max-overlap %d -f %d -z -d %s >>%s 2>&1' %\
            (output_forward_paired_filename,
                 output_reverse_paired_filename,
                 args.min_paired_end_reads_overlap,
                 max_amplicon_len,
                 avg_read_length,
                 OUTPUT_DIRECTORY,log_filename)
#            cmd='flash %s %s --min-overlap %d -f %d -r %d -s %d  -z -d %s >>%s 2>&1' %\
#            (output_forward_paired_filename,
#                 output_reverse_paired_filename,
#                 args.min_paired_end_reads_overlap,
#                 len_amplicon,avg_read_length,
#                 std_fragment_length,
#                 OUTPUT_DIRECTORY,log_filename)

            FLASH_STATUS=sb.call(cmd,shell=True)
            if FLASH_STATUS:
                raise FlashException('Flash failed to run, please check the log file.')

            info('Done!')

            flash_hist_filename=_jp('out.hist')
            flash_histogram_filename=_jp('out.histogram')
            flash_not_combined_1_filename=_jp('out.notCombined_1.fastq.gz')
            flash_not_combined_2_filename=_jp('out.notCombined_2.fastq.gz')

            processed_output_filename=_jp('out.extendedFrags.fastq.gz')


        if args.left_adapter_umi_trim_seq is not None:
            #first, flip all reads into alignment with (the first) amplicon
            temp_flip_filename = jp(processed_output_filename + '.temp_flipped_forward.fq')
            flip_command = "python flipReadsForward.py --fastq %s --fastq_out %s --amplicon %s"%(processed_output_filename,temp_flip_filename,amplicons[0])
            FLIP_STATUS=sb.call(flip_command,shell=True)
            if FLIP_STATUS:
                raise UMIException('UMI deduplication flipping failed to run, please check the log file.')

            #next, deduplicate and trim adapter sequences
            left_trim_seq_string = args.left_adapter_umi_trim_seq if args.left_adapter_umi_trim_seq is not None else ""
            right_trim_seq_string = args.right_adapter_umi_trim_seq if args.right_adapter_umi_trim_seq is not None else ""
            dedupUMI_command = "python dedupTrimUMI.py.py --fastq %s --fastq_out %s %s %s"%(temp_flip_filename,processed_output_filename,left_adapter_umi_trim_seq_string,right_adapter_umi_trim_seq_string)
            UMI_STATUS=sb.call(flip_command,shell=True)
            if FLIP_STATUS:
                raise UMIException('UMI deduplication failed to run, please check the log file.')

        #count reads
        N_READS_INPUT=get_n_reads_fastq(args.fastq_r1)
        N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)
        if N_READS_AFTER_PREPROCESSING == 0:
            raise NoReadsAfterQualityFiltering('No reads in input or no reads survived the average or single bp quality filtering.')



        info('Aligning sequences...')

        ####INITIALIZE CACHE####
        variantCache = {}

        #operates on variantCache
        alnStats = process_fastq(processed_output_filename,variantCache,ref_names,refs,args)

        info('Done!')

        info('Quantifying indels/substitutions...')

        #ANALYZE ALIGNMENTS

        ###initialize
        N_TOTAL = 0


        class_counts = {}

        counts_total = {}
        counts_modified = {}
        counts_unmodified = {}

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

        #we don't touch the exons at all, the read can be still modified though..
        counts_non_modified_non_frameshift = {}

        counts_splicing_sites_modified = {}

        ################
        class_counts = {} # number of reads in each class e.g. "ref1_UNMODIFIED" -> 50

        alleles_list = [] #will be turned into df with rows with information for each variant (allele)

        #for each reference, the following are computed individually
        all_insertion_count_vectors = {} #all insertions (including masked bases)
        all_insertion_left_count_vectors = {} #all insertions (including masked bases)
        all_deletion_count_vectors = {}
        all_substitution_count_vectors = {}
        all_indelsub_count_vectors = {}
        all_substitution_base_vectors = {}
        all_base_count_vectors = {} #number of times each base is seen

        insertion_count_vectors = {} #insertions that are in the target region
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

        ###iterate through variants
        for variant in variantCache:
            #skip variant if there were none observed
            variantCount = variantCache[variant][0]
            if (variantCount == 0):
                continue
            N_TOTAL += variantCount

            #create count for classes piechart
            class_names = [] #list of classes -- these get concatted later
            aln_reference_names = [] #list of references
            for payloadInd in range(1,len(variantCache[variant])):
                class_name = variantCache[variant][payloadInd]['closest_aln_name']+"_"+variantCache[variant][payloadInd]['classification']
                class_names.append(class_name)
                aln_reference_names.append(variantCache[variant][payloadInd]['closest_aln_name'])
            class_name = "&".join(class_names)
            aln_reference_name = "&".join(aln_reference_names)

            if class_name not in class_counts:
                    class_counts[class_name] = 0
            class_counts[class_name]+=variantCount

            #iterate through payloads -- if a read aligned equally-well to two references, it could have more than one payload
            for payloadInd in range(1,len(variantCache[variant])):
                variantPayload = variantCache[variant][payloadInd]
                refName = variantPayload['closest_aln_name']
                if refName not in ref_names: #got 'Ambiguous'
                    continue
                counts_total[refName] += variantCount
                if variantPayload['classification'] == 'MODIFIED':
                    counts_modified[refName] += variantCount
                else:
                    counts_unmodified[refName] += variantCount

                this_effective_len = variantPayload['closest_aln_length']


                this_has_insertions = False
                all_insertion_count_vectors[refName][variantPayload['all_insertion_positions']]+=variantCount
                all_insertion_left_count_vectors[refName][variantPayload['all_insertion_left_positions']]+=variantCount
                all_indelsub_count_vectors[refName][variantPayload['all_insertion_positions']]+=variantCount
                if not args.ignore_insertions:
                    inserted_n_lists[refName].extend([variantPayload['insertion_n']]*variantCount)
                    insertion_count_vectors[refName][variantPayload['insertion_positions']]+=variantCount
                    indelsub_count_vectors[refName][variantPayload['insertion_positions']]+=variantCount
                    this_effective_len = this_effective_len + variantPayload['insertion_n']
                    if variantPayload['insertion_n'] > 0:
                         counts_insertion[refName] += variantCount
                         this_has_insertions = True


                this_has_deletions = False
                all_deletion_count_vectors[refName][variantPayload['all_deletion_positions']]+=variantCount
                all_indelsub_count_vectors[refName][variantPayload['all_deletion_positions']]+=variantCount
                if not args.ignore_deletions:
                    deleted_n_lists[refName].extend([variantPayload['deletion_n']]*variantCount)
                    deletion_count_vectors[refName][variantPayload['deletion_positions']]+=variantCount
                    indelsub_count_vectors[refName][variantPayload['deletion_positions']]+=variantCount
                    this_effective_len = this_effective_len - variantPayload['deletion_n']
                    if variantPayload['deletion_n'] > 0:
                         counts_deletion[refName] += variantCount
                         this_has_deletions = True

                effective_len_lists[refName].extend([this_effective_len]*variantCount)

                this_has_substitutions = False
                all_substitution_count_vectors[refName][variantPayload['all_substitution_positions']] += variantCount
                all_indelsub_count_vectors[refName][variantPayload['all_substitution_positions']] += variantCount
                if not args.ignore_substitutions:
                    substituted_n_lists[refName].extend([variantPayload['substitution_n']] * variantCount)
                    substitution_count_vectors[refName][variantPayload['substitution_positions']] += variantCount
                    indelsub_count_vectors[refName][variantPayload['substitution_positions']] += variantCount
                    if variantPayload['substitution_n'] > 0:
                         counts_substitution[refName] += variantCount
                         this_has_substitutions = True

                    nucs = ['A','T','C','G','N']
                    for nuc in nucs:
                        isNuc = variantPayload['all_substitution_values'] == ord(nuc)
                        if(np.sum(isNuc) > 0):
                            locs = np.array(variantPayload['all_substitution_positions'])[isNuc]
                            all_substitution_base_vectors[refName + "_" + nuc ][locs] += variantCount


                if this_has_deletions:
                    if this_has_insertions:
                        if this_has_substitutions:
                            counts_insertion_and_deletion_and_substitution[refName] += variantCount
                        else:
                            counts_insertion_and_deletion[refName] += variantCount
                    else:
                        if this_has_substitutions:
                            counts_deletion_and_substitution[refName] += variantCount
                        else:
                            counts_only_deletion[refName] += variantCount
                else: #no deletions
                    if this_has_insertions:
                        if this_has_substitutions:
                            counts_insertion_and_substitution[refName] += variantCount
                        else:
                            counts_only_insertion[refName] += variantCount
                    else:
                        if this_has_substitutions:
                            counts_only_substitution[refName] += variantCount

                #set all_base_count_vectors
                aln_seq = variantPayload['aln_seq']
                ref_pos = variantPayload['ref_positions']
                for i in range(len(aln_seq)):
                    if ref_pos[i] < 0:
                        continue
                    nuc = aln_seq[i]
                    all_base_count_vectors[refName + "_" + nuc][ref_pos[i]] += variantCount

                exon_positions = refs[refName]['exon_positions']
                splicing_positions = refs[refName]['splicing_positions']
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



                length_modified_positions_exons=[]
                current_read_exons_modified = False
                current_read_spliced_modified = False

                for idx_ins,(ins_start,ins_end) in enumerate(insertion_coordinates):
                    insertion_length_vectors[refName][ins_start]+=(insertion_sizes[idx_ins]*variantCount)
                    insertion_length_vectors[refName][ins_end]+=(insertion_sizes[idx_ins]*variantCount)

                    if refs[refName]['contains_coding_seq']:
                        if set(exon_positions).intersection((ins_start, ins_end)): # check that we are inserting in one exon
                            set1 = set(exon_positions).intersection((ins_start, ins_end))
                            length_modified_positions_exons.append(insertion_sizes[idx_ins])
                            current_read_exons_modified = True

                for idx_del, (del_start,del_end) in enumerate(deletion_coordinates):
                    deletion_length_vectors[refName][range(del_start,del_end)] += (deletion_sizes[idx_del]*variantCount)

                if refs[refName]['contains_coding_seq']:
                    del_positions_to_append = sorted(set(exon_positions).intersection(set(deletion_positions)))
                    if del_positions_to_append:
                        # Always use the low include upper not
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
                        counts_splicing_sites_modified[refName] += variantCount

                    # if modified check if frameshift
                    if current_read_exons_modified:

                        if not length_modified_positions_exons:
                            # there are no indels
                            counts_modified_non_frameshift[ref_name] += variantCount
                            hists_inframe[refName][0] += variantCount
                        else:

                            effective_length = sum(length_modified_positions_exons)

                            if (effective_length % 3) == 0:
                                counts_modified_non_frameshift[ref_name] += variantCount
                                hists_inframe[refName][effective_length] += variantCount
                            else:
                                counts_modified_frameshift[ref_name] += variantCount
                                hists_frameshift[refName][effective_length] += variantCount

                    # the indels and subtitutions are outside the exon/s  so we don't care!
                    else:
                        counts_non_modified_non_frameshift[ref_name] += variantCount
                        insertion_count_vectors_noncoding[ref_name][insertion_positions] += variantCount
                        deletion_count_vectors_noncoding[ref_name][deletion_positions] += variantCount
                        substitution_count_vectors_noncoding[ref_name][substitution_positions] += variantCount

                alleleRow = {'#Reads':variantCount,
                           'Aligned_Sequence':variantPayload['aln_seq'],
                           'Reference_Sequence':variantPayload['aln_ref'],
                           'n_inserted':variantPayload['insertion_n'],
                           'n_deleted':variantPayload['deletion_n'],
                           'n_mutated':variantPayload['substitution_n'],
                           'Reference_Name':variantPayload['closest_aln_name'],
                           'Read_Status':variantPayload['classification'],
                           'Aligned_Reference_Names':aln_reference_name,
                           'Aligned_Reference_Scores':str(variantPayload['aln_scores']),
                           'ref_positions':variantPayload['ref_positions']
                }
                alleles_list.append(alleleRow)


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

        info('Done!')


        #order class_counts
        decorated_class_counts = []
        for class_count_name in class_counts:
            thisRefInd = 100
            thisIsMod = 1
            for idx,refName in enumerate(ref_names):
                if class_count_name.startswith(refName):
                    thisRefInd = idx
                    break
            if "UNMODIFIED" in class_count_name:
                thisIsMod = 0
            decorated_class_counts.append((thisRefInd,thisIsMod,class_count_name))
        decorated_class_counts.sort()
        class_counts_order = [class_count_name for thisRefInd,thisIsMod,class_count_name in decorated_class_counts]

        if N_TOTAL == 0:
            raise NoReadsAlignedException('Error: No alignments were found')

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


        all_insertion_pct_vectors = {} #all insertions/tot (including masked bases)
        all_deletion_pct_vectors = {}
        all_substitution_pct_vectors = {}
        all_indelsub_pct_vectors = {}

        insertion_pct_vectors = {} #insertions that are in the target region
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
                min_cut=min(refs[ref_name]['cut_points'])
                max_cut=max(refs[ref_name]['cut_points'])
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

        ref_info_file_name = _jp('CRISPResso_reference_info.txt')
        ref_info_file = open(ref_info_file_name,'w')
        refString = ( 'name' + "\t" +
            'sequence' + "\t" +
            'sequence_length' + "\t" +
            'min_aln_score' + "\t" +
            'cut_points' + "\t" +
            'gap_incentive' + "\t" +
            'offset_plots' + "\t" +
            'sgRNA_intervals' + "\t" +
            'sgRNA_sequences' + "\t" +
            'contains_guide' + "\t" +
            'contains_coding_seq' + "\t" +
            'exon_positions' + "\t" +
            'exon_intervals' + "\t" +
            'splicing_positions' + "\t" +
            'include_idxs' + "\t" +
            'exclude_idxs' + "\t" +
            'plot_idxs' + "\t" +
            'idx_cloned_from' + "\n")
        ref_info_file.write(refString)
        np.set_printoptions(linewidth=1000**1000) #no line breaks
        for ref_name in ref_names:
            refString = ( refs[ref_name]['name'] + "\t" +
                str(refs[ref_name]['sequence']) + "\t" +
                str(refs[ref_name]['sequence_length']) + "\t" +
                str(refs[ref_name]['min_aln_score']) + "\t" +
                str(refs[ref_name]['cut_points']) + "\t" +
                str(refs[ref_name]['gap_incentive']) + "\t" +
                str(refs[ref_name]['offset_plots']) + "\t" +
                str(refs[ref_name]['sgRNA_intervals']) + "\t" +
                str(refs[ref_name]['sgRNA_sequences']) + "\t" +
                str(refs[ref_name]['contains_guide']) + "\t" +
                str(refs[ref_name]['contains_coding_seq']) + "\t" +
                str(refs[ref_name]['exon_positions']) + "\t" +
                str(refs[ref_name]['exon_intervals']) + "\t" +
                str(refs[ref_name]['splicing_positions']) + "\t" +
                str(refs[ref_name]['include_idxs']) + "\t" +
                str(refs[ref_name]['exclude_idxs']) + "\t" +
                str(refs[ref_name]['plot_idxs']) + "\t" +
                str(refs[ref_name]['idx_cloned_from']) + "\n")
            ref_info_file.write(refString)
        ref_info_file.close()



        info('Making Plots...')
        ###############################################################################################################################################
        n_refs = len(ref_names)
        #helper function .. if there is only one reference, don't print the name on the top of every plot
        def get_plot_title_with_ref_name(plotTitle,refName):
            if n_refs > 1:
                return (plotTitle + ": " + refName)
            return plotTitle
        #(1)plot effective length
        for ref_name in ref_names:
            ref_len = refs[ref_name]['sequence_length']
            xmin = refs[ref_name]['xmin']
            xmax = refs[ref_name]['xmax']
            min_cut = refs[ref_name]['min_cut']
            max_cut = refs[ref_name]['max_cut']
            hdensity = refs[ref_name]['hdensity']
            hlengths = refs[ref_name]['hlengths']
            center_index = refs[ref_name]['center_index']

            plt.figure(figsize=(8.3,8))
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
            ax.set_yticklabels(['%.1f%% (%.0f)' % (pct,pct/100*N_TOTAL) for pct in y_label_values])
            ax.set_xlabel('Indel size (bp)')
            #lgd=plt.legend(['No indel','Indel'])
            lgd=ax.legend(['No indel','Indel'],loc='center', bbox_to_anchor=(0.5, -0.22),ncol=1, fancybox=True, shadow=True)
            lgd.legendHandles[0].set_height(3)
            lgd.legendHandles[1].set_height(3)

            plt.savefig(_jp('1.'+ref_name+'.Indel_size_distribution.pdf'),bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('1.'+ref_name+'.Indel_size_distribution.png'),bbox_inches='tight')
            plt.close()




        ###############################################################################################################################################

        ###############################################################################################################################################

        #(2) a piechart of classes
        labels = []
        sizes = []
        for class_name in class_counts_order:
            labels.append(class_name + "\n(" + str(class_counts[class_name]) + " reads)")
            sizes.append(100*class_counts[class_name]/float(N_TOTAL))

        fig=plt.figure(figsize=(12*1.5,12*1.5))
        ax = plt.subplot(111)
        patches, texts, autotexts =ax.pie(sizes, labels=labels,\
                        autopct='%1.1f%%')

        plt.axis('off')

        proptease = fm.FontProperties()
        proptease.set_size('x-large')
        plt.setp(autotexts, fontproperties=proptease)
        plt.setp(texts, fontproperties=proptease)
        plt.axis("equal")
        plt.savefig(_jp('2a.Alignment_Pie_Chart.pdf'),pad_inches=1,bbox_inches='tight')
        if args.save_also_png:
            plt.savefig(_jp('2a.Alignment_Pie_Chart.png'),pad_inches=1,bbox_inches='tight')
        plt.close()

        #(2b) a barchart of classes
        fig=plt.figure(figsize=(12*1.5,14.5*1.5))
        ax = plt.subplot(111)

        rects = ax.bar(np.arange(len(sizes)),sizes,color='silver')
        #label each bar
        for rect in rects:
            height = rect.get_height()
            if len(sizes) > 4:
                ax.text(rect.get_x() + rect.get_width()/2, height + 0.05,"%.2f%%"%height, ha='center', va='bottom',fontsize=12)
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
            plt.setp(ax.get_xticklabels(), fontsize=12, rotation='vertical',multialignment='right')
        plt.ylim(0,max(sizes)*1.1)
        plt.tight_layout()

        plt.savefig(_jp('2b.Alignment_Barplot.pdf'),pad_inches=1,bbox_inches='tight')
        if args.save_also_png:
            plt.savefig(_jp('2b.Alignment_Barplot.png'),pad_inches=1,bbox_inches='tight')
        plt.close()


        ###############################################################################################################################################


        ###############################################################################################################################################

        #(3) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);

        for ref_name in ref_names:
            y_values_mut = refs[ref_name]['y_values_mut']
            x_bins_mut = refs[ref_name]['x_bins_mut']
            y_values_ins = refs[ref_name]['y_values_ins']
            x_bins_ins = refs[ref_name]['x_bins_ins']
            y_values_del = refs[ref_name]['y_values_del']
            x_bins_del = refs[ref_name]['x_bins_del']

            n_this_category = counts_total[ref_name]
            if n_this_category < 1:
                continue

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

            plt.savefig(_jp('3.'+ref_name+'.Insertion_Deletion_Substitutions_size_hist.pdf'),bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('3.'+ref_name+'.Insertion_Deletion_Substitutions_size_hist.png'),bbox_inches='tight')
            plt.close()



        #(4) another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);

        #Indels location Plots
        for ref_name in ref_names:
            len_amplicon = refs[ref_name]['sequence_length']
            cut_points = refs[ref_name]['cut_points']
            offset_plots = refs[ref_name]['offset_plots']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
            include_idxs_list = list(refs[ref_name]['include_idxs'])

            n_this_category = counts_total[ref_name]
            if n_this_category < 1:
                continue

            n_this_category_modified = 0
            modifiedName = ref_name + "_MODIFIED"
            if modifiedName in class_counts:
                n_this_category_modified = class_counts[modifiedName]

            plt.figure(figsize=(10,10))

            y_max=max(all_indelsub_count_vectors[ref_name])*1.2

            #shade quantification window
            if len(include_idxs_list) > 1:
                lastStart = include_idxs_list[0]
                lastIdx = include_idxs_list[0]
                for idx in range(1,len(include_idxs_list)):
                    if include_idxs_list[idx] == lastIdx + 1:
                        lastIdx = include_idxs_list[idx]
                    else:
                        p = matplotlib.patches.Rectangle((lastStart,0),(lastIdx-lastStart)+1,y_max,facecolor='lightgray')
                        plt.gca().add_patch(p) #gca = get current axis
                        lastStart = include_idxs_list[idx]
                        lastIdx = include_idxs_list[idx]
                p = matplotlib.patches.Rectangle((lastStart,0),(lastIdx-lastStart)+1,y_max,facecolor='lightgray',label='Computed region')
                plt.gca().add_patch(p)

            plt.plot(all_indelsub_count_vectors[ref_name],'r',lw=3,label=get_plot_title_with_ref_name('Combined Insertions/Deletions/Substitutions',ref_name))
             #plt.hold(True)

            if cut_points:
                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                        plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                    else:
                        plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')


                for idx,sgRNA_int in enumerate(sgRNA_intervals):
                    if idx==0:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                    else:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


            lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.23),ncol=1, fancybox=True, shadow=True)
            ylabel_values = np.arange(0,1,1.0/6.0)
            if y_max > 0:
                y_label_values=np.arange(0,y_max,y_max/6.0)
            if len(ref_names) == 1:
                plt.ylabel('Sequences: % Total ( no. )')
                plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])
            else:
                plt.ylabel('Sequences: % Total ( % '+ref_name+', no. )')
                plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])
            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

            plt.title(get_plot_title_with_ref_name('Mutation position distribution',ref_name))
            plt.xlabel('Reference amplicon position (bp)')
            plt.ylim(0,max(1,y_max))
            plt.xlim(xmax=len_amplicon-1)
            plt.savefig(_jp('4a.'+ref_name+'.Combined_Insertion_Deletion_Substitution_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('4a.'+ref_name+'.Combined_Insertion_Deletion_Substitution_Locations.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
            plt.close()

#            print("subs: " + refName + ":"+ str(all_substitution_count_vectors[refName]))
            plt.figure(figsize=(10,10))

            #shade quantification window
            if len(include_idxs_list) > 1:
                lastStart = include_idxs_list[0]
                lastIdx = include_idxs_list[0]
                for idx in range(1,len(include_idxs_list)):
                    if include_idxs_list[idx] == lastIdx + 1:
                        lastIdx = include_idxs_list[idx]
                    else:
                        p = matplotlib.patches.Rectangle((lastStart,0),(lastIdx-lastStart)+1,y_max,facecolor='lightgray')
                        plt.gca().add_patch(p) #gca = get current axis
                        lastStart = include_idxs_list[idx]
                        lastIdx = include_idxs_list[idx]
                p = matplotlib.patches.Rectangle((lastStart,0),(lastIdx-lastStart)+1,y_max,facecolor='lightgray',label='Computed region')
                plt.gca().add_patch(p)

            plt.plot(all_insertion_count_vectors[refName],'r',lw=3,label='Insertions')
            #plt.hold(True)
            plt.plot(all_deletion_count_vectors[refName],'m',lw=3,label='Deletions')
            plt.plot(all_substitution_count_vectors[refName],'g',lw=3,label='Substitutions')

            y_max=max(max(all_insertion_count_vectors[refName]),max(all_deletion_count_vectors[refName]),max(all_substitution_count_vectors[refName]))*1.2


            if cut_points:
                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                        plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                    else:
                        plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',lw=2,label='_nolegend_')

                for idx,sgRNA_int in enumerate(sgRNA_intervals):
                    if idx==0:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                    else:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

            lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
            y_label_values = np.arange(0,1,1.0/6.0)
            if y_max > 0:
                y_label_values=np.arange(0,y_max,y_max/6.0)
            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

            plt.xlabel('Reference amplicon position (bp)')
            if len(ref_names) == 1:
                plt.ylabel('Sequences: % Total ( no. )')
                #plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(N_MODIFIED)*100, n_reads) for n_reads in y_label_values])
                plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])
            else:
                plt.ylabel('Sequences: % Total ( % '+ref_name+', no. )')
                plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])

            plt.ylim(0,max(1,y_max))
            plt.xlim(xmax=len_amplicon-1)

            plt.title(get_plot_title_with_ref_name('Mutation position distribution',ref_name))
            plt.savefig(_jp('4b.'+ref_name+'.Insertion_Deletion_Substitution_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('4b.'+ref_name+'.Insertion_Deletion_Substitution_Locations.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
            plt.close()

            plt.figure(figsize=(10,10))

            y_max=max(max(insertion_count_vectors[ref_name]),max(deletion_count_vectors[ref_name]),max(substitution_count_vectors[ref_name]),1)*1.2

            #shade quantification window
            if len(include_idxs_list) > 1:
                lastStart = include_idxs_list[0]
                lastIdx = include_idxs_list[0]
                for idx in range(1,len(include_idxs_list)):
                    if include_idxs_list[idx] == lastIdx + 1:
                        lastIdx = include_idxs_list[idx]
                    else:
                        p = matplotlib.patches.Rectangle((lastStart,0),(lastIdx-lastStart)+1,y_max,facecolor='lightgray')
                        plt.gca().add_patch(p) #gca = get current axis
                        lastStart = include_idxs_list[idx]
                        lastIdx = include_idxs_list[idx]
                p = matplotlib.patches.Rectangle((lastStart,0),(lastIdx-lastStart)+1,y_max,facecolor='lightgray',label='Computed region')
                plt.gca().add_patch(p)



            plt.plot(insertion_count_vectors[ref_name],'r',linewidth=3,label='Insertions')
            plt.plot(deletion_count_vectors[ref_name],'m',linewidth=3,label='Deletions')
            plt.plot(substitution_count_vectors[ref_name],'g',linewidth=3,label='Substitutions')


            if cut_points:
                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                        plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                    else:
                        plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='_nolegend_')


                for idx,sgRNA_int in enumerate(sgRNA_intervals):
                    if idx==0:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                    else:
                        plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')


            lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
            y_label_values = np.arange(0,1,1.0/6.0)
            if y_max > 0:
                y_label_values=np.arange(0,y_max,y_max/6.0).astype(int)
            plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])
            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

            plt.xlabel('Reference amplicon position (bp)')
            if len(ref_names) == 1:
                plt.ylabel('Sequences: % Total ( no. )')
                plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/float(N_TOTAL)*100,n_reads) for n_reads in y_label_values])
            else:
                plt.ylabel('Sequences: % Total ( % '+ref_name+', no. )')
                plt.yticks(y_label_values,['%.1f%% (%.1f%% , %d)' % (n_reads/float(N_TOTAL)*100,n_reads/float(n_this_category)*100, n_reads) for n_reads in y_label_values])


            plt.ylim(0,max(1,y_max))
            plt.xlim(xmax=len_amplicon-1)
            plt.title(get_plot_title_with_ref_name('Mutation position distribution',ref_name))
            plt.savefig(_jp('4c.'+ref_name+'.Masked_Insertion_Deletion_Substitution_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('4c.'+ref_name+'.Masked_Insertion_Deletion_Substitution_Locations.png'),bbox_extra_artists=(lgd,), bbox_inches='tight',pad=1)
            plt.close()

            #Position dependent indels plot
            fig=plt.figure(figsize=(24,10))
            ax1=fig.add_subplot(1,2,1)
            markerline, stemlines, baseline=ax1.stem(insertion_length_vectors[ref_name])
            plt.setp(markerline, 'markerfacecolor', 'r', 'markersize', 8)
            plt.setp(baseline, 'linewidth', 0)
            plt.setp(stemlines, 'color', 'r','linewidth',3)
            #plt.hold(True)
            y_max=max(insertion_length_vectors[ref_name])*1.2
            if cut_points:

                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                        ax1.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                    else:
                        ax1.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='_nolegend_')

            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
            plt.xlabel('Reference amplicon position (bp)')
            plt.ylabel('Average insertion length')
            plt.ylim(0,max(1,y_max))
            plt.xlim(xmax=len_amplicon-1)
            ax1.set_title(get_plot_title_with_ref_name('Position dependent insertion size',ref_name))
            plt.tight_layout()

            ax2=fig.add_subplot(1,2,2)
            markerline, stemlines, baseline=ax2.stem(deletion_length_vectors[ref_name])
            plt.setp(markerline, 'markerfacecolor', 'm', 'markersize', 8)
            plt.setp(baseline, 'linewidth', 0)
            plt.setp(stemlines, 'color', 'm','linewidth',3)
            #plt.hold(True)
            y_max=max(deletion_length_vectors[ref_name])*1.2
            if cut_points:

                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                        ax2.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                    else:
                        ax2.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='_nolegend_')

            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
            plt.xlabel('Reference amplicon position (bp)')
            plt.ylabel('Average deletion length')

            plt.ylim(ymin=0,ymax=max(1,y_max))
            plt.xlim(xmax=len_amplicon-1)
            ax2.set_title(get_plot_title_with_ref_name('Position dependent deletion size', ref_name))

            plt.tight_layout()


            plt.savefig(_jp('4d.'+ref_name+'.Position_dependent_average_indel_size.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('4d.'+ref_name+'.Position_dependent_average_indel_size.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close()


        ###############################################################################################################################################


        ###############################################################################################################################################
        #(5, 6) frameshift analyses plots
        for ref_name in ref_names:
            n_this_category = counts_total[ref_name]
            if n_this_category < 1:
                continue
            if (refs[ref_name]['contains_coding_seq']): #PERFORM FRAMESHIFT ANALYSIS
                #make frameshift plots
                len_amplicon = refs[ref_name]['sequence_length']
                cut_points = refs[ref_name]['cut_points']
                offset_plots = refs[ref_name]['offset_plots']
                sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
                MODIFIED_FRAMESHIFT = counts_modified_frameshift[ref_name]
                MODIFIED_NON_FRAMESHIFT = counts_modified_non_frameshift[ref_name]
                NON_MODIFIED_NON_FRAMESHIFT = counts_non_modified_non_frameshift[ref_name]
                SPLICING_SITES_MODIFIED = counts_splicing_sites_modified[ref_name]
                exon_intervals = refs[ref_name]['exon_intervals']
                count_total = counts_total[ref_name]
                count_modified = counts_modified[ref_name]
                count_unmodified = counts_unmodified[ref_name]

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
                ax2.plot([0,len_amplicon],[0,0],'-k',linewidth=2,label=ref_name+' sequence')
                #plt.hold(True)

                for idx,exon_interval in enumerate(exon_intervals):
                    if idx==0:
                        ax2.plot(exon_interval,[0,0],'-',linewidth=10,c=(0,0,1,0.5),label='Coding sequence/s',solid_capstyle='butt')
                    else:
                        ax2.plot(exon_interval,[0,0],'-',linewidth=10,c=(0,0,1,0.5),label='_nolegend_',solid_capstyle='butt')

                if cut_points:
                   ax2.plot(cut_points+offset_plots,np.zeros(len(cut_points)),'vr', ms=25,label='Predicted Cas9 cleavage site/s')

                plt.legend(bbox_to_anchor=(0, 0, 1., 0),  ncol=1, mode="expand", borderaxespad=0.,numpoints=1)
                plt.xlim(0,len_amplicon)
                plt.axis('off')

                proptease = fm.FontProperties()
                proptease.set_size('xx-large')
                plt.setp(autotexts, fontproperties=proptease)
                plt.setp(texts, fontproperties=proptease)
                plt.savefig(_jp('5.'+ref_name+'.Frameshift_In-frame_mutations_pie_chart.pdf'),pad_inches=1,bbox_inches='tight')
                if args.save_also_png:
                    plt.savefig(_jp('5.'+ref_name+'.Frameshift_In-frame_mutations_pie_chart.png'),pad_inches=1,bbox_inches='tight')
                plt.close()


                 #profiles-----------------------------------------------------------------------------------
                fig=plt.figure(figsize=(22,10))
                ax1=fig.add_subplot(2,1,1)
                x,y=map(np.array,zip(*[a for a in hists_frameshift[ref_name].iteritems()]))
                if sum(hists_frameshift[ref_name].values()) != 0:
                    y=y/float(sum(hists_frameshift[ref_name].values()))*100
                ax1.bar(x-0.5,y)
                ax1.set_xlim(-30.5,30.5)
                ax1.set_frame_on(False)
                ax1.set_xticks([idx for idx in range(-30,31) if idx % 3])
                ax1.tick_params(which='both',      # both major and minor ticks are affected
                   bottom='off',      # ticks along the bottom edge are off
                   top='off',         # ticks along the top edge are off
                   labelbottom='on') # labels along the bottom edge are off)
                ax1.yaxis.tick_left()
                xmin, xmax = ax1.get_xaxis().get_view_interval()
                ymin, ymax = ax1.get_yaxis().get_view_interval()
                ax1.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if idx % 3]],rotation='vertical')
                plt.title(get_plot_title_with_ref_name('Frameshift profile',ref_name))
                ax1.tick_params(axis='both', which='major', labelsize=32)
                ax1.tick_params(axis='both', which='minor', labelsize=32)
                plt.tight_layout()
                plt.ylabel('%')

                ax2=fig.add_subplot(2,1,2)
                x,y=map(np.array,zip(*[a for a in hists_inframe[ref_name].iteritems()]))
                if sum(hists_inframe[ref_name].values()) > 0:
                    y=y/float(sum(hists_inframe[ref_name].values()))*100
                ax2.bar(x-0.5,y,color=(0,1,1,0.2))
                ax2.set_xlim(-30.5,30.5)
                ax2.set_frame_on(False)
                ax2.set_xticks([idx for idx in range(-30,31) if (idx % 3 ==0) ])
                ax2.tick_params(which='both',      # both major and minor ticks are affected
                   bottom='off',      # ticks along the bottom edge are off
                   top='off',         # ticks along the top edge are off
                   labelbottom='on') # labels along the bottom edge are off)
                ax2.yaxis.tick_left()
                xmin, xmax = ax2.xaxis.get_view_interval()
                ymin, ymax = ax2.yaxis.get_view_interval()
                ax2.set_xticklabels([str(idx)  for idx in [idx for idx in range(-30,31) if (idx % 3==0)]],rotation='vertical')
                plt.title(get_plot_title_with_ref_name('In-frame profile',ref_name))
                plt.tight_layout()
                plt.ylabel('%')
                ax2.tick_params(axis='both', which='major', labelsize=32)
                ax2.tick_params(axis='both', which='minor', labelsize=32)
                plt.tight_layout()

                plt.savefig(_jp('6.'+ref_name+'.Frameshift_In-frame_mutation_profiles.pdf'),pad_inches=1,bbox_inches='tight')
                if args.save_also_png:
                    plt.savefig(_jp('6.'+ref_name+'.Frameshift_In-frame_mutation_profiles.png'),pad_inches=1,bbox_inches='tight')
                plt.close()

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
                plt.savefig(_jp('8.'+ref_name+'.Potential_Splice_Sites_pie_chart.pdf'),pad_inches=1,bbox_inches='tight')
                if args.save_also_png:
                    plt.savefig(_jp('8.'+ref_name+'.Potential_Splice_Sites_pie_chart.png'),pad_inches=1,bbox_inches='tight')
                plt.close()

                #non coding
                plt.figure(figsize=(10,10))
                plt.plot(insertion_count_vectors_noncoding[ref_name],'r',linewidth=3,label='Insertions')
                #plt.hold(True)
                plt.plot(deletion_count_vectors_noncoding[ref_name],'m',linewidth=3,label='Deletions')
                plt.plot(substitution_count_vectors_noncoding[ref_name],'g',linewidth=3,label='Substitutions')

                y_max=max(max(insertion_count_vectors_noncoding[ref_name]),max(deletion_count_vectors_noncoding[ref_name]),max(substitution_count_vectors_noncoding[ref_name]))*1.2

                if cut_points:

                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                                plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='Predicted cleavage position')
                        else:
                                plt.plot([cut_point+offset_plots[idx],cut_point+offset_plots[idx]],[0,y_max],'--k',linewidth=2,label='_nolegend_')

                        for idx,sgRNA_int in enumerate(sgRNA_intervals):
                            if idx==0:
                               plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='sgRNA',solid_capstyle='butt')
                            else:
                               plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],linewidth=10,c=(0,0,0,0.15),label='_nolegend_',solid_capstyle='butt')

                lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.28),ncol=1, fancybox=True, shadow=True)
                plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )

                plt.xlabel('Reference amplicon position (bp)')
                plt.ylabel('Sequences (no.)')
                plt.ylim(0,max(1,y_max))
                plt.xlim(xmax=len_amplicon-1)
                plt.title(get_plot_title_with_ref_name('Noncoding mutation position distribution',ref_name))
                plt.savefig(_jp('7.'+ref_name+'.Insertion_Deletion_Substitution_Locations_Noncoding.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                if args.save_also_png:
                    plt.savefig(_jp('7.'+ref_name+'.Insertion_Deletion_Substitution_Locations_Noncoding.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close()

        def count_alternate_alleles(sub_base_vectors,ref_name,ref_sequence,ref_total_aln_reads):
            #create vectors with all allele frequencies -- not just the substitution (the reference allele will not be 0)
            count_sub_base_vectors = {}
            count_sub_base_vectors['A'] = list(sub_base_vectors[ref_name+"_A"])
            count_sub_base_vectors['C'] = list(sub_base_vectors[ref_name+"_C"])
            count_sub_base_vectors['G'] = list(sub_base_vectors[ref_name+"_G"])
            count_sub_base_vectors['T'] = list(sub_base_vectors[ref_name+"_T"])
            count_sub_base_vectors['N'] = list(sub_base_vectors[ref_name+"_N"])

            #count the total number of times each substitution occurs
            alt_nuc_counts = {}
            alph = ['A','C','G','T','N']
            for a in alph:
                alt_nuc_counts[a] = {}
                for b in alph:
                    alt_nuc_counts[a][b] = 0

            for idx,c in enumerate(ref_sequence):
                tot_sub_at_idx = 0
                for a in alph:
                    sub = sub_base_vectors[ref_name+"_" + a][idx]
                    alt_nuc_counts[c][a] += sub
                    tot_sub_at_idx += sub

            df_subs = pd.DataFrame([count_sub_base_vectors["A"],count_sub_base_vectors["C"],count_sub_base_vectors["G"],count_sub_base_vectors["T"],count_sub_base_vectors["N"]])
            df_subs.index = ['A','C','G','T','N']
            df_subs.columns = list(ref_sequence)
            return (df_subs,alt_nuc_counts)

            ############

        ######PLOT substitution INFORMATION

        if not args.crispresso1_mode:
            for ref_name in ref_names:
                n_this_category = counts_total[ref_name]
                if n_this_category < 1:
                    continue

                ref_len = refs[ref_name]['sequence_length']
                ref_seq = refs[ref_name]['sequence']
                ref_include_idx = refs[ref_name]['include_idxs']
                ref_plot_idxs = refs[ref_name]['plot_idxs']


                tot_aln_reads = counts_total[ref_name]
                mask_ref_seq = [list(ref_seq)[x] for x in ref_include_idx]

                ##nucleotide counts
                df_nuc_freq = pd.DataFrame([base_count_vectors[ref_name+"_A"],base_count_vectors[ref_name+"_C"],base_count_vectors[ref_name+"_G"],base_count_vectors[ref_name+"_T"],base_count_vectors[ref_name+"_N"],base_count_vectors[ref_name+'_-']])
                df_nuc_freq.index = ['A','C','G','T','N','-']
                df_nuc_freq.columns = mask_ref_seq
                #print table showing nuc frequencies (sum to total alleles) (in target region)
                df_nuc_freq.to_csv(_jp(ref_name + '.target_nucleotide_frequency_table.txt'),sep='\t',header=True,index=True)

                df_nuc_pct = df_nuc_freq.divide(tot_aln_reads)
                df_nuc_pct.to_csv(_jp(ref_name + '.target_nucleotide_percentage_table.txt'),sep='\t',header=True,index=True)

                df_nuc_freq_all = pd.DataFrame([all_base_count_vectors[ref_name+"_A"],all_base_count_vectors[ref_name+"_C"],all_base_count_vectors[ref_name+"_G"],all_base_count_vectors[ref_name+"_T"],all_base_count_vectors[ref_name+"_N"],all_base_count_vectors[ref_name+'_-']])
                df_nuc_freq_all.index = ['A','C','G','T','N','-']
                df_nuc_freq_all.columns = list(ref_seq)
                #print table showing nuc frequencies (sum to total alleles) (in entire region)
                df_nuc_freq_all.to_csv(_jp(ref_name + '.nucleotide_frequency_table.txt'),sep='\t',header=True,index=True)

                df_nuc_pct_all = df_nuc_freq_all.divide(tot_aln_reads)
                df_nuc_pct_all.to_csv(_jp(ref_name + '.nucleotide_percentage_table.txt'),sep='\t',header=True,index=True)

                #substitution frequencies
                df_sub_freq,alt_nuc_counts = count_alternate_alleles(
                    sub_base_vectors = substitution_base_vectors,
                    ref_name = ref_name,
                    ref_sequence = mask_ref_seq,
                    ref_total_aln_reads = tot_aln_reads
                    )

                #print table showing sub frequencies
                df_sub_freq.to_csv(_jp(ref_name + '.Target_Substitution_Frequency_table.txt'),sep='\t',header=True,index=True)

                df_sub_freq_all,alt_nuc_counts_all = count_alternate_alleles(
                    sub_base_vectors = all_substitution_base_vectors,
                    ref_name = ref_name,
                    ref_sequence = ref_seq,
                    ref_total_aln_reads = tot_aln_reads
                    )


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
                nuc_df_for_plot.insert(0,'Batch',ref_name)
                mod_df_for_plot = modification_percentage_summary_df.copy()
                mod_df_for_plot.insert(0,'Batch',ref_name)
                CRISPRessoPlot.plot_nucleotide_quilt(nuc_df_for_plot,mod_df_for_plot,_jp('11.'+ref_name + '.nucleotide_percentage_quilt'),args.save_also_png,sgRNA_intervals=refs[ref_name]['sgRNA_intervals'])

                if args.base_editor_mode:
                    CRISPRessoPlot.plot_subs_across_ref(
                        ref_len = ref_len,
                        ref_seq = ref_seq,
                        ref_name = ref_name,
                        ref_count = tot_aln_reads,
                        all_substitution_base_vectors = all_substitution_base_vectors,

                        plot_title = get_plot_title_with_ref_name('Substitution frequency', ref_name),
                        fig_filename_root= _jp('10a.'+ref_name+'.Substitution_Frequencies'),
                        save_also_png = args.save_also_png
                        )

                    #plot all substitution rates in mask
                    CRISPRessoPlot.plot_sub_freqs(
                        alt_nuc_counts = alt_nuc_counts,
                        plot_title = get_plot_title_with_ref_name('Substitution frequency in mask', ref_name),
                        fig_filename_root = _jp('10b.'+ref_name+'.Target_Substitution_Frequency_Barplot'),
                        save_also_png = args.save_also_png
                        )

                    #plot all substitution rates in entire region
                    CRISPRessoPlot.plot_sub_freqs(
                        alt_nuc_counts = alt_nuc_counts_all,
                        plot_title = get_plot_title_with_ref_name('Substitution frequency in entire amplicon', ref_name),
                        fig_filename_root = _jp('10c.'+ref_name+'.Substitution_Frequency_Barplot'),
                        save_also_png = args.save_also_png
                        )

                #print table showing all nuc frequencies (sum to total alleles) (in entire region)
                df_sub_freq_all.to_csv(_jp(ref_name + '.substitution_frequency_table.txt'),sep='\t',header=True,index=True)

#                CRISPRessoPlot.plot_nuc_freqs(
#                    df_nuc_freq = df_nuc_freq,
#                    tot_aln_reads = tot_aln_reads,
#                    plot_title = get_plot_title_with_ref_name('Nucleotide Frequencies',ref_name),
#                    fig_filename_root = _jp('14a.'+ref_name+'.nucleotide_frequency'),
#                    save_also_png = args.save_also_png
#                    )

                CRISPRessoPlot.plot_log_nuc_freqs(
                    df_nuc_freq = df_nuc_freq,
                    tot_aln_reads = tot_aln_reads,
                    plot_title = get_plot_title_with_ref_name('Log2 Nucleotide Frequencies',ref_name),
                    fig_filename_root = _jp('14a.'+ref_name+'.log2_Nucleotide_Frequency'),
                    save_also_png = args.save_also_png
                    )

                plot_ref_seq = ''.join([ref_seq[i] for i in ref_plot_idxs])
                plot_nuc_pcts = df_nuc_pct_all.iloc[:,ref_plot_idxs]
                plot_nuc_freqs = df_nuc_freq_all.iloc[:,ref_plot_idxs]

                CRISPRessoPlot.plot_conversion_at_sel_nucs(
                    df_subs = plot_nuc_pcts,
                    ref_name = ref_name,
                    ref_sequence = plot_ref_seq,
                    plot_title = get_plot_title_with_ref_name('Substitution Frequencies at Target Nucleotides',ref_name),
                    conversion_nuc_from = args.conversion_nuc_from,
                    fig_filename_root = _jp('10f.'+ref_name+'.Target_Selected_Conversion'),
                    save_also_png = args.save_also_png
                    )
                CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_nosub(
                    df_subs = plot_nuc_pcts,
                    ref_name = ref_name,
                    ref_sequence = plot_ref_seq,
                    plot_title = get_plot_title_with_ref_name('Substitution Frequencies at Target Nucleotides',ref_name),
                    conversion_nuc_from = args.conversion_nuc_from,
                    fig_filename_root = _jp('10g.'+ref_name+'.Target_Selected_Conversion_no_ref'),
                    save_also_png = args.save_also_png
                    )
                from_nuc_indices = [pos for pos, char in enumerate(list(plot_nuc_pcts.columns.values)) if char == args.conversion_nuc_from]
                just_sel_nuc_pcts = plot_nuc_pcts.iloc[:,from_nuc_indices].copy() #only nucleotides targeted by base editing
                just_sel_nuc_pcts.columns = [char + str(pos+1) for pos,char in enumerate(list(just_sel_nuc_pcts.columns.values))]
                just_sel_nuc_freqs = plot_nuc_freqs.iloc[:,from_nuc_indices].copy()
                just_sel_nuc_freqs.columns = [char + str(pos+1) for pos,char in enumerate(list(just_sel_nuc_freqs.columns.values))]
                just_sel_nuc_pcts.to_csv(_jp(ref_name + '.target_selected_nucleotide_percentage_table.txt'),sep='\t',header=True,index=True)
                just_sel_nuc_freqs.to_csv(_jp(ref_name + '.target_selected_nucleotide_frequency_table.txt'),sep='\t',header=True,index=True)


        ##new plots alleles around cut_sites

        for ref_name in ref_names:
            n_this_category = counts_total[ref_name]
            if n_this_category < 1:
                continue

            sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
            cut_points = refs[ref_name]['cut_points']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']

            for sgRNA,cut_point in zip(sgRNA_sequences,cut_points):
                df_allele_around_cut=CRISPRessoShared.get_dataframe_around_cut(df_alleles.loc[df_alleles['Reference_Name'] == ref_name],cut_point,args.offset_around_cut_to_plot)

                    #write alleles table to file
                df_allele_around_cut.to_csv(_jp('%s.Alleles_frequency_table_around_cut_site_for_%s.txt' % (ref_name,sgRNA)),sep='\t',header=True)

                ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-args.offset_around_cut_to_plot+1:cut_point+args.offset_around_cut_to_plot+1]
                CRISPRessoPlot.plot_alleles_table(ref_seq_around_cut,df_alleles=df_allele_around_cut,fig_filename_root=_jp('9.%s.Alleles_frequency_table_around_cut_site_for_%s' % (ref_name,sgRNA)), MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot,MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot,SAVE_ALSO_PNG=args.save_also_png,base_editor_mode=args.base_editor_mode,sgRNA_intervals=sgRNA_intervals)

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

            if args.left_adapter_umi_trim_seq is not None:
                files_to_remove +=temp_flip_filename

            for file_to_remove in files_to_remove:
                try:
                    if os.path.islink(file_to_remove):
                        os.unlink(file_to_remove)
                    else:
                        os.remove(file_to_remove)
                except:
                    warn('Skipping:%s' %file_to_remove)

             #write effect vectors as plain text files
        info('Saving processed data...')
        def save_vector_to_file(vector,name):
            np.savetxt(_jp('%s.txt' %name), np.vstack([(np.arange(len(vector))+1),vector]).T, fmt=['%d','%.18e'],delimiter='\t', newline='\n', header='amplicon position\teffect',footer='', comments='# ')

        def save_count_vectors_to_file(vectors,vectorNames,refSeq,name):
            outfile = open(_jp('%s.txt'%name),"w")
            outfile.write("Sequence\t"+"\t".join(list(refSeq))+"\n") #first row: reference sequence
            for vector,vectorName in zip(vectors,vectorNames):
                outfile.write(vectorName +"\t" + "\t".join([str(x) for x in vector]) + "\n") #next, vectors are printed
            outfile.close()

        #write alleles table
        #crispresso1Cols = ["Aligned_Sequence","Reference_Sequence","NHEJ","UNMODIFIED","HDR","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        #df_alleles.ix[:,crispresso1Cols].to_csv(_jp('Alleles_frequency_table.txt'),sep='\t',header=True,index=None)
        #crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
#        crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads","Aligned_Reference_Names","Aligned_Reference_Scores"]
        crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        df_alleles.ix[:,crispresso2Cols].to_csv(_jp('Alleles_frequency_table.txt'),sep='\t',header=True,index=None)



        with open(_jp('Quantification_of_editing_frequency.txt'),'w+') as outfile:
            outfile.write("Quantification of editing frequency:\n")
            for ref_name in ref_names:
                n_unmod = counts_unmodified[ref_name]
                n_mod = counts_modified[ref_name]

                n_insertion = counts_insertion[ref_name]
                n_deletion = counts_deletion[ref_name]
                n_substitution = counts_substitution[ref_name]

                outfile.write("%s: Unmodified: %d Modified: %d\n" % (ref_name,n_unmod,n_mod))
                outfile.write("(%d reads with insertions, %d reads with deletions, %d reads with substitutions)\n" % (n_insertion,n_deletion,n_substitution))

            outfile.write('Total Aligned:%d reads ' % N_TOTAL)

        with open(_jp('CRISPResso_quantification_of_editing_frequency.txt'),'w+') as outfile:
            outfile.write('Reference\tTotal\tUnmodified\tModified\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions\n')
            for ref_name in ref_names:
                n_tot = counts_total[ref_name]
                n_unmod = counts_unmodified[ref_name]
                n_mod = counts_modified[ref_name]

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

                vals = [ref_name]
                vals.extend([str(x) for x in [n_tot,n_unmod,n_mod,n_insertion,n_deletion,n_substitution,n_only_insertion,n_only_deletion,n_only_substitution,n_insertion_and_deletion,n_insertion_and_substitution,n_deletion_and_substitution,n_insertion_and_deletion_and_substitution]])
                outfile.write("\t".join(vals) + "\n")



        #write statistics
        with open(_jp('Mapping_statistics.txt'),'w+') as outfile:
            outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d\n' % (N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_TOTAL))

        with open(_jp('CRISPResso_mapping_statistics.txt'),'w+') as outfile:
            outfile.write('READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n')
            outfile.write("\t".join([str(x) for x in[N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_TOTAL,alnStats['N_COMPUTED_ALN'],alnStats['N_CACHED_ALN'],alnStats['N_COMPUTED_NOTALN'],alnStats['N_CACHED_NOTALN']]]) + "\n")


        for ref_name in ref_names:

            n_this_category = counts_total[ref_name]
            if n_this_category < 1:
                continue

            save_vector_to_file(insertion_pct_vectors[ref_name],ref_name+'.effect_vector_insertion')
            save_vector_to_file(deletion_pct_vectors[ref_name],ref_name+'.effect_vector_deletion')
            save_vector_to_file(substitution_pct_vectors[ref_name],ref_name+'.effect_vector_substitution')
            save_vector_to_file(indelsub_pct_vectors[ref_name],ref_name+'.effect_vector_combined')

            #save mods in computation window
            save_count_vectors_to_file([insertion_count_vectors[ref_name],
                        deletion_count_vectors[ref_name],
                        substitution_count_vectors[ref_name],
                        indelsub_count_vectors[ref_name],
                        [counts_total[ref_name]]*refs[ref_name]['sequence_length']],
                        ['Insertions','Deletions','Substitutions','All_modifications','Total'],
                            refs[ref_name]['sequence'],ref_name+'.target_modification_count_vectors')

            #save all mods
            save_count_vectors_to_file([all_insertion_count_vectors[ref_name],
                        all_insertion_left_count_vectors[ref_name],
                        all_deletion_count_vectors[ref_name],
                        all_substitution_count_vectors[ref_name],
                        all_indelsub_count_vectors[ref_name],
                        [counts_total[ref_name]]*refs[ref_name]['sequence_length']],
                        ['Insertions','Insertions_Left','Deletions','Substitutions','All_modifications','Total'],
                            refs[ref_name]['sequence'],ref_name+'.modification_count_vectors')

            if (refs[ref_name]['contains_coding_seq']): #PERFORM FRAMESHIFT ANALYSIS
                MODIFIED_FRAMESHIFT = counts_modified_frameshift[ref_name]
                MODIFIED_NON_FRAMESHIFT = counts_modified_non_frameshift[ref_name]
                NON_MODIFIED_NON_FRAMESHIFT = counts_non_modified_non_frameshift[ref_name]
                SPLICING_SITES_MODIFIED = counts_splicing_sites_modified[ref_name]
                with open(_jp(ref_name+'.Frameshift_analysis.txt'),'w+') as outfile:
                        outfile.write('Frameshift analysis:\n\tNoncoding mutation:%d reads\n\tIn-frame mutation:%d reads\n\tFrameshift mutation:%d reads\n' %(NON_MODIFIED_NON_FRAMESHIFT, MODIFIED_NON_FRAMESHIFT ,MODIFIED_FRAMESHIFT))

                with open(_jp('Splice_sites_analysis.txt'),'w+') as outfile:
                        outfile.write('Splice sites analysis:\n\tUnmodified:%d reads\n\tPotential splice sites modified:%d reads\n' %(n_reads- SPLICING_SITES_MODIFIED, SPLICING_SITES_MODIFIED))

                save_vector_to_file(insertion_pct_vectors_noncoding[ref_name],ref_name+'.effect_vector_insertion_noncoding')
                save_vector_to_file(deletion_pct_vectors_noncoding[ref_name],ref_name+'.effect_vector_deletion_noncoding')
                save_vector_to_file(substitution_pct_vectors_noncoding[ref_name],ref_name+'.effect_vector_substitution_noncoding')

            len_amplicon = refs[ref_name]['sequence_length']
            cut_points = refs[ref_name]['cut_points']
            offset_plots = refs[ref_name]['offset_plots']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
            if cut_points:
                cp.dump(sgRNA_intervals, open( _jp(ref_name+'.sgRNA_intervals.pickle'), 'wb' ) )

            if sgRNA_intervals:
                cp.dump( cut_points, open( _jp(ref_name+'.cut_points.pickle'), 'wb' ) )

            if offset_plots.any():
                cp.dump(offset_plots,open( _jp(ref_name+'.offset_plots.pickle'), 'wb' ) )

            hdensity = refs[ref_name]['hdensity']
            hlengths = refs[ref_name]['hlengths']
            center_index = refs[ref_name]['center_index']

            y_values_mut = refs[ref_name]['y_values_mut']
            x_bins_mut = refs[ref_name]['x_bins_mut']
            y_values_ins = refs[ref_name]['y_values_ins']
            x_bins_ins = refs[ref_name]['x_bins_ins']
            y_values_del = refs[ref_name]['y_values_del']
            x_bins_del = refs[ref_name]['x_bins_del']

            pd.DataFrame(np.vstack([hlengths,hdensity]).T,columns=['indel_size','fq']).to_csv(_jp(ref_name+'.indel_histogram.txt'),index=None,sep='\t')
            pd.DataFrame(np.vstack([x_bins_ins[:-1],y_values_ins]).T,columns=['ins_size','fq']).to_csv(_jp(ref_name+'.insertion_histogram.txt'),index=None,sep='\t')
            pd.DataFrame(np.vstack([-x_bins_del[:-1],y_values_del]).T,columns=['del_size','fq']).to_csv(_jp(ref_name+'.deletion_histogram.txt'),index=None,sep='\t')
            pd.DataFrame(np.vstack([x_bins_mut[:-1],y_values_mut]).T,columns=['sub_size','fq']).to_csv(_jp(ref_name+'.substitution_histogram.txt'),index=None,sep='\t')


            if args.dump:
                info('Dumping all the processed data...')
                np.savez(_jp('%s.effect_vector_insertion'%ref_name),insertion_pct_vectors[ref_name])
                np.savez(_jp('%s.effect_vector_deletion'%ref_name),deletion_pct_vectors[ref_name])
                np.savez(_jp('%s.effect_vector_substitution'%ref_name),substitution_pct_vectors[ref_name])

                np.savez(_jp('%s.effect_vector_combined'%ref_name),indelsub_pct_vectors[ref_name])

                np.savez(_jp('%s.position_dependent_vector_avg_insertion_size'%ref_name),insertion_length_vectors[ref_name])
                np.savez(_jp('%s.position_dependent_vector_avg_deletion_size'%ref_name),deletion_length_vectors[ref_name])



        info('All Done!')
        print(CRISPRessoShared.get_crispresso_footer())

        sys.exit(0)



    except NTException as e:
        print_stacktrace_if_debug()
        error('Alphabet error, please check your input.\n\nERROR: %s' % e)
        sys.exit(1)
    except SgRNASequenceException as e:
        print_stacktrace_if_debug()
        error('sgRNA error, please check your input.\n\nERROR: %s' % e)
        sys.exit(2)

    except TrimmomaticException as e:
        print_stacktrace_if_debug()
        error('Trimming error, please check your input.\n\nERROR: %s' % e)
        sys.exit(4)
    except FlashException as e:
        print_stacktrace_if_debug()
        error('Merging error, please check your input.\n\nERROR: %s' % e)
        sys.exit(5)
    except BadParameterException as e:
        print_stacktrace_if_debug()
        error('Parameter error, please check your input.\n\nERROR: %s' % e)
        sys.exit(6)
    except NoReadsAlignedException as e:
        print_stacktrace_if_debug()
        error('Alignment error, please check your input.\n\nERROR: %s' % e)
        sys.exit(7)
    except AutoException as e:
        print_stacktrace_if_debug()
        error('Autorun error. This sample cannot be run in auto mode.\n\nERROR: %s' % e)
        sys.exit(8)

    except ExonSequenceException as e:
        print_stacktrace_if_debug()
        error('Coding sequence error, please check your input.\n\nERROR: %s' % e)
        sys.exit(11)
    except DuplicateSequenceIdException as e:
        print_stacktrace_if_debug()
        error('Fastq file error, please check your input.\n\nERROR: %s' % e)
        sys.exit(12)
    except NoReadsAfterQualityFiltering as e:
        print_stacktrace_if_debug()
        error('Filtering error, please check your input.\n\nERROR: %s' % e)
        sys.exit(13)
    except Exception as e:
        print_stacktrace_if_debug()
        error('Unexpected error, please check your input.\n\nERROR: %s' % e)
        sys.exit(-1)
