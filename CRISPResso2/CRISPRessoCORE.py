#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''

import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import argparse
from collections import defaultdict
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, wait
import errno
import gzip
import json
import zipfile
import os
import re
import subprocess as sb
import traceback


from CRISPResso2 import CRISPRessoCOREResources
from CRISPResso2 import CRISPRessoReport
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoPlot
from CRISPResso2 import CRISPResso2Align
from CRISPResso2 import CRISPRessoMultiProcessing

from datetime import datetime
present = datetime.now()
#d1 = datetime.strptime('21/07/2019','%d/%m/%Y')
#if present > d1:
#    print('\nYour version of CRISPResso2 is out of date. Please download a new version.\n')
#    sys.exit(1)

import logging

logging.basicConfig(
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info

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
     p = sb.Popen(cmd, shell=True, stdout=sb.PIPE)
     return int(p.communicate()[0].strip())

def get_n_reads_fastq(fastq_filename):
    p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < \"%s\" | wc -l" % fastq_filename, shell=True, stdout=sb.PIPE)
    return int(float(p.communicate()[0])/4.0)

def get_n_reads_bam(bam_filename,bam_chr_loc=""):
    cmd = "samtools view -c " + bam_filename + " " + bam_chr_loc
    p = sb.Popen(cmd, shell=True, stdout=sb.PIPE)
    try:
        retval = int(float(p.communicate()[0]))
    except ValueError:
        raise CRISPRessoShared.InstallationException('Error when running the command:' + cmd + '\nCheck that samtools is installed correctly.')
    return retval

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

def get_pe_scaffold_search(prime_edited_ref_sequence, prime_editing_pegRNA_extension_seq, prime_editing_pegRNA_scaffold_seq, prime_editing_pegRNA_scaffold_min_match_length):
    """
    For prime editing, determines the scaffold string to search for (the shortest substring of args.prime_editing_pegRNA_scaffold_seq not in the prime-edited reference sequence)
    params:
     prime_edited_ref_sequence: reference sequence of the prime-edited sequence
     prime_editing_extension_seq: RNA sequence of extension sequence
     prime_editing_pegRNA_scaffold_seq: sequence of the scaffold sequence
     prime_editing_pegRNA_scaffold_min_match_length: minimum number of bases required to match between scaffold and read to count as scaffold-incorporated
    returns:
     tuple of(
         index of location in ref to find scaffold seq if it exists
         shortest dna sequence to identify scaffold sequence
     )
    """
    info('Processing pegRNA scaffold sequence...')
    #first, define the sequence we are looking for (extension plus the first base(s) of the scaffold)
    scaffold_dna = CRISPRessoShared.reverse_complement(prime_editing_pegRNA_scaffold_seq)

    extension_seq_dna_top_strand = prime_editing_pegRNA_extension_seq.upper().replace('U', 'T')
    prime_editing_extension_seq_dna = CRISPRessoShared.reverse_complement(extension_seq_dna_top_strand)

    scaffold_start_loc = prime_edited_ref_sequence.index(prime_editing_extension_seq_dna)+len(prime_editing_extension_seq_dna)

    #next find min length scaffold that when combined with min extension sequence is not in edited sequence
    len_scaffold_to_use = prime_editing_pegRNA_scaffold_min_match_length # length of the scaffold sequence to be added to the extension sequence, start at 1bp
    scaffold_dna_search = prime_editing_extension_seq_dna + scaffold_dna[0:len_scaffold_to_use]
    while scaffold_dna_search in prime_edited_ref_sequence:
        if len_scaffold_to_use > len(scaffold_dna):
            raise CRISPRessoShared.BadParameterException('The DNA scaffold provided is found in the unedited reference sequence. Please provide a longer scaffold sequence.')
        len_scaffold_to_use += 1
        scaffold_dna_search = prime_editing_extension_seq_dna + scaffold_dna[0:len_scaffold_to_use]

    info('Searching for scaffold-templated reads with the sequence: \'' + str(scaffold_dna[0:len_scaffold_to_use]) +'\' starting at position '+ str(scaffold_start_loc) + ' in reads that align to the prime-edited sequence')
    return (scaffold_start_loc, scaffold_dna[0:len_scaffold_to_use])

def get_new_variant_object(args, fastq_seq, refs, ref_names, aln_matrix, pe_scaffold_dna_info):
    """
    Gets the payload object for a read that hasn't been seen in the cache yet
    params:
     args: CRISPResso2 args
     fastq_seq: read sequence to align
     refs: dict with info for all refs
     ref_names: list of ref names
     aln_matrix: alignment matrix for needleman wunsch
     pe_scaffold_dna_info: for prime-editing tuple of(
         index of location in ref to find scaffold seq if it exists
         shortest dna sequence to identify scaffold sequence
         )

    returns:
     variant payload
    """

    aln_scores = []
    best_match_score = -1
    best_match_s1s = []
    best_match_s2s = []
    best_match_names = []
    best_match_strands = []
    ref_aln_details = []
    for idx, ref_name in enumerate(ref_names):
        #get alignment and score from cython
        #score = 100 * #matchedBases / length(including gaps)
        seed_i = 0
        found_forward_count = 0
        found_reverse_count = 0
        aln_strand = '+'
        while seed_i < args.aln_seed_count and seed_i < len(refs[ref_name]['fw_seeds']):
            if refs[ref_name]['fw_seeds'][seed_i] in fastq_seq: #is forward
                found_forward_count += 1
            if refs[ref_name]['rc_seeds'][seed_i] in fastq_seq: #is rc
                found_reverse_count += 1
            seed_i += 1
        if found_forward_count > args.aln_seed_min and found_reverse_count == 0:
            fws1, fws2, fwscore=CRISPResso2Align.global_align(fastq_seq, refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            s1 = fws1
            s2 = fws2
            score = fwscore
        elif found_forward_count == 0 and found_reverse_count > args.aln_seed_min:
            rvs1, rvs2, rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            s1 = rvs1
            s2 = rvs2
            score = rvscore
            aln_strand = '-'
        else:
            fws1, fws2, fwscore=CRISPResso2Align.global_align(fastq_seq, refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            rvs1, rvs2, rvscore=CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq_seq), refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            s1 = fws1
            s2 = fws2
            score = fwscore
            if (rvscore > fwscore):
                s1 = rvs1
                s2 = rvs2
                score = rvscore
                aln_strand = '-'

#                print "for " + ref_name + " got fws1: " + str(fws1) + " and fws2: " + str(fws2) + " score: " +str(fwscore)
        aln_scores.append(score)
        ref_aln_details.append((ref_name, s1, s2, score))

        #reads are matched to the reference to which they best align. The 'min_aln_score' is calculated using only the changes in 'include_idxs'
        if score > best_match_score and score > refs[ref_name]['min_aln_score']:
            best_match_score = score
            best_match_s1s = [s1]
            best_match_s2s = [s2]
            best_match_names = [ref_name]
            best_match_strands = [aln_strand]
        elif score == best_match_score:
            best_match_s1s.append(s1)
            best_match_s2s.append(s2)
            best_match_names.append(ref_name)
            best_match_strands.append(aln_strand)

    if best_match_score > 0:
        new_variant = {}
        new_variant['count'] = 1
        new_variant['aln_ref_names'] = best_match_names
        new_variant['aln_scores'] = aln_scores
        new_variant['ref_aln_details'] = ref_aln_details
        new_variant['best_match_score'] = best_match_score
        class_names = []

        for idx in range(len(best_match_names)):
            best_match_name = best_match_names[idx]

            if args.use_legacy_insertion_quantification:
                payload = CRISPRessoCOREResources.find_indels_substitutions_legacy(best_match_s1s[idx], best_match_s2s[idx], refs[best_match_name]['include_idxs'])
            else:
                payload = CRISPRessoCOREResources.find_indels_substitutions(best_match_s1s[idx], best_match_s2s[idx], refs[best_match_name]['include_idxs'])

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
            payload['aln_strand'] = best_match_strands[idx]

            new_variant['variant_'+best_match_name] = payload

        new_variant['class_name'] = "&".join(class_names)

    else:
        new_variant = {}
        new_variant['count'] = 1
        new_variant['aln_scores'] = aln_scores
        new_variant['ref_aln_details'] = ref_aln_details
        new_variant['best_match_score'] = best_match_score
        return new_variant #return new variant with best match score of 0, but include the scores of insufficient alignments

    #handle ambiguous alignments
    if len(best_match_names) > 1:
        if args.assign_ambiguous_alignments_to_first_reference: #if ambiguous, and this flag is set, just assign to the first amplicon
            new_variant['class_name'] = class_names[0]
            new_variant['aln_ref_names'] = [best_match_names[0]]
        elif not args.expand_ambiguous_alignments: #got 'Ambiguous' -- don't count toward total (e.g. indels at each position for the ref)
            new_variant['class_name'] = 'AMBIGUOUS'

    #search for prime editing scaffold
    #if it's there, assign this to be assigned to the reference '"Scaffold-incorporated"'
    if args.prime_editing_pegRNA_scaffold_seq and 'Prime-edited' in best_match_names: #any scaffold extensions must be closer to the prime-edited sequence
        pe_read_possible_scaffold_loc = new_variant['variant_Prime-edited']['ref_positions'].index(pe_scaffold_dna_info[0]-1) + 1
        if new_variant['variant_Prime-edited']['aln_seq'][pe_read_possible_scaffold_loc:(pe_read_possible_scaffold_loc+len(pe_scaffold_dna_info[1]))] == pe_scaffold_dna_info[1]:
#            print('comparingHERE ' + new_variant['variant_Prime-edited']['aln_seq'][pe_read_possible_scaffold_loc:(pe_read_possible_scaffold_loc+len(pe_scaffold_dna_info[1])+5)] + ' from ' + new_variant['variant_Prime-edited']['aln_seq'] + ' and ' + new_variant['variant_Prime-edited']['aln_ref'])
            new_variant['aln_ref_names'] = ["Scaffold-incorporated"]
            new_variant['class_name'] = "Scaffold-incorporated"
            old_payload = deepcopy(new_variant['variant_Prime-edited']) #keep prime-edited allele and alignment
            old_payload['ref_name'] = "Scaffold-incorporated"
            new_variant['variant_'+"Scaffold-incorporated"] = old_payload

    return new_variant

def get_greater_qual_nuc(nuc1,qual1,nuc2,qual2):
    """
    Gets the base with the higher quality if nuc1 != nuc2, otherwise nuc1 if they are the same

    params:
    nuc1: nucleotide 1
    quality1: quality of nuc1 in phred score
    nuc2: nucleotide 2
    quality2: quality of nuc2 in phred score

    returns
    nuc: the nucleotide with the greater quality if nuc1 != nuc2
    """
    if nuc1 == nuc2:
        return nuc1
    elif ord(qual1) >= ord(qual2):
        return nuc1
    else:
        return nuc2

def get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2):
    """
    Creates a consensus alignment from alignments of two reads

    params:
    aln1_seq: read alignment of r1
    aln1_ref: ref alignment of r1
    qual1: quality scores of bases in r1
    aln2_seq: read alignment of r2
    aln2_ref: ref alignment of r2
    qual2: quality scores of bases in r2

    returns:
    aln_seq: consensus read alignment
    ref_seq: consensus ref alignment
    score: alignment homology score
    """

    # three sets of indices
    aln_ind_in_r1 = 0  # indexes in r1 aln1_seq (including gaps)
    aln_ind_in_r2 = 0
    ind_in_r1 = 0  # indexes in non-gapped positions in r1 aln1_seq to match to quailty
    ind_in_r2 = 0
    ref_ind_in_r1 = 0  # indexes in non-gapped positions in r1 ref aln1_ref to match to r2 ref aln2_ref
    ref_ind_in_r2 = 0

    #positions where we start and stop having information for each read
    ind_r1_start = len(aln1_seq) - len(aln1_seq.lstrip("-"))
    ind_r2_start = len(aln2_seq) - len(aln2_seq.lstrip("-"))
    ind_r1_stop = len(aln1_seq.rstrip("-")) - 1
    ind_r2_stop = len(aln2_seq.rstrip("-")) - 1

    final_aln = ""
    final_ref = ""


    # iterate over all positions
    # aln1_seq and aln1_ref have the same length, so it doesn't matter which one we choose
    while aln_ind_in_r1 < len(aln1_ref) or aln_ind_in_r2 < len(aln2_ref):
#        print('aln_ind_in_r1: ' + str(aln_ind_in_r1) + ' aln_ind_in_r2: ' + str(aln_ind_in_r2))
#        print('ind_in_r1: ' + str(ind_in_r1) + ' ind_in_r1: ' + str(ind_in_r2))
        if aln_ind_in_r1 >= ind_r1_start and aln_ind_in_r1 <= ind_r1_stop:
            if aln_ind_in_r2 >= ind_r2_start and aln_ind_in_r2 <= ind_r2_stop:
                #both are in range
                this_nuc = get_greater_qual_nuc(aln1_seq[aln_ind_in_r1], qual1[ind_in_r1], aln2_seq[aln_ind_in_r2], qual2[ind_in_r2])
                final_aln += this_nuc
                final_ref += aln1_ref[ref_ind_in_r1]
            else:
                #only r1 is in range
                this_nuc = aln1_seq[aln_ind_in_r1]
                final_aln += this_nuc
                final_ref += aln1_ref[ref_ind_in_r1]
        elif aln_ind_in_r2 >= ind_r2_start and aln_ind_in_r2 <= ind_r2_stop:
            #only r2 is in range
            this_nuc = aln2_seq[aln_ind_in_r2]
            final_aln += this_nuc
            final_ref += aln1_ref[ref_ind_in_r1]
        else:
            final_aln += "N"
            final_ref += aln1_ref[ref_ind_in_r1]

        if aln_ind_in_r1 < len(aln1_seq):
            if aln1_seq[aln_ind_in_r1] != '-':
                ind_in_r1 += 1
            if aln1_ref[aln_ind_in_r1] != '-':
                ref_ind_in_r1 += 1
        aln_ind_in_r1 += 1

        if aln_ind_in_r2 < len(aln2_seq):
            if aln2_seq[aln_ind_in_r2] != '-':
                ind_in_r2 += 1
            if aln2_ref[aln_ind_in_r2] != '-':
                ref_ind_in_r2 += 1
        aln_ind_in_r2 += 1



    final_homology_score = 0
    for i in range(len(final_ref)):
        if final_ref[i] == final_aln[i]:
            final_homology_score += 1

    return final_aln, final_ref, int(100*final_homology_score/float(len(final_ref)))


def get_new_variant_object_from_paired(args, fastq1_seq, fastq2_seq, fastq1_qual, fastq2_qual, refs, ref_names, aln_matrix, pe_scaffold_dna_info):
    """
    Gets the payload object for a read that hasn't been seen in the cache yet
    This operates for paired-end reads by aligning the reads, then if there is a discrepancy at any position, the value from the higher-quality read is used
    params:
     args: CRISPResso2 args
     fastq1_seq: r1 read sequence to align
     fastq2_seq: r2 read sequence to align
     fastq1_qual: r1 read quality
     fastq2_qual: r2 read quality
     refs: dict with info for all refs
     ref_names: list of ref names
     aln_matrix: alignment matrix for needleman wunsch
     pe_scaffold_dna_info: for prime-editing tuple of(
         index of location in ref to find scaffold seq if it exists
         shortest dna sequence to identify scaffold sequence
         )

    returns:
     variant payload
    """

    aln_scores = []
    best_match_score = -1
    best_match_s1s = []
    best_match_s2s = []
    best_match_names = []
    ref_aln_details = []
    for idx, ref_name in enumerate(ref_names):
        #get alignment and score from cython
        #score = 100 * #matchedBases / length(including gaps)
        seed_i = 0
        found_forward_count = 0
        found_reverse_count = 0
        while seed_i < args.aln_seed_count and seed_i < len(refs[ref_name]['fw_seeds']):
            if refs[ref_name]['fw_seeds'][seed_i] in fastq1_seq or refs[ref_name]['fw_seeds'][seed_i] in fastq2_seq: #is forward
                found_forward_count += 1
            if refs[ref_name]['rc_seeds'][seed_i] in fastq1_seq or refs[ref_name]['rc_seeds'][seed_i] in fastq2_seq: #is rc
                found_reverse_count += 1
            seed_i += 1
        if found_forward_count > args.aln_seed_min and found_reverse_count == 0:
            r1_fws1, r1_fws2, r1_fwscore = CRISPResso2Align.global_align(fastq1_seq, refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            r2_fws1, r2_fws2, r2_fwscore = CRISPResso2Align.global_align(fastq2_seq, refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            s1, s2, score = get_consensus_alignment_from_pairs(r1_fws1, r1_fws2, fastq1_qual, r2_fws1, r2_fws2, fastq2_qual)
        elif found_forward_count == 0 and found_reverse_count > args.aln_seed_min:
            r1_rvs1, r1_rvs2, r1_rvscore = CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq1_seq), refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            r2_rvs1, r2_rvs2, r2_rvscore = CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq2_seq), refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            s1, s2, score = get_consensus_alignment_from_pairs(r1_rvs1, r1_rvs2, fastq1_qual, r2_rvs1, r2_rvs2, fastq2_qual)
            s1 = rvs1
            s2 = rvs2
            score = rvscore
        else:
            r1_fws1, r1_fws2, r1_fwscore = CRISPResso2Align.global_align(fastq1_seq, refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            r2_fws1, r2_fws2, r2_fwscore = CRISPResso2Align.global_align(fastq2_seq, refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            fws1, fws2, fwscore = get_consensus_alignment_from_pairs(r1_fws1, r1_fws2, fastq1_qual, r2_fws1, r2_fws2, fastq2_qual)
            r1_rvs1, r1_rvs2, r1_rvscore = CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq1_seq), refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            r2_rvs1, r2_rvs2, r2_rvscore = CRISPResso2Align.global_align(CRISPRessoShared.reverse_complement(fastq2_seq), refs[ref_name]['sequence'], matrix=aln_matrix, gap_incentive=refs[ref_name]['gap_incentive'], gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
            rvs1, rvs2, rvscore = get_consensus_alignment_from_pairs(r1_rvs1, r1_rvs2, fastq1_qual, r2_rvs1, r2_rvs2, fastq2_qual)

            s1 = fws1
            s2 = fws2
            score = fwscore
            if (rvscore > fwscore):
                s1 = rvs1
                s2 = rvs2
                score = rvscore

#                print "for " + ref_name + " got fws1: " + str(fws1) + " and fws2: " + str(fws2) + " score: " +str(fwscore)
        aln_scores.append(score)
        ref_aln_details.append((ref_name, s1, s2, score))

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

    if best_match_score > 0:
        new_variant = {}
        new_variant['count'] = 1
        new_variant['aln_ref_names'] = best_match_names
        new_variant['aln_scores'] = aln_scores
        new_variant['ref_aln_details'] = ref_aln_details
        new_variant['best_match_score'] = best_match_score
        class_names = []

        for idx in range(len(best_match_names)):
            best_match_name = best_match_names[idx]

            if args.use_legacy_insertion_quantification:
                payload = CRISPRessoCOREResources.find_indels_substitutions_legacy(best_match_s1s[idx], best_match_s2s[idx], refs[best_match_name]['include_idxs'])
            else:
                payload = CRISPRessoCOREResources.find_indels_substitutions(best_match_s1s[idx], best_match_s2s[idx], refs[best_match_name]['include_idxs'])

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

            new_variant['variant_'+best_match_name] = payload

        new_variant['class_name'] = "&".join(class_names)

    else:
        new_variant = {}
        new_variant['count'] = 1
        new_variant['aln_scores'] = aln_scores
        new_variant['ref_aln_details'] = ref_aln_details
        new_variant['best_match_score'] = best_match_score
        return new_variant #return new variant with best match score of 0, but include the scores of insufficient alignments

    #handle ambiguous alignments
    if len(best_match_names) > 1 and not args.expand_ambiguous_alignments: #got 'Ambiguous' -- don't count toward total (e.g. indels at each position for the ref)
            new_variant['class_name'] = 'AMBIGUOUS'

    #search for prime editing scaffold
    #if it's there, assign this to be assigned to the reference '"Scaffold-incorporated"'
    if args.prime_editing_pegRNA_scaffold_seq and 'Prime-edited' in best_match_names: #any scaffold extensions must be closer to the prime-edited sequence
        pe_read_possible_scaffold_loc = new_variant['variant_Prime-edited']['ref_positions'].index(pe_scaffold_dna_info[0]-1) + 1
        if new_variant['variant_Prime-edited']['aln_seq'][pe_read_possible_scaffold_loc:(pe_read_possible_scaffold_loc+len(pe_scaffold_dna_info[1]))] == pe_scaffold_dna_info[1]:
#            print('comparingHERE ' + new_variant['variant_Prime-edited']['aln_seq'][pe_read_possible_scaffold_loc:(pe_read_possible_scaffold_loc+len(pe_scaffold_dna_info[1])+5)] + ' from ' + new_variant['variant_Prime-edited']['aln_seq'] + ' and ' + new_variant['variant_Prime-edited']['aln_ref'])
            new_variant['aln_ref_names'] = ["Scaffold-incorporated"]
            new_variant['class_name'] = "Scaffold-incorporated"
            old_payload = deepcopy(new_variant['variant_Prime-edited']) #keep prime-edited allele and alignment
            old_payload['ref_name'] = "Scaffold-incorporated"
            new_variant['variant_'+"Scaffold-incorporated"] = old_payload

    return new_variant

def process_fastq(fastq_filename, variantCache, ref_names, refs, args):
    """process_fastq processes each of the reads contained in a fastq file, given a cache of pre-computed variants
        fastqIn: name of fastq (e.g. output of FLASH)
            This file can be gzipped or plain text

        variantCache: dict with keys: sequence
            dict with keys:
                'count' : number of time sequence was observed
                'aln_ref_names' : names of reference it was aligned to
                'aln_scores' : score of alignment to each reference
                'aln_details' # details (seq1, seq2, score) of alignment to each other reference sequence
                'class_name' : string with class names it was aligned to
                'best_match_score' : score of best match (0 if no alignments matched above amplicon threshold)
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
                # 'sgRNA_plot_cut_points'#whether cut points should be plotted (not plotted for base editing, prime editing flap)
                # 'sgRNA_plot_idxs' #indices along reference sequence for which to plot the allele plot (allele frequency plot around sgRNA)
                # 'sgRNA_names' #names of sgRNAs (in case there are multiple matches for a single sgRNA)
                # 'sgRNA_mismatches' #indices along guide that are mismatched against a 'flexiguide_seq'
                # 'sgRNA_orig_sequences' #original sgRNA sequences as passed in as parameters (e.g. not including flexiguide changes or case changes)
                # 'contains_coding_seq'
                # 'exon_positions'
                # 'exon_intervals'
                # 'exon_len_mods': the modification to the original exon length (if we copied the exon positions from another reference, this reference could introduce an indel, resulting in a non-zero length modification)
                # 'splicing_positions'
                # 'include_idxs' # sorted numpy array
                # 'exclude_idxs'
                # 'plot_idxs' #sorted numpy array
                # 'plot_name' #unique plotting name
                # 'idx_cloned_from' #if this reference didn't contain a guide (or exon sequence), it was aligned to 'idx_cloned_from' reference, and cut_points, plot_cut_points, gap_incentive, sgRNA_intervals, inculde_idx, ane exon information were cloned from it (at the appropriate indices)
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

    aln_matrix_loc = os.path.join(_ROOT, args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    pe_scaffold_dna_info = (0, None) #scaffold start loc, scaffold seq to search
    if args.prime_editing_pegRNA_scaffold_seq != "":
        pe_scaffold_dna_info = get_pe_scaffold_search(refs['Prime-edited']['sequence'], args.prime_editing_pegRNA_extension_seq, args.prime_editing_pegRNA_scaffold_seq, args.prime_editing_pegRNA_scaffold_min_match_length)

    not_aln = {} #cache for reads that don't align

    if fastq_filename.endswith('.gz'):
        fastq_handle = gzip.open(fastq_filename, 'rt')
    else:
        fastq_handle=open(fastq_filename)

    while(fastq_handle.readline()):
        #read through fastq in sets of 4
        fastq_seq = fastq_handle.readline().strip()
        fastq_plus = fastq_handle.readline()
        fastq_qual = fastq_handle.readline()

        if (N_TOT_READS % 10000 == 0):
            info("Processing reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))

        N_TOT_READS+=1

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
            new_variant = get_new_variant_object(args, fastq_seq, refs, ref_names, aln_matrix, pe_scaffold_dna_info)
            if new_variant['best_match_score'] <= 0:
                N_COMPUTED_NOTALN+=1
                not_aln[fastq_seq] = 1
            else:
                N_COMPUTED_ALN+=1
                variantCache[fastq_seq] = new_variant

    fastq_handle.close()

    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))
    aln_stats = {"N_TOT_READS" : N_TOT_READS,
               "N_CACHED_ALN" : N_CACHED_ALN,
               "N_CACHED_NOTALN" : N_CACHED_NOTALN,
               "N_COMPUTED_ALN" : N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN" : N_COMPUTED_NOTALN
               }
    return(aln_stats)


def process_paired_fastq(fastq1_filename, fastq2_filename, variantCache, ref_names, refs, args):
    """Processes paired reads by aligning each read to the reference sequence
        This method avoids the use of flash to merge paired-end reads

        fastq1_filename: input fastq read 1 file
        fastq2_filename: input fastq read 2 file

        variantCache: dict with keys: sequence dict with keys (described in process_fastq)
        ref_names: list of reference names
        refs: dictionary of sequences name>ref object
        args: crispresso2 args
        """

    N_TOT_READS = 0
    N_CACHED_ALN = 0 # read was found in cache
    N_CACHED_NOTALN = 0 #read was found in 'not aligned' cache
    N_COMPUTED_ALN = 0 # not in cache, aligned to at least 1 sequence with min cutoff
    N_COMPUTED_NOTALN = 0 #not in cache, not aligned to any sequence with min cutoff

    aln_matrix_loc = os.path.join(_ROOT, args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    pe_scaffold_dna_info = (0, None) #scaffold start loc, scaffold seq to search
    if args.prime_editing_pegRNA_scaffold_seq != "":
        pe_scaffold_dna_info = get_pe_scaffold_search(refs['Prime-edited']['sequence'], args.prime_editing_pegRNA_extension_seq, args.prime_editing_pegRNA_scaffold_seq, args.prime_editing_pegRNA_scaffold_min_match_length)

    not_aln = {} #cache for reads that don't align

    if fastq1_filename.endswith('.gz'):
        fastq1_handle = gzip.open(fastq1_filename, 'rt')
    else:
        fastq1_handle=open(fastq1_filename)

    if fastq2_filename.endswith('.gz'):
        fastq2_handle = gzip.open(fastq2_filename, 'rt')
    else:
        fastq2_handle=open(fastq2_filename)

    while(fastq1_handle.readline()):
        #read through fastq in sets of 4
        fastq1_seq = fastq1_handle.readline().strip()
        fastq1_plus = fastq1_handle.readline()
        fastq1_qual = fastq1_handle.readline()

        fastq2_id = fastq2_handle.readline()
        fastq2_seq = fastq2_handle.readline().strip()
        fastq2_plus = fastq2_handle.readline()
        fastq2_qual = fastq2_handle.readline()

        if (N_TOT_READS % 10000 == 0):
            info("Processing reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))

        N_TOT_READS+=1

        #if the sequence has been seen and can't be aligned, skip it
        #cache the sequence of both r1 and r2 sequences as lookup_fastq_seq
        lookup_fastq_seq = fastq1_seq + "+" + fastq2_seq
        if (lookup_fastq_seq in not_aln):
            N_CACHED_NOTALN += 1
            continue
        #if the sequence is already associated with a variant in the variant cache, pull it out
        if (lookup_fastq_seq in variantCache):
            N_CACHED_ALN+=1
            variantCache[lookup_fastq_seq]['count'] += 1

        #otherwise, create a new variant object, and put it in the cache
        else:
            new_variant = get_new_variant_object_from_paired(args, fastq1_seq, fastq2_seq, fastq1_qual, fastq2_qual, refs, ref_names, aln_matrix, pe_scaffold_dna_info)
            if new_variant['best_match_score'] <= 0:
                N_COMPUTED_NOTALN+=1
                not_aln[lookup_fastq_seq] = 1
            else:
                N_COMPUTED_ALN+=1
                variantCache[lookup_fastq_seq] = new_variant

    fastq1_handle.close()
    fastq2_handle.close()

    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))
    aln_stats = {"N_TOT_READS" : N_TOT_READS,
               "N_CACHED_ALN" : N_CACHED_ALN,
               "N_CACHED_NOTALN" : N_CACHED_NOTALN,
               "N_COMPUTED_ALN" : N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN" : N_COMPUTED_NOTALN
               }
    return(aln_stats)

def process_bam(bam_filename, bam_chr_loc, output_bam, variantCache, ref_names, refs, args):
    """
    process_bam processes each of the reads contained in a bam file, given a cache of pre-computed variants
        bam_filename: name of bam (e.g. output of bowtie2)
        bam_chr_loc: positions in the input bam to read from
        output_bam: name out bam to write to (includes crispresso alignment info)

        variantCache: dict with keys: sequence dict with keys (described in process_fastq)
        ref_names: list of reference names
        refs: dictionary of sequences name>ref object
        args: crispresso2 args
    """

    N_TOT_READS = 0
    N_CACHED_ALN = 0 # read was found in cache
    N_CACHED_NOTALN = 0 #read was found in 'not aligned' cache
    N_COMPUTED_ALN = 0 # not in cache, aligned to at least 1 sequence with min cutoff
    N_COMPUTED_NOTALN = 0 #not in cache, not aligned to any sequence with min cutoff

    aln_matrix_loc = os.path.join(_ROOT, args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    pe_scaffold_dna_info = (0, None) #scaffold start loc, scaffold sequence
    if args.prime_editing_pegRNA_scaffold_seq != "":
        pe_scaffold_dna_info = get_pe_scaffold_search(refs['Prime-edited']['sequence'], args.prime_editing_pegRNA_extension_seq, args.prime_editing_pegRNA_scaffold_seq, args.prime_editing_pegRNA_scaffold_min_match_length)

    not_aln = {} #cache for reads that don't align

    output_sam = output_bam+".sam"
    with open(output_sam, "w") as sam_out:
        #first, write header to sam
        proc = sb.Popen(['samtools', 'view', '-H', bam_filename], stdout=sb.PIPE, stderr=sb.PIPE, encoding='utf-8')
        for line in proc.stdout:
            sam_out.write(line)

        crispresso_cmd_to_write = ' '.join(sys.argv)
        sam_out.write('@PG\tID:crispresso2\tPN:crispresso2\tVN:'+CRISPRessoShared.__version__+'\tCL:"'+crispresso_cmd_to_write+'"\n')

        if bam_chr_loc != "":
            proc = sb.Popen(['samtools', 'view', bam_filename, bam_chr_loc], stdout=sb.PIPE, encoding='utf-8')
        else:
            proc = sb.Popen(['samtools', 'view', bam_filename], stdout=sb.PIPE, encoding='utf-8')
        for sam_line in proc.stdout:
            sam_line_els = sam_line.rstrip().split("\t")
            fastq_seq = sam_line_els[9]

            if (N_TOT_READS % 10000 == 0):
                info("Processing reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))

            N_TOT_READS+=1

            #if the sequence has been seen and can't be aligned, skip it
            if (fastq_seq in not_aln):
                N_CACHED_NOTALN += 1
                sam_line_els.append(not_aln[fastq_seq]) #Crispresso2 alignment: NA
                sam_out.write("\t".join(sam_line_els)+"\n")
                continue
            #if the sequence is already associated with a variant in the variant cache, pull it out
            if (fastq_seq in variantCache):
                N_CACHED_ALN+=1
                variantCache[fastq_seq]['count'] += 1
                #sam_line_els[5] = variantCache[fastq_seq]['sam_cigar']
                sam_line_els.append(variantCache[fastq_seq]['crispresso2_annotation'])
                sam_out.write("\t".join(sam_line_els)+"\n")

            #otherwise, create a new variant object, and put it in the cache
            else:
                new_variant = get_new_variant_object(args, fastq_seq, refs, ref_names, aln_matrix, pe_scaffold_dna_info)
                if new_variant['best_match_score'] <= 0:
                    N_COMPUTED_NOTALN+=1
                    crispresso_sam_optional_fields = "c2:Z:ALN=NA" +\
                            " ALN_SCORES=" + ('&'.join([str(x) for x in new_variant['aln_scores']])) +\
                            " ALN_DETAILS=" + ('&'.join([','.join([str(y) for y in x]) for x in new_variant['ref_aln_details']]))
                    sam_line_els.append(crispresso_sam_optional_fields)
                    not_aln[fastq_seq] = crispresso_sam_optional_fields
                    sam_out.write("\t".join(sam_line_els)+"\n")
                else:
                    N_COMPUTED_ALN+=1

                    class_names = []
                    ins_inds = []
                    del_inds = []
                    sub_inds = []
                    edit_strings = []

    #                for idx, best_match_name in enumerate(best_match_names):
                    for idx, best_match_name in enumerate(new_variant['aln_ref_names']):
                        payload=new_variant['variant_'+best_match_name]

                        del_inds.append([str(x[0][0])+"("+str(x[1])+")" for x in zip(payload['deletion_coordinates'], payload['deletion_sizes'])])

                        ins_vals = []
                        for ins_coord,ins_size in zip(payload['insertion_coordinates'],payload['insertion_sizes']):
                            ins_start = payload['ref_positions'].index(ins_coord[0])
                            ins_vals.append(payload['aln_seq'][ins_start+1:ins_start+1+ins_size])
                        ins_inds.append([str(x[0][0])+"("+str(x[1])+"+"+x[2]+")" for x in zip(payload['insertion_coordinates'], payload['insertion_sizes'], ins_vals)])

                        sub_inds.append(payload['substitution_positions'])
                        edit_strings.append('D'+str(int(payload['deletion_n']))+';I'+str(int(payload['insertion_n']))+';S'+str(int(payload['substitution_n'])))


                    crispresso_sam_optional_fields = "c2:Z:ALN="+("&".join(new_variant['aln_ref_names'])) +\
                            " CLASS="+new_variant['class_name']+\
                            " MODS="+("&".join(edit_strings))+\
                            " DEL="+("&".join([';'.join(x) for x in del_inds])) +\
                            " INS="+("&".join([';'.join(x) for x in ins_inds])) +\
                            " SUB=" + ("&".join([';'.join([str(y) for y in x]) for x in sub_inds])) +\
                            " ALN_REF=" + ('&'.join([new_variant['variant_'+name]['aln_ref'] for name in new_variant['aln_ref_names']])) +\
                            " ALN_SEQ=" + ('&'.join([new_variant['variant_'+name]['aln_seq'] for name in new_variant['aln_ref_names']]))
                    sam_line_els.append(crispresso_sam_optional_fields)

                    #cigar strings are in reference to the given amplicon, not to the genomic sequence to which this read is aligned..
                    #first_variant = new_variant['variant_'+new_variant['aln_ref_names'][0]]
                    #sam_cigar = ''.join(CRISPRessoShared.unexplode_cigar(''.join([CRISPRessoShared.CIGAR_LOOKUP[x] for x in zip(first_variant['aln_seq'],first_variant['aln_ref'])])))
                    #sam_line_els[5] = sam_cigar
                    #new_variant['sam_cigar'] = sam_cigar
                    new_variant['crispresso2_annotation'] = crispresso_sam_optional_fields

                    sam_out.write("\t".join(sam_line_els)+"\n")

                    variantCache[fastq_seq] = new_variant

    output_sam = output_bam+".sam"
    cmd = 'samtools view -Sb '+output_sam + '>'+output_bam + ' && samtools index ' + output_bam
    bam_status=sb.call(cmd, shell=True)
    if bam_status:
            raise CRISPRessoShared.BadParameterException(
                'Bam creation failed. Command used: {0}'.format(cmd),
            )

    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))
    aln_stats = {"N_TOT_READS": N_TOT_READS,
               "N_CACHED_ALN": N_CACHED_ALN,
               "N_CACHED_NOTALN": N_CACHED_NOTALN,
               "N_COMPUTED_ALN": N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN": N_COMPUTED_NOTALN,
               }
    return(aln_stats)

def process_fastq_write_out(fastq_input, fastq_output, variantCache, ref_names, refs, args):
    """
    process_fastq_write_out processes each of the reads contained in a fastq input file, given a cache of pre-computed variants. All reads are read in, analyzed, and written to output with annotation
        fastq_input: input fastq
        fastq_output: fastq to write out
        variantCache: dict with keys: sequence dict with keys (described in process_fastq)

        ref_names: list of reference names
        refs: dictionary of sequences name>ref object
        args: crispresso2 args
    """

    N_TOT_READS = 0
    N_CACHED_ALN = 0 # read was found in cache
    N_CACHED_NOTALN = 0 #read was found in 'not aligned' cache
    N_COMPUTED_ALN = 0 # not in cache, aligned to at least 1 sequence with min cutoff
    N_COMPUTED_NOTALN = 0 #not in cache, not aligned to any sequence with min cutoff

    aln_matrix_loc = os.path.join(_ROOT, args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    pe_scaffold_dna_info = (0, None) #scaffold start loc, scaffold sequence
    if args.prime_editing_pegRNA_scaffold_seq != "":
        pe_scaffold_dna_info = get_pe_scaffold_search(refs['Prime-edited']['sequence'], args.prime_editing_pegRNA_extension_seq, args.prime_editing_pegRNA_scaffold_seq, args.prime_editing_pegRNA_scaffold_min_match_length)
    not_aln = {} #cache for reads that don't align
    not_aln[''] = "" #add empty sequence to the not_aln in case the fastq has an extra newline at the end

    if fastq_input.endswith('.gz'):
        fastq_input_handle=gzip.open(fastq_input, 'rt')
    else:
        fastq_input_handle=open(fastq_input)

    fastq_out_handle = gzip.open(fastq_output, 'wt')

    fastq_id = fastq_input_handle.readline()
    while(fastq_id):
        #read through fastq in sets of 4
        fastq_seq = fastq_input_handle.readline().strip()
        fastq_plus = fastq_input_handle.readline().strip()
        fastq_qual = fastq_input_handle.readline()

        if (N_TOT_READS % 10000 == 0):
            info("Processing reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))

        N_TOT_READS+=1

        #if the sequence has been seen and can't be aligned, skip it
        if fastq_seq in not_aln:
            N_CACHED_NOTALN += 1
            fastq_out_handle.write(fastq_id+fastq_seq+"\n"+fastq_plus+not_aln[fastq_seq]+"\n"+fastq_qual) #not_aln[fastq_seq] is alignment: NA
        elif fastq_seq in variantCache: #if the sequence is already associated with a variant in the variant cache, pull it out
            N_CACHED_ALN+=1
            variantCache[fastq_seq]['count'] += 1
            fastq_out_handle.write(fastq_id+fastq_seq+"\n"+fastq_plus+variantCache[fastq_seq]['crispresso2_annotation']+"\n"+fastq_qual)

        #otherwise, create a new variant object, and put it in the cache
        else:
            new_variant = get_new_variant_object(args, fastq_seq, refs, ref_names, aln_matrix, pe_scaffold_dna_info)
            if new_variant['best_match_score'] <= 0:
                N_COMPUTED_NOTALN+=1
                crispresso2_annotation = " ALN=NA" +\
                        " ALN_SCORES=" + ('&'.join([str(x) for x in new_variant['aln_scores']])) +\
                        " ALN_DETAILS=" + ('&'.join([','.join([str(y) for y in x]) for x in new_variant['ref_aln_details']]))
                not_aln[fastq_seq] = crispresso2_annotation
                fastq_out_handle.write(fastq_id+fastq_seq+"\n"+fastq_plus+crispresso2_annotation+"\n"+fastq_qual)
            else:
                N_COMPUTED_ALN+=1
                ins_inds = []
                del_inds = []
                sub_inds = []
                edit_strings = []

#                for idx, best_match_name in enumerate(best_match_names):
                for idx, best_match_name in enumerate(new_variant['aln_ref_names']):
                    payload=new_variant['variant_'+best_match_name]

                    del_inds.append([str(x[0][0])+"("+str(x[1])+")" for x in zip(payload['deletion_coordinates'], payload['deletion_sizes'])])

                    ins_vals = []
                    for ins_coord,ins_size in zip(payload['insertion_coordinates'],payload['insertion_sizes']):
                        ins_start = payload['ref_positions'].index(ins_coord[0])
                        ins_vals.append(payload['aln_seq'][ins_start+1:ins_start+1+ins_size])
                    ins_inds.append([str(x[0][0])+"("+str(x[1])+"+"+x[2]+")" for x in zip(payload['insertion_coordinates'], payload['insertion_sizes'], ins_vals)])

                    sub_inds.append(payload['substitution_positions'])
                    edit_strings.append('D'+str(int(payload['deletion_n']))+';I'+str(int(payload['insertion_n']))+';S'+str(int(payload['substitution_n'])))

                crispresso2_annotation = " ALN="+("&".join(new_variant['aln_ref_names'])) +\
                        " ALN_SCORES=" + ('&'.join([str(x) for x in new_variant['aln_scores']])) +\
                        " ALN_DETAILS=" + ('&'.join([','.join([str(y) for y in x]) for x in new_variant['ref_aln_details']])) +\
                        " CLASS="+new_variant['class_name']+\
                        " MODS="+("&".join(edit_strings))+\
                        " DEL="+("&".join([';'.join(x) for x in del_inds])) +\
                        " INS="+("&".join([';'.join(x) for x in ins_inds])) +\
                        " SUB=" + ("&".join([';'.join([str(y) for y in x]) for x in sub_inds])) +\
                        " ALN_REF=" + ('&'.join([new_variant['variant_'+name]['aln_ref'] for name in new_variant['aln_ref_names']])) +\
                        " ALN_SEQ=" + ('&'.join([new_variant['variant_'+name]['aln_seq'] for name in new_variant['aln_ref_names']]))

                new_variant['crispresso2_annotation'] = crispresso2_annotation

                fastq_out_handle.write(fastq_id+fastq_seq+"\n"+fastq_plus+crispresso2_annotation+"\n"+fastq_qual)

                variantCache[fastq_seq] = new_variant

        #last step of loop = read next line
        fastq_id = fastq_input_handle.readline()
    fastq_input_handle.close()
    fastq_out_handle.close()

    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))
    aln_stats = {"N_TOT_READS": N_TOT_READS,
               "N_CACHED_ALN": N_CACHED_ALN,
               "N_CACHED_NOTALN": N_CACHED_NOTALN,
               "N_COMPUTED_ALN": N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN": N_COMPUTED_NOTALN,
               }
    return(aln_stats)

def process_single_fastq_write_bam_out(fastq_input, bam_output, bam_header, variantCache, ref_names, refs, args):
    """
    process_fastq_write_out processes each of the reads contained in a fastq input file, given a cache of pre-computed variants. All reads are read in, analyzed, and written to output with annotation

    **Note that in this mode, refs must contain the following values for each reference:
        aln_chr: location where reads will be aligned - this is the reference sequence
        aln_start: alignment start position of reference sequence
        aln_end: alignment end position of reference sequence
        aln_strand: strand reference aligned to

    **Note that we expect one input fastq file (R1 and R2 merged for paired sequencing). At some point, we may want to accept paired ends and record each read separately.

        fastq_input: input fastq
        bam_output: bam to write out
        bam_header: bam header to write to output bam file (eg including crispresso version and command used)
        variantCache: dict with keys: sequence dict with keys (described in process_fastq)

        ref_names: list of reference names
        refs: dictionary of sequences name>ref object
        args: crispresso2 args
    """

    N_TOT_READS = 0
    N_CACHED_ALN = 0  # read was found in cache
    N_CACHED_NOTALN = 0  # read was found in 'not aligned' cache
    N_COMPUTED_ALN = 0  # not in cache, aligned to at least 1 sequence with min cutoff
    N_COMPUTED_NOTALN = 0  # not in cache, not aligned to any sequence with min cutoff

    aln_matrix_loc = os.path.join(_ROOT, args.needleman_wunsch_aln_matrix_loc)
    CRISPRessoShared.check_file(aln_matrix_loc)
    aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

    pe_scaffold_dna_info = (0, None)  # scaffold start loc, scaffold sequence
    if args.prime_editing_pegRNA_scaffold_seq != "":
        pe_scaffold_dna_info = get_pe_scaffold_search(refs['Prime-edited']['sequence'], args.prime_editing_pegRNA_extension_seq, args.prime_editing_pegRNA_scaffold_seq, args.prime_editing_pegRNA_scaffold_min_match_length)
    not_aln = {}  # cache for reads that don't align
    not_aln[''] = ""  # add empty sequence to the not_aln in case the fastq has an extra newline at the end

    if fastq_input.endswith('.gz'):
        fastq_input_handle = gzip.open(fastq_input, 'rt')
    else:
        fastq_input_handle = open(fastq_input)

    sam_out = bam_output+".sam"
    sam_out_handle = open(sam_out, 'wt')

    # write sam output header
    sam_out_handle.write(bam_header)

    fastq_id = fastq_input_handle.readline().strip()[1:]
    while(fastq_id):
        # read through fastq in sets of 4
        fastq_seq = fastq_input_handle.readline().strip()
        fastq_plus = fastq_input_handle.readline().strip()
        fastq_qual = fastq_input_handle.readline().strip()

        if (N_TOT_READS % 10000 == 0):
            info("Processing reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))

        N_TOT_READS += 1

        # if the sequence has been seen and can't be aligned, skip it
        if fastq_seq in not_aln:
            N_CACHED_NOTALN += 1
            new_sam_entry = not_aln[fastq_seq][:]
            new_sam_entry[0] = fastq_id
            new_sam_entry[10] = fastq_qual
            sam_out_handle.write("\t".join(new_sam_entry)+"\n")  # not_aln[fastq_seq] is alignment: NA
        elif fastq_seq in variantCache:  # if the sequence is already associated with a variant in the variant cache, pull it out
            N_CACHED_ALN += 1
            variantCache[fastq_seq]['count'] += 1
            new_sam_entry = variantCache[fastq_seq]['sam_entry'][:]
            new_sam_entry[0] = fastq_id
            new_sam_entry[10] = fastq_qual
            sam_out_handle.write("\t".join(new_sam_entry)+"\n")  # write cached alignment with modified read id and qual

        # otherwise, create a new variant object, and put it in the cache
        else:
            new_variant = get_new_variant_object(args, fastq_seq, refs, ref_names, aln_matrix, pe_scaffold_dna_info)
            if new_variant['best_match_score'] <= 0:
                N_COMPUTED_NOTALN += 1
                new_sam_entry = [
                    fastq_id,  # read id
                    '4',             # flag = unmapped 0x4
                    '*',             # aln chr
                    '0',             # aln loc
                    '0',             # quality
                    '0',             # cigar
                    '*',             # next
                    '0',             # next
                    '0',             # tlen
                    fastq_seq,       # seq
                    fastq_qual       # qual
                ]
                sam_out_handle.write("\t".join(new_sam_entry)+"\n")  # write cached alignment with modified read id and qual
                crispresso_sam_optional_fields = "c2:Z:ALN=NA" +\
                        " ALN_SCORES=" + ('&'.join([str(x) for x in new_variant['aln_scores']])) +\
                        " ALN_DETAILS=" + ('&'.join([','.join([str(y) for y in x]) for x in new_variant['ref_aln_details']]))
                new_sam_entry.append(crispresso_sam_optional_fields)
                not_aln[fastq_seq] = new_sam_entry
                sam_out_handle.write("\t".join(new_sam_entry)+"\n")  # write cached alignment with modified read id and qual
            else:
                N_COMPUTED_ALN += 1
                ins_inds = []
                del_inds = []
                sub_inds = []
                edit_strings = []

#                for idx, best_match_name in enumerate(best_match_names):
                for idx, best_match_name in enumerate(new_variant['aln_ref_names']):
                    payload = new_variant['variant_'+best_match_name]

                    del_inds.append([str(x[0][0])+"("+str(x[1])+")" for x in zip(payload['deletion_coordinates'], payload['deletion_sizes'])])

                    ins_vals = []
                    for ins_coord, ins_size in zip(payload['insertion_coordinates'], payload['insertion_sizes']):
                        ins_start = payload['ref_positions'].index(ins_coord[0])
                        ins_vals.append(payload['aln_seq'][ins_start+1:ins_start+1+ins_size])
                    ins_inds.append([str(x[0][0])+"("+str(x[1])+"+"+x[2]+")" for x in zip(payload['insertion_coordinates'], payload['insertion_sizes'], ins_vals)])

                    sub_inds.append(payload['substitution_positions'])
                    edit_strings.append('D'+str(int(payload['deletion_n']))+';I'+str(int(payload['insertion_n']))+';S'+str(int(payload['substitution_n'])))

                first_ref_name = new_variant['aln_ref_names'][0]
                first_payload = new_variant['variant_' + first_ref_name]
                exploded_cigar = ''.join([CRISPRessoShared.CIGAR_LOOKUP[x] for x in zip(first_payload['aln_seq'], first_payload['aln_ref'])])
                sam_cigar_els = CRISPRessoShared.unexplode_cigar(exploded_cigar)
                sam_cigar = ''.join(sam_cigar_els)

                sam_flag = '0'
                if first_payload['aln_strand'] == '-':
                    sam_flag = '16'

                this_aln_loc = refs[first_ref_name]['aln_start']
#                #if cigar starts with deletions
#                if sam_cigar_els[0].endswith('D'):
#                    num = int(sam_cigar_els[0][0:-1])
#                    this_aln_loc += num
#                elif sam_cigar_els[0].endswith('I'):
#                    num = int(sam_cigar_els[0][0:-1])
#                    this_aln_loc -= num
                seq_to_write = fastq_seq
                qual_to_write = fastq_qual
                cigar_to_write = sam_cigar
                if refs[first_ref_name]['aln_strand'] == '-':
                    if sam_flag == '16':
                        sam_flag = '0'
                    else:
                        sam_flag = '16'
                    cigar_to_write = ''.join(sam_cigar_els[::-1])
                    seq_to_write = CRISPRessoShared.reverse_complement(fastq_seq)
                    qual_to_write = fastq_qual[::-1]

                new_sam_entry = [
                    fastq_id,                                   # read id
                    sam_flag,                                   # flag = mapped
                    refs[first_ref_name]['aln_chr'],            # aln chr
                    str(this_aln_loc),                          # aln loc
                    str(int(new_variant['best_match_score'])),  # quality
                    cigar_to_write,                             # cigar
                    '*',                                        # next
                    '0',                                        # next
                    '0',                                        # tlen
                    seq_to_write,                               # seq
                    qual_to_write                               # qual
                ]

                crispresso_sam_optional_fields = "c2:Z:ALN="+("&".join(new_variant['aln_ref_names'])) +\
                        " ALN_SCORES=" + ('&'.join([str(x) for x in new_variant['aln_scores']])) +\
                        " ALN_DETAILS=" + ('&'.join([','.join([str(y) for y in x]) for x in new_variant['ref_aln_details']])) +\
                        " CLASS="+new_variant['class_name']+\
                        " MODS="+("&".join(edit_strings))+\
                        " DEL="+("&".join([';'.join(x) for x in del_inds])) +\
                        " INS="+("&".join([';'.join(x) for x in ins_inds])) +\
                        " SUB=" + ("&".join([';'.join([str(y) for y in x]) for x in sub_inds])) +\
                        " ALN_REF=" + ('&'.join([new_variant['variant_'+name]['aln_ref'] for name in new_variant['aln_ref_names']])) +\
                        " ALN_SEQ=" + ('&'.join([new_variant['variant_'+name]['aln_seq'] for name in new_variant['aln_ref_names']]))

                new_sam_entry.append(crispresso_sam_optional_fields)
                new_variant['sam_entry'] = new_sam_entry

                sam_out_handle.write("\t".join(new_sam_entry)+"\n")  # write cached alignment with modified read id and qual

                variantCache[fastq_seq] = new_variant

        #last step of loop = read next line
        fastq_id = fastq_input_handle.readline().strip()[1:]
    fastq_input_handle.close()
    sam_out_handle.close()

    sort_and_index_cmd = 'samtools sort ' + sam_out + ' -o ' + bam_output + ' && samtools index ' + bam_output
    sort_bam_status = sb.call(sort_and_index_cmd, shell=True)
    if sort_bam_status:
        raise CRISPRessoShared.BadParameterException(
            'Bam sort failed. Command used: {0}'.format(sort_and_index_cmd)
        )

    if not args.debug:
        os.remove(sam_out)

    info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))
    aln_stats = {"N_TOT_READS": N_TOT_READS,
               "N_CACHED_ALN": N_CACHED_ALN,
               "N_CACHED_NOTALN": N_CACHED_NOTALN,
               "N_COMPUTED_ALN": N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN": N_COMPUTED_NOTALN,
               }
    return(aln_stats)


def split_interleaved_fastq(fastq_filename, output_filename_r1, output_filename_r2):
    if fastq_filename.endswith('.gz'):
        fastq_handle = gzip.open(fastq_filename, 'rt')
    else:
        fastq_handle=open(fastq_filename)

    #we cannot use with on gzip with python 2.6 :(
    try:
        fastq_splitted_outfile_r1 = gzip.open(output_filename_r1, 'wt')
        fastq_splitted_outfile_r2 = gzip.open(output_filename_r2, 'wt')
        [fastq_splitted_outfile_r1.write(line) if (i % 8 < 4) else fastq_splitted_outfile_r2.write(line) for i, line in enumerate(fastq_handle)]
    except:
        raise CRISPRessoShared.BadParameterException('Error in splitting read pairs from a single file')

    return output_filename_r1, output_filename_r2


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
        description = ['~~~CRISPResso 2~~~', '-Analysis of genome editing outcomes from deep sequencing data-']
        header = CRISPRessoShared.get_crispresso_header(description=description, header_str=None)
        print(header)

        arg_parser = CRISPRessoShared.getCRISPRessoArgParser()
        args = arg_parser.parse_args()

        aln_matrix_loc = os.path.join(_ROOT, "EDNAFULL")
        CRISPRessoShared.check_file(aln_matrix_loc)
        aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)

        n_processes = 1
        if args.n_processes == "max":
            n_processes = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes = int(args.n_processes)

        #check files and get output name
        if args.fastq_r1:
            CRISPRessoShared.check_file(args.fastq_r1)
            if args.fastq_r2:
                CRISPRessoShared.check_file(args.fastq_r2)
        elif args.bam_input:
            CRISPRessoShared.check_file(args.bam_input)
        else:
            arg_parser.print_help()
            raise CRISPRessoShared.BadParameterException('Please provide input data for analysis e.g. using the --fastq_r1 parameter.')

        if args.assign_ambiguous_alignments_to_first_reference and args.expand_ambiguous_alignments:
            arg_parser.print_help()
            raise CRISPRessoShared.BadParameterException('The options --assign_ambiguous_alignments_to_first_reference and --expand_ambiguous_alignments are mutually exclusive and cannot be both set.')

        if args.amplicon_seq is None and args.auto is False:
            arg_parser.print_help()
            raise CRISPRessoShared.BadParameterException('Please provide an amplicon sequence for analysis using the --amplicon_seq parameter.')

        if (args.needleman_wunsch_gap_open > 0):
            raise CRISPRessoShared.BadParameterException("Needleman Wunsch gap open penalty must be <= 0")
        if (args.needleman_wunsch_gap_extend > 0):
            raise CRISPRessoShared.BadParameterException("Needleman Wunsch gap extend penalty must be <= 0")


        #create output directory
        get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq', '').replace('.gz', '').replace('.fq', '')
        get_name_from_bam=lambda  x: os.path.basename(x).replace('.bam', '')

        #normalize name and remove not allowed characters
        if not args.name:
            if args.fastq_r2!='':
                database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1), get_name_from_fasta(args.fastq_r2))
            elif args.fastq_r1 != '':
                database_id='%s' % get_name_from_fasta(args.fastq_r1)
            elif args.bam_input != '':
                database_id='%s' % get_name_from_bam(args.bam_input)
        else:
            clean_name=CRISPRessoShared.slugify(args.name)
            if args.name!= clean_name:
                warn('The specified name %s contained invalid characters and was changed to: %s' % (args.name, clean_name))
            database_id=clean_name

        clean_file_prefix = ""
        if args.file_prefix != "":
            clean_file_prefix = CRISPRessoShared.slugify(args.file_prefix)
            if not clean_file_prefix.endswith("."):
                clean_file_prefix += "."

        OUTPUT_DIRECTORY='CRISPResso_on_%s' % database_id

        if args.output_folder:
            OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY, clean_file_prefix + filename) #handy function to put a file in the output directory

        crispresso2_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        log_filename=_jp('CRISPResso_RUNNING_LOG.txt')
        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        crispresso2_info['running_info']['name'] = database_id

        crispresso_cmd_to_write = ' '.join(sys.argv)
        if args.write_cleaned_report:
            cmd_copy = sys.argv[:]
            cmd_copy[0] = 'CRISPResso'
            for i in range(len(cmd_copy)):
                if os.sep in cmd_copy[i]:
                    cmd_copy[i] = os.path.basename(cmd_copy[i])

            crispresso_cmd_to_write = ' '.join(cmd_copy) #clean command doesn't show the absolute path to the executable or other files
        crispresso2_info['running_info']['command_used'] = crispresso_cmd_to_write

        try:
            os.makedirs(OUTPUT_DIRECTORY)
            info('Creating Folder %s' % OUTPUT_DIRECTORY)
#            info('Done!') #crispresso2 doesn't announce that the folder is created... save some electricity here
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        finally:
            logger.addHandler(logging.FileHandler(log_filename))

            with open(log_filename, 'w+') as outfile:
                outfile.write('CRISPResso version %s\n[Command used]:\n%s\n\n[Execution log]:\n' %(CRISPRessoShared.__version__, crispresso_cmd_to_write))

        files_to_remove = [] #these files will be deleted at the end of the run

        if args.no_rerun:
            if os.path.exists(crispresso2_info_file):
                previous_run_data = CRISPRessoShared.load_crispresso_info(OUTPUT_DIRECTORY)
                if previous_run_data['running_info']['version'] == CRISPRessoShared.__version__:
                    args_are_same = True
                    for arg in vars(args):
                        if arg == "no_rerun":
                            continue
                        if arg not in vars(previous_run_data['running_info']['args']):
                            info('Comparing current run to previous run: old run had argument ' + str(arg) + ' \nRerunning.')
                            args_are_same = False
                        elif str(getattr(previous_run_data['running_info']['args'], arg)) != str(getattr(args, arg)):
                            info('Comparing current run to previous run:\n\told argument ' + str(arg) + ' = ' + str(getattr(previous_run_data['running_info']['args'], arg)) + '\n\tnew argument: ' + str(arg) + ' = ' + str(getattr(args, arg)) + '\nRerunning.')
                            args_are_same = False

                    if args_are_same:
                        info('Analysis already completed on %s!'%previous_run_data['running_info']['end_time_string'])
                        sys.exit(0)
                else:
                    info('The no_rerun flag is set, but this analysis will be rerun because the existing run was performed using an old version of CRISPResso (' + str(previous_run_data['running_info']['version']) + ').')


        def rreplace(s, old, new):
            li = s.rsplit(old)
            return new.join(li)
        # if bam input, make sure bam is sorted and indexed
        if args.bam_input:
            if args.bam_chr_loc != "":  # we only need an index if we're accessing specific positions later
                if os.path.exists(rreplace(args.bam_input, ".bam", ".bai")):
                    info('Index file for input .bam file exists, skipping generation.')
                elif os.path.exists(args.bam_input+'.bai'):
                    info('Index file for input .bam file exists, skipping generation.')
                else:
                    info('Creating index file for input .bam file...')
                    bam_input_file = _jp(os.path.basename(args.bam_input))+".sorted.bam"
                    sb.call('samtools sort -o '+bam_input_file+' ' + args.bam_input, shell=True)
                    sb.call('samtools index %s ' % (bam_input_file), shell=True)
                    files_to_remove.append(bam_input_file)
                    files_to_remove.append(bam_input_file+".bai")
                    args.bam_input = bam_input_file
            crispresso2_info['running_info']['bam_input'] = args.bam_input
            bam_output = _jp('CRISPResso_output.bam')
            crispresso2_info['running_info']['bam_output'] = bam_output
            output_sam = bam_output+".sam"
            files_to_remove.append(output_sam)

            if args.fastq_output:
                raise CRISPRessoShared.BadParameterException('bam_input is not compatable with fastq_output! Please either use bam_input or fastq_output.')

        if args.fastq_output:
            fastq_output = _jp('CRISPResso_output.fastq.gz')
            crispresso2_info['fastq_output'] = fastq_output
            info('Writing fastq output file: ' + fastq_output)
        if args.bam_output:
            if args.fastq_output:
                raise CRISPRessoShared.BadParameterException('bam_output is not compatable with fastq_output! Please either use bam_output or fastq_output.')
            bam_output = _jp('CRISPResso_output.bam')
            crispresso2_info['bam_output'] = bam_output
            info('Writing bam output file: ' + bam_output)


        #### ASSERT GUIDE(S)
        guides = []
        if args.guide_seq:
            for current_guide_seq in args.guide_seq.split(','):
                wrong_nt=CRISPRessoShared.find_wrong_nt(current_guide_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The sgRNA sequence contains bad characters:%s'  % ' '.join(wrong_nt))
                guides.append(current_guide_seq)

        #each guide has a name, quantification window centers (relative to the end of the guide), and quantification window sizes
        guide_names = [""]*len(guides)
        if args.guide_name:
            guide_name_arr = args.guide_name.split(",")
            if len(guide_name_arr) > len(guides):
                raise CRISPRessoShared.BadParameterException("More guide names were given than guides. Guides: %d Guide names: %d"%(len(guides), len(guide_name_arr)))
            for idx, guide_name in enumerate(guide_name_arr):
                if guide_name != "":
                    guide_names[idx] = guide_name

        guide_qw_centers = CRISPRessoShared.set_guide_array(args.quantification_window_center, guides, 'guide quantification center')
        guide_qw_sizes = CRISPRessoShared.set_guide_array(args.quantification_window_size, guides, 'guide quantification size')

        guide_plot_cut_points = [True]*len(guides) #whether to plot cut point -- prime editing flap cut points aren't plotted
        if (args.base_editor_output):
            guide_plot_cut_points = [False]*len(guides) #whether to plot cut point -- base editor and prime editing flap cut points aren't plotted


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
            number_of_reads_to_consider = 5000

            auto_fastq_r1 = args.fastq_r1 #paths to fastq files for performing auto functions
            auto_fastq_r2 = args.fastq_r2
            if args.bam_input != "": #if input is a bam, create temp files with reads for processing here
                p = sb.Popen("samtools view -c -f 1 " + args.bam_input + " " + args.bam_chr_loc, shell=True, stdout=sb.PIPE)
                n_reads_paired = int(float(p.communicate()[0]))
                if n_reads_paired > 0:
                    auto_fastq_r1 = _jp(os.path.basename(args.bam_input)+".r1.fastq")
                    auto_fastq_r2 = _jp(os.path.basename(args.bam_input)+".r2.fastq")
                    files_to_remove.append(auto_fastq_r1)
                    files_to_remove.append(auto_fastq_r2)
                    sb.call('samtools view -h ' + args.bam_input + ' ' + args.bam_chr_loc + ' | head -n ' + str(number_of_reads_to_consider+100) + ' | ' + \
                        'samtools bam2fq -1 ' + auto_fastq_r1 + ' -2 ' + auto_fastq_r2 + ' -0 /dev/null -s /dev/null -n 2> /dev/null > ', shell=True) #-0 READ_OTHER, -s singleton, -n don't append /1 and /2 to read name
                else:
                    auto_fastq_r1 = _jp(os.path.basename(args.bam_input)+".r1.fastq")
                    auto_fastq_r2 = ""
                    files_to_remove.append(auto_fastq_r1)
                    sb.call('samtools view -h ' + args.bam_input + ' ' + args.bam_chr_loc + ' | head -n ' + str(number_of_reads_to_consider+100) + ' | samtools bam2fq -  2> /dev/null > ' + auto_fastq_r1, shell=True)

            amplicon_seq_arr = CRISPRessoShared.guess_amplicons(auto_fastq_r1, auto_fastq_r2, number_of_reads_to_consider, args.flash_command, args.max_paired_end_reads_overlap, args.min_paired_end_reads_overlap,
                aln_matrix, args.needleman_wunsch_gap_open, args.needleman_wunsch_gap_extend)
            amp_dummy = ['']
            amp_dummy.extend(list(range(2, len(amplicon_seq_arr)+1)))
            amplicon_name_arr = ['Inferred'+str(x) for x in amp_dummy]

            if args.amplicon_seq is not None: #if user provides an amplicon sequence, add it to the guessed amplicons
                amp_names = args.amplicon_name.split(",")
                for idx, amp_seq in enumerate(args.amplicon_seq.split(",")):
                    if amp_seq not in amplicon_seq_arr:
                        amplicon_seq_arr.append(amp_seq)
                        if idx < len(amp_names):
                            amplicon_name_arr.append(amp_names[idx])
                        else:
                            amplicon_name_arr.append('User-defined amplicon' + idx)

            if len(guides) == 0:
                for amplicon_seq in amplicon_seq_arr:
                    (potential_guide, is_base_editor) = CRISPRessoShared.guess_guides(amplicon_seq, auto_fastq_r1, auto_fastq_r2, number_of_reads_to_consider, args.flash_command,
                        args.max_paired_end_reads_overlap, args.min_paired_end_reads_overlap, aln_matrix, args.needleman_wunsch_gap_open, args.needleman_wunsch_gap_extend)
                    if potential_guide is not None and potential_guide not in guides:
                        guides.append(potential_guide)
                        guide_names.append('Guessed sgRNA')
                        guide_qw_centers.append(int(args.quantification_window_center.split(",")[0]))
                        guide_qw_sizes.append(int(args.quantification_window_size.split(",")[0]))
                        guide_plot_cut_points.append(True)

            amplicon_min_alignment_score_arr = []
            plural_string = ""
            if len(amplicon_seq_arr) > 1:
                plural_string = "s"
            info("Auto-detected %d reference amplicon%s"%(len(amplicon_seq_arr), plural_string))

            if args.debug:
                for idx, seq in enumerate(amplicon_seq_arr):
                    info('Detected amplicon ' + str(idx) + ":" + str(seq))

            if len(guides) > 1:
                plural_string = "s"
            info("Auto-detected %d guide%s"%(len(guides), plural_string))
            if args.debug:
                for idx, seq in enumerate(guides):
                    info('Detected guide ' + str(idx) + ":" + str(seq))
            if amplicon_seq_arr == 0:
                raise CRISPRessoShared.BadParameterException(
                    "Cannot automatically infer amplicon sequence.",
                )

        else: #not auto
            amplicon_seq_arr = args.amplicon_seq.split(",")
            amplicon_name_arr = args.amplicon_name.split(",")
            #split on commas, only accept empty values
            amplicon_min_alignment_score_arr = [float(x) for x in args.amplicon_min_alignment_score.split(",") if x]

            for idx, amp_seq in enumerate(amplicon_seq_arr):
                if len(amp_seq) == 0:
                    raise CRISPRessoShared.BadParameterException(
                        "Amplicon {0} has length 0. Please check the --amplicon_seq parameter.".format(
                            idx + 1,
                        ),
                    )

                this_name = 'Amplicon'+str(idx)
                if idx >= len(amplicon_name_arr):
                    amplicon_name_arr.append(this_name)

                this_min_aln_score = args.default_min_aln_score
                if idx >= len(amplicon_min_alignment_score_arr):
                    amplicon_min_alignment_score_arr.append(this_min_aln_score)

        if args.expected_hdr_amplicon_seq != "":
            amplicon_seq_arr.append(args.expected_hdr_amplicon_seq)
            amplicon_name_arr.append('HDR')
            amplicon_min_alignment_score_arr.append(args.default_min_aln_score)

        amplicon_quant_window_coordinates_arr = ['']*len(amplicon_seq_arr)
        if args.quantification_window_coordinates is not None:
            #fill in quant window coordinates if they are given for each amplicon, because we're about to add some amplicons
            for idx, coords in enumerate(args.quantification_window_coordinates.strip("'").strip('"').split(",")):
                if coords != "":
                    amplicon_quant_window_coordinates_arr[idx] = coords


        #Prime editing
        prime_editing_extension_seq_dna = "" #global var for the editing extension sequence for the scaffold quantification below
        prime_editing_edited_amp_seq = ""
        if args.prime_editing_pegRNA_extension_seq != "":
#            if args.amplicon_name_arr eq "": #considering changing the name of the default amplicon from 'Reference' to 'Unedited'
#                amplicon_name_arr[0] = 'Unedited'
            extension_seq_dna_top_strand = args.prime_editing_pegRNA_extension_seq.upper().replace('U', 'T')
            wrong_nt=CRISPRessoShared.find_wrong_nt(extension_seq_dna_top_strand)
            if wrong_nt:
                raise CRISPRessoShared.NTException('The prime editing pegRNA extension sequence contains bad characters:%s'  % ' '.join(wrong_nt))
            prime_editing_extension_seq_dna = CRISPRessoShared.reverse_complement(extension_seq_dna_top_strand)

            #check to make sure the pegRNA spacer seq is in the RTT/extension seq
            #this is critical because we need to know where to check for scaffold incorporation. The pegRNA spacer should be FW and the RTT/extension should RC compared to the reference
            #If down the road we want users to be able to give this flexibly, we need to update the scaffold incorporation search in get_new_variant_object as well as in the definition of the quant window
            if args.prime_editing_pegRNA_spacer_seq == "":
                raise CRISPRessoShared.BadParameterException('The prime editing pegRNA spacer sequence (--prime_editing_pegRNA_spacer_seq) is required for prime editing analysis.')
            pegRNA_spacer_seq = args.prime_editing_pegRNA_spacer_seq.upper().replace('U', 'T')

            #check that the pegRNA aligns to the reference (and not the RC)
            amp_incentive = np.zeros(len(amplicon_seq_arr[0])+1,dtype=int)
            f1,f2,fw_score=CRISPResso2Align.global_align(pegRNA_spacer_seq,amplicon_seq_arr[0].upper(),matrix=aln_matrix,gap_incentive=amp_incentive,gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
            r1,r2,rv_score=CRISPResso2Align.global_align(pegRNA_spacer_seq,CRISPRessoShared.reverse_complement(amplicon_seq_arr[0].upper()),matrix=aln_matrix,gap_incentive=amp_incentive,gap_open=args.needleman_wunsch_gap_open,gap_extend=args.needleman_wunsch_gap_extend,)
            if rv_score > fw_score:
                raise CRISPRessoShared.BadParameterException('The prime editing pegRNA spacer sequence appears to be given in the 3\'->5\' order. The prime editing pegRNA spacer sequence (--prime_editing_pegRNA_spacer_seq) must be given in the RNA 5\'->3\' order.')

            ref_incentive = np.zeros(len(prime_editing_extension_seq_dna)+1, dtype=int)
            f1, f2, fw_score=CRISPResso2Align.global_align(pegRNA_spacer_seq, prime_editing_extension_seq_dna, matrix=aln_matrix, gap_incentive=ref_incentive, gap_open=args.needleman_wunsch_gap_open, gap_extend=0,)
            r1, r2, rv_score=CRISPResso2Align.global_align(pegRNA_spacer_seq, extension_seq_dna_top_strand, matrix=aln_matrix, gap_incentive=ref_incentive, gap_open=args.needleman_wunsch_gap_open, gap_extend=0,)
            if rv_score > fw_score:
                raise CRISPRessoShared.BadParameterException("The pegRNA spacer aligns to the pegRNA extension sequence in 3'->5' direction. The prime editing pegRNA spacer sequence (--prime_editing_pegRNA_spacer_seq) must be given in the RNA 5'->3' order, and the pegRNA extension sequence (--prime_editing_pegRNA_extension_seq) must be given in the 5'->3' order. In other words, the pegRNA spacer sequence should be found in the given reference sequence, and the reverse complement of the pegRNA extension sequence should be found in the reference sequence.")

            #setting refs['Prime-edited']['sequence']
            #first, align the extension seq to the reference amplicon
            #we're going to consider the first reference only (so if multiple alleles exist at the editing position, this may get messy)
            best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(prime_editing_extension_seq_dna, amplicon_seq_arr[0],aln_matrix,args.needleman_wunsch_gap_open,0)
            new_ref = s2[0:best_aln_start] + prime_editing_extension_seq_dna + s2[best_aln_end:]
            if args.debug:
                info('Alignment between extension sequence and reference sequence: \n' + s1 + '\n' + s2)

            if args.prime_editing_override_prime_edited_ref_seq != "":
                if args.debug:
                    info('Using provided prime_editing_override_prime_edited_ref_seq\nDetected: ' + new_ref + '\nProvided: ' + args.prime_editing_override_prime_edited_ref_seq)
                new_ref = args.prime_editing_override_prime_edited_ref_seq
                if args.prime_editing_pegRNA_extension_seq not in new_ref and CRISPRessoShared.reverse_complement(args.prime_editing_pegRNA_extension_seq) not in new_ref:
                    raise CRISPRessoShared.BadParameterException('The provided extension sequence is not in the prime edited override reference sequence!')

            if new_ref in amplicon_seq_arr:
                raise CRISPRessoShared.BadParameterException('The calculated prime-edited amplicon is the same as the reference sequence.')
            amplicon_seq_arr.append(new_ref)
            if 'Prime-edited' in amplicon_name_arr:
                raise CRISPRessoShared.BadParameterException("An amplicon named 'Prime-edited' must not be provided.")
            amplicon_name_arr.append('Prime-edited')
            amplicon_quant_window_coordinates_arr.append('')
            prime_editing_edited_amp_seq = new_ref

        #finished prime editing bit

        for idx, this_seq in enumerate(amplicon_seq_arr):
            wrong_nt=CRISPRessoShared.find_wrong_nt(this_seq)
            if wrong_nt:
                this_name = amplicon_name_arr[idx]
                raise CRISPRessoShared.NTException('Reference amplicon sequence %d (%s) contains invalid characters: %s'%(idx, this_name, ' '.join(wrong_nt)))



        def get_prime_editing_guides(this_amp_seq, this_amp_name, ref0_seq, prime_edited_seq, prime_editing_extension_seq_dna, prime_editing_pegRNA_extension_quantification_window_size, nicking_qw_center, nicking_qw_size,aln_matrix,needleman_wunsch_gap_open,needleman_wunsch_gap_extend):
            """
            gets prime editing guide sequences for this amplicon
            this_amp_seq : sequence of this amplicon
            this_amp_name : name of this amplicon
            ref0_seq : sequence of the 0th amplicon (the wildtype amplicon)
            prime_editing_edited_amp_seq : sequence of the edited amplicon
            prime_editing_extension_seq_dna : the extension sequence on the DNA strand
            prime_editing_pegRNA_extension_quantification_window_size : qw size for extension sequence (starts at end of extension sequence)
            nicking_qw_center: for the nicking guides, what is the quantification center (usually int(args.quantification_window_center.split(",")[0]))
            nicking_qw_size: for the nicking guides, what is the quantification size (usually int(args.quantification_window_size.split(",")[0]))
            aln_matrix: matrix specifying alignment substitution scores in the NCBI format
            needleman_wunsch_gap_open: alignment penalty assignment used to determine similarity of two sequences
            needleman_wunsch_gap_extend: alignment penalty assignment used to determine similarity of two sequences
            """
            pe_guides = []
            pe_orig_guide_seqs = []
            pe_guide_mismatches = []
            pe_guide_names = []
            pe_guide_qw_centers = []
            pe_guide_qw_sizes = []
            pe_guide_plot_cut_points = []

            #editing extension aligns to the prime-edited sequence only
            #if this is the prime edited sequence, add it directly
            if this_amp_seq == prime_editing_edited_amp_seq:
                best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(prime_editing_extension_seq_dna, prime_editing_edited_amp_seq,aln_matrix,needleman_wunsch_gap_open,0)
                pe_guides.append(best_aln_seq)
                pe_orig_guide_seqs.append(args.prime_editing_pegRNA_extension_seq)
                pe_guide_mismatches.append(best_aln_mismatches)
                pe_guide_names.append('PE Extension')
                pe_guide_qw_centers.append(0)
                pe_guide_qw_sizes.append(prime_editing_pegRNA_extension_quantification_window_size)
                pe_guide_plot_cut_points.append(False)
            #otherwise, clone the coordinates from the prime_editing_edited_amp_seq
            else:
                best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(prime_editing_extension_seq_dna, prime_editing_edited_amp_seq,aln_matrix,needleman_wunsch_gap_open,0)
                match = re.search(best_aln_seq, prime_editing_edited_amp_seq)
                pe_start_loc = match.start()
                pe_end_loc = match.end()
                coords_l, coords_r = CRISPRessoShared.get_alignment_coordinates(to_sequence=this_amp_seq, from_sequence=prime_editing_edited_amp_seq,aln_matrix=aln_matrix,needleman_wunsch_gap_open=needleman_wunsch_gap_open,needleman_wunsch_gap_extend=needleman_wunsch_gap_extend)
                new_seq = this_amp_seq[coords_l[pe_start_loc]:coords_r[pe_end_loc]]
                pe_guides.append(new_seq)
                pe_orig_guide_seqs.append(args.prime_editing_pegRNA_extension_seq)
                rev_coords_l, rev_coords_r = CRISPRessoShared.get_alignment_coordinates(to_sequence=prime_editing_edited_amp_seq, from_sequence=this_amp_seq,aln_matrix=aln_matrix,needleman_wunsch_gap_open=needleman_wunsch_gap_open,needleman_wunsch_gap_extend=needleman_wunsch_gap_extend)

                this_mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(this_amp_seq, prime_editing_edited_amp_seq, pe_start_loc, pe_end_loc, coords_l, coords_r, rev_coords_l, rev_coords_r)
                this_mismatches += [coords_l[i] for i in best_aln_mismatches] #add mismatches to original sequence
                pe_guide_mismatches.append(this_mismatches)

                pe_guide_names.append('PE Extension')
                pe_guide_qw_centers.append(0)
                pe_guide_qw_sizes.append(prime_editing_pegRNA_extension_quantification_window_size)
                pe_guide_plot_cut_points.append(False)


            #now handle the pegRNA spacer seq
            if args.prime_editing_pegRNA_spacer_seq == "":
                raise CRISPRessoShared.BadParameterException('The prime editing pegRNA spacer sequence (--prime_editing_pegRNA_spacer_seq) is required for prime editing analysis.')
            pegRNA_spacer_seq = args.prime_editing_pegRNA_spacer_seq.upper().replace('U', 'T')
            wrong_nt=CRISPRessoShared.find_wrong_nt(pegRNA_spacer_seq)
            if wrong_nt:
                raise CRISPRessoShared.NTException('The prime editing pegRNA spacer sgRNA sequence contains bad characters:%s'  % ' '.join(wrong_nt))
            if pegRNA_spacer_seq not in amplicon_seq_arr[0] and CRISPRessoShared.reverse_complement(pegRNA_spacer_seq) not in amplicon_seq_arr[0]:
                raise CRISPRessoShared.BadParameterException('The given prime editing pegRNA spacer is not found in the reference sequence')

            #spacer is found in the first amplicon (unmodified ref), may be modified in the other amplicons
            #if this is the first sequence, add it directly
            if this_amp_seq == ref0_seq:
                best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(pegRNA_spacer_seq, ref0_seq,aln_matrix,needleman_wunsch_gap_open,0)
                pe_guides.append(best_aln_seq)
                pe_orig_guide_seqs.append(args.prime_editing_pegRNA_spacer_seq)
                pe_guide_mismatches.append(best_aln_mismatches)
                pe_guide_names.append('PE spacer sgRNA')
                pe_guide_qw_centers.append(nicking_qw_center)
                pe_guide_qw_sizes.append(nicking_qw_size)
                pe_guide_plot_cut_points.append(True)
            #otherwise, clone the coordinates from the ref0 amplicon
            else:
                best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(pegRNA_spacer_seq, ref0_seq,aln_matrix,needleman_wunsch_gap_open,0)
                match = re.search(best_aln_seq, ref0_seq)
                r0_start_loc = match.start()
                r0_end_loc = match.end()
                coords_l, coords_r = CRISPRessoShared.get_alignment_coordinates(to_sequence=this_amp_seq, from_sequence=ref0_seq,aln_matrix=aln_matrix,needleman_wunsch_gap_open=needleman_wunsch_gap_open,needleman_wunsch_gap_extend=needleman_wunsch_gap_extend)
                new_seq = this_amp_seq[coords_l[r0_start_loc]:coords_r[r0_end_loc]]
                pe_guides.append(new_seq)
                pe_orig_guide_seqs.append(args.prime_editing_pegRNA_spacer_seq)
                rev_coords_l, rev_coords_r = CRISPRessoShared.get_alignment_coordinates(to_sequence=ref0_seq, from_sequence=this_amp_seq,aln_matrix=aln_matrix,needleman_wunsch_gap_open=needleman_wunsch_gap_open,needleman_wunsch_gap_extend=needleman_wunsch_gap_extend)
                this_mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(this_amp_seq, ref0_seq, r0_start_loc, r0_end_loc, coords_l, coords_r, rev_coords_l, rev_coords_r)
                this_mismatches += [coords_l[i] for i in best_aln_mismatches] #add mismatches to original sequence
                pe_guide_mismatches.append(this_mismatches)

                pe_guide_names.append('PE spacer sgRNA')
                nicking_center_ref0 = r0_end_loc + nicking_qw_center #if there are indels in this amplicon between the end of the guide and the nicking center, adjust the center
                nicking_center_this_amp_seq = rev_coords_r[nicking_center_ref0] - coords_r[r0_end_loc]
                pe_guide_qw_centers.append(nicking_center_this_amp_seq)
                pe_guide_qw_sizes.append(nicking_qw_size)
                pe_guide_plot_cut_points.append(True)

            #nicking guide
            if args.prime_editing_nicking_guide_seq:
                nicking_guide_seq = args.prime_editing_nicking_guide_seq.upper().replace('U', 'T')
                wrong_nt=CRISPRessoShared.find_wrong_nt(nicking_guide_seq)
                if wrong_nt:
                    raise CRISPRessoShared.NTException('The prime editing nicking sgRNA sequence contains bad characters:%s'  % ' '.join(wrong_nt))

                #nicking guide is found in the reverse_complement of the first amplicon, may be modified in the other amplicons
                if this_amp_seq == ref0_seq:
                    rc_ref0_seq = CRISPRessoShared.reverse_complement(ref0_seq)
                    best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(nicking_guide_seq, rc_ref0_seq,aln_matrix,needleman_wunsch_gap_open,0)
                    if nicking_guide_seq not in rc_ref0_seq:
                        warn('The given prime editing nicking guide is not found in the reference sequence. Using the best match: ' + str(best_aln_seq))
                    pe_guides.append(best_aln_seq)
                    pe_orig_guide_seqs.append(args.prime_editing_nicking_guide_seq)
                    pe_guide_mismatches.append(best_aln_mismatches)
                    pe_guide_names.append('PE nicking sgRNA')
                    pe_guide_qw_centers.append(nicking_qw_center)
                    pe_guide_qw_sizes.append(nicking_qw_size)
                    pe_guide_plot_cut_points.append(True)
                #otherwise, clone the coordinates from the ref0 amplicon
                else:
                    rc_ref0_seq = CRISPRessoShared.reverse_complement(ref0_seq)
                    rc_this_amp_seq = CRISPRessoShared.reverse_complement(this_amp_seq)
                    best_aln_seq, best_aln_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2 = CRISPRessoShared.get_best_aln_pos_and_mismatches(nicking_guide_seq, rc_ref0_seq,aln_matrix,needleman_wunsch_gap_open,0)
                    match = re.search(best_aln_seq, rc_ref0_seq)
                    r0_start_loc = match.start()
                    r0_end_loc = match.end()
                    coords_l, coords_r = CRISPRessoShared.get_alignment_coordinates(to_sequence=rc_this_amp_seq, from_sequence=rc_ref0_seq,aln_matrix=aln_matrix,needleman_wunsch_gap_open=needleman_wunsch_gap_open,needleman_wunsch_gap_extend=needleman_wunsch_gap_extend)
                    new_seq = rc_this_amp_seq[coords_l[r0_start_loc]:coords_r[r0_end_loc]]
                    pe_guides.append(new_seq)
                    pe_orig_guide_seqs.append(args.prime_editing_nicking_guide_seq)
                    rev_coords_l, rev_coords_r = CRISPRessoShared.get_alignment_coordinates(to_sequence=rc_ref0_seq, from_sequence=rc_this_amp_seq,aln_matrix=aln_matrix,needleman_wunsch_gap_open=needleman_wunsch_gap_open,needleman_wunsch_gap_extend=needleman_wunsch_gap_extend)
                    this_mismatches = CRISPRessoShared.get_sgRNA_mismatch_vals(rc_this_amp_seq, rc_ref0_seq, r0_start_loc, r0_end_loc, coords_l, coords_r, rev_coords_l, rev_coords_r)
                    this_mismatches += [coords_l[i] for i in best_aln_mismatches] #add mismatches to original sequence

                    pe_guide_mismatches.append(this_mismatches)
                    pe_guide_names.append('PE nicking sgRNA')
                    nicking_center_ref0 = r0_end_loc + nicking_qw_center
                    nicking_center_this_amp_seq = rev_coords_r[nicking_center_ref0] - coords_r[r0_end_loc]
                    pe_guide_qw_centers.append(nicking_center_this_amp_seq)
                    pe_guide_qw_sizes.append(nicking_qw_size)
                    pe_guide_plot_cut_points.append(True)

            return(pe_guides, pe_orig_guide_seqs, pe_guide_mismatches, pe_guide_names, pe_guide_qw_centers, pe_guide_qw_sizes, pe_guide_plot_cut_points)
        #end prime editing guide function

        #now that we're done with adding possible guides and amplicons, go through each amplicon and compute quantification windows
        found_guide_seq = [False]*len(guides)
        found_coding_seq = [False]*len(coding_seqs)

        max_amplicon_len = 0 #for flash
        min_amplicon_len = 99**99 #for flash

        for idx, seq in enumerate(amplicon_seq_arr):
            this_seq = seq.strip().upper()
            this_seq_length = len(this_seq)
            if this_seq_length > max_amplicon_len:
                max_amplicon_len = this_seq_length
            if this_seq_length < min_amplicon_len:
                min_amplicon_len = this_seq_length

            this_name = 'Amplicon'+str(idx)
            if idx < len(amplicon_name_arr):
                this_name = amplicon_name_arr[idx]

            this_min_aln_score = args.default_min_aln_score
            if idx < len(amplicon_min_alignment_score_arr):
                this_min_aln_score = amplicon_min_alignment_score_arr[idx]

            this_quant_window_coordinates = None
            if amplicon_quant_window_coordinates_arr[idx] != "":
                this_quant_window_coordinates = amplicon_quant_window_coordinates_arr[idx]

            this_guides = guides[:]
            this_orig_guide_seqs = guides[:]
            this_guide_mismatches = [[]]*len(guides) #all guides have no mismatches up to this point
            this_guide_names = guide_names[:]
            this_guide_qw_centers = guide_qw_centers[:]
            this_guide_qw_sizes = guide_qw_sizes[:]
            this_guide_plot_cut_points = guide_plot_cut_points[:]

            #for this amplicon, add flexi guides (that can contain mismatches)
            if args.flexiguide_seq:
                flexi_guides = []
                flexi_guide_orig_seqs = []
                flexi_guide_mismatches = []
                all_flexiguide_names = args.flexiguide_name.split(",")
                flexi_guide_names = []
                for idx, guide in enumerate(args.flexiguide_seq.split(",")):
                    #for all amps in forward and reverse complement amps:
                    for amp_seq in [this_seq, CRISPRessoShared.reverse_complement(this_seq)]:
                        ref_incentive = np.zeros(len(amp_seq)+1, dtype=int)
                        s1, s2, score=CRISPResso2Align.global_align(guide, amp_seq, matrix=aln_matrix, gap_incentive=ref_incentive, gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
                        potential_guide = s1.strip("-")
                        if abs(len(potential_guide) - len(guide)) < 2: #if length of putative guide is off by less than 2, keep it (allows 1 gap)
                            loc = s1.find(potential_guide)
                            potential_ref = amp_seq[loc:loc+len(potential_guide)]
                            #realign to test for number of mismatches
                            ref_incentive = np.zeros(len(potential_ref)+1, dtype=int)
                            sub_s1, sub_s2, sub_score=CRISPResso2Align.global_align(guide, potential_ref, matrix=aln_matrix, gap_incentive=ref_incentive, gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend,)
                            mismatches = []
                            for i in range(len(sub_s1)):
                                if sub_s1[i] != sub_s2[i]:
                                    mismatches.append(i)

                            if sub_score > args.flexiguide_homology:
                                flexi_guides.append(potential_ref)
                                flexi_guide_orig_seqs.append(guide)
                                flexi_guide_mismatches.append(mismatches)
                                this_flexiguide_name = ''
                                if idx < len (all_flexiguide_names) and all_flexiguide_names[idx] != '':
                                    this_flexiguide_name = all_flexiguide_names[idx]
                                flexi_guide_names.append(this_flexiguide_name)

                flexi_guide_count = 0
                for i in range(len(flexi_guides)):
                    flexi_guide = flexi_guides[i]
                    if flexi_guide not in guides:
                        flexi_guide_count += 1
                        this_guides.append(flexi_guide)
                        this_orig_guide_seqs.append(flexi_guide_orig_seqs[i])
                        this_guide_mismatches.append(flexi_guide_mismatches[i])
                        this_guide_names.append(flexi_guide_names[i])
                        this_guide_qw_centers.append(int(args.quantification_window_center.split(",")[0]))
                        this_guide_qw_sizes.append(int(args.quantification_window_size.split(",")[0]))
                        if args.base_editor_output:
                            this_guide_plot_cut_points.append(False)
                        else:
                            this_guide_plot_cut_points.append(True)
                info('Added %d guides with flexible matching\n\tOriginal flexiguides: %s\n\tFound guides: %s\n\tMismatch locations: %s'%(flexi_guide_count, str(args.flexiguide_seq.split(",")), str(flexi_guides), str(flexi_guide_mismatches)))

            if args.prime_editing_pegRNA_extension_seq:
                nicking_qw_center = int(args.quantification_window_center.split(",")[0])
                nicking_qw_size = int(args.quantification_window_size.split(",")[0])
                pe_guides, pe_orig_guide_seqs, pe_guide_mismatches, pe_guide_names, pe_guide_qw_centers, pe_guide_qw_sizes, pe_guide_plot_cut_points = get_prime_editing_guides(this_seq, this_name, amplicon_seq_arr[0],
                        prime_editing_edited_amp_seq, prime_editing_extension_seq_dna, args.prime_editing_pegRNA_extension_quantification_window_size, nicking_qw_center, nicking_qw_size, aln_matrix,args.needleman_wunsch_gap_open,args.needleman_wunsch_gap_extend)
                this_guides.extend(pe_guides)
                this_orig_guide_seqs.extend(pe_orig_guide_seqs)
                this_guide_mismatches.extend(pe_guide_mismatches)
                this_guide_names.extend(pe_guide_names)
                this_guide_qw_centers.extend(pe_guide_qw_centers)
                this_guide_qw_sizes.extend(pe_guide_qw_sizes)
                this_guide_plot_cut_points.extend(pe_guide_plot_cut_points)


            # Calculate cut sites for this reference
            (this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_sgRNA_names, this_include_idxs,
                this_exclude_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(this_seq, this_guides, this_guide_mismatches, this_guide_names, this_guide_qw_centers,
                this_guide_qw_sizes, this_quant_window_coordinates, args.exclude_bp_from_left, args.exclude_bp_from_right, args.plot_window_size, this_guide_plot_cut_points, args.discard_guide_positions_overhanging_amplicon_edge)

            this_orig_guide_lookup = {} #dict of new seq to original (input) sequence
            for idx, guide in enumerate(this_guides):
                this_orig_guide_lookup[guide.upper()] = this_orig_guide_seqs[idx]
            this_sgRNA_orig_seqs = []
            for seq in this_sgRNA_sequences:
                this_sgRNA_orig_seqs.append(this_orig_guide_lookup[seq])

            this_contains_guide = False
            if len(this_sgRNA_sequences) > 0:
                this_contains_guide = True

            for guide_idx, guide_seq in enumerate(guides):
                #this_sgRNA_sequences contains all good guides contained in this amplicon
                #don't consider flexiguides or other guides in this_sgRNA_sequences which may be specific to this amplicon
                if guide_seq.upper() in this_sgRNA_sequences:
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

            this_gap_incentive = np.zeros(this_seq_length+1, dtype=int)
            for cut_point in this_sgRNA_cut_points:
                this_gap_incentive[cut_point+1] = args.needleman_wunsch_gap_incentive

            seq_rc = CRISPRessoShared.reverse_complement(this_seq)
            seeds = []
            rc_seeds = []
            seedStarts = list(range(args.exclude_bp_from_left, this_seq_length-args.exclude_bp_from_right-args.aln_seed_len, args.aln_seed_count)) #define all possible seed starts
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

            refObj = {'name': this_name,
                   'sequence': this_seq,
                   'sequence_length': this_seq_length,
                   'min_aln_score': this_min_aln_score,
                   'gap_incentive': this_gap_incentive,
                   'sgRNA_cut_points': this_sgRNA_cut_points,
                   'sgRNA_intervals': this_sgRNA_intervals,
                   'sgRNA_sequences': this_sgRNA_sequences,
                   'sgRNA_plot_cut_points': this_sgRNA_plot_cut_points,
                   'sgRNA_plot_idxs': this_sgRNA_plot_idxs,
                   'sgRNA_names': this_sgRNA_names,
                   'sgRNA_mismatches': this_sgRNA_mismatches,
                   'sgRNA_orig_sequences': this_sgRNA_orig_seqs,
                   'contains_guide': this_contains_guide,
                   'contains_coding_seq': this_contains_coding_seq,
                   'exon_positions': this_exon_positions,
                   'exon_len_mods': this_exon_len_mods,
                   'exon_intervals': this_exon_intervals,
                   'splicing_positions': this_splicing_positions,
                   'include_idxs': this_include_idxs,
                   'exclude_idxs': this_exclude_idxs,
                   'idx_cloned_from': None,
                   'fw_seeds': seeds,
                   'rc_seeds': rc_seeds,
                   'aln_genome': None,
                   'aln_chr': None,
                   'aln_start': None,
                   'aln_end': None,
                   'aln_strand': None
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
                    raise CRISPRessoShared.ExonSequenceException('The coding subsequence %d (%s) provided is not contained in any amplicon sequence!\n\nPlease check your input!' % (idx, coding_seqs[idx]))

        #clone cut points and include idx from first reference where those are set (also exons)
        clone_ref_name = None
        clone_ref_idx = None
        clone_has_cut_points = False
        clone_has_exons = False
        for idx, ref_name in enumerate(ref_names):
            cut_points = refs[ref_name]['sgRNA_cut_points']
            exon_positions = refs[ref_name]['exon_positions']
            if cut_points or (idx < len(amplicon_quant_window_coordinates_arr) and amplicon_quant_window_coordinates_arr[idx] != ""):
                if len(ref_names) > 1 and args.debug:
                    info("Using cut points from %s as template for other references"%ref_name)
                clone_ref_name = ref_name
                clone_ref_idx = idx
                clone_has_cut_points = True
                if exon_positions:
                    clone_has_exons = True
                break
            if exon_positions:
                if len(ref_names) > 1 and args.debug:
                    info("Using exons positions from %s as template for other references"%ref_name)
                clone_ref_name = ref_name
                clone_ref_idx = idx
                clone_has_exon_positions = True
                if cut_points:
                    clone_has_cut_points = True
                break

        if clone_ref_name is not None:
            for ref_name in ref_names:
                if clone_ref_name == ref_name:
                    continue
                cut_points = refs[ref_name]['sgRNA_cut_points']
                sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
                exon_positions = refs[ref_name]['exon_positions']

                needs_cut_points = False
                needs_sgRNA_intervals = False
                needs_exon_positions = False

                #if HDR, copy everything from ref 1 so the quantification window will be accurate
                if ref_name == "HDR" and args.expected_hdr_amplicon_seq != "":
                    needs_cut_points = True
                    needs_sgRNA_intervals = True
                    needs_exon_positions = True
                else:
                    if cut_points:
                        if len(ref_names) > 1 and args.debug:
                            info("Reference '%s' has cut points defined: %s. Not inferring."%(ref_name, cut_points))
                    else:
                        needs_cut_points = True

                    if sgRNA_intervals:
                        if len(ref_names) > 1 and args.debug:
                            info("Reference '%s' has sgRNA_intervals defined: %s. Not inferring."%(ref_name, sgRNA_intervals))
                    else:
                        needs_sgRNA_intervals = True

                    if exon_positions:
                        if len(exon_positions) > 1 and args.debug:
                            info("Reference '%s' has exon_positions defined: %s. Not inferring."%(ref_name, exon_positions))
                    else:
                        needs_exon_positions = True


                if not needs_cut_points and not needs_sgRNA_intervals and not needs_exon_positions and args.quantification_window_coordinates is None:
                    continue

                fws1, fws2, fwscore=CRISPResso2Align.global_align(refs[ref_name]['sequence'], refs[clone_ref_name]['sequence'], matrix=aln_matrix, gap_open=args.needleman_wunsch_gap_open, gap_extend=args.needleman_wunsch_gap_extend, gap_incentive=refs[clone_ref_name]['gap_incentive'])
                if fwscore < 60:
                    continue

                if (needs_sgRNA_intervals or needs_cut_points) and clone_has_cut_points and args.debug:
                    info("Reference '%s' has NO cut points or sgRNA intervals idxs defined. Inferring from '%s'."%(ref_name, clone_ref_name))
                if needs_exon_positions and clone_has_exons and args.debug:
                    info("Reference '%s' has NO exon_positions defined. Inferring from '%s'."%(ref_name, clone_ref_name))
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
                    this_gap_incentive = np.zeros(refs[ref_name]['sequence_length']+1, dtype=int)
                    for cut_point in this_cut_points:
                        this_gap_incentive[cut_point + 1] = args.needleman_wunsch_gap_incentive

                    this_sgRNA_intervals = []
                    this_sgRNA_mismatches = []
                    this_sgRNA_plot_idxs = []
                    for idx, (sgRNA_interval_start, sgRNA_interval_end) in enumerate(refs[clone_ref_name]['sgRNA_intervals']):
                        this_sgRNA_intervals.append((s1inds[sgRNA_interval_start], s1inds[sgRNA_interval_end]))

                        sgRNA_cut_point_this = s1inds[refs[clone_ref_name]['sgRNA_cut_points'][idx]]
                        sgRNA_seq_clone = refs[clone_ref_name]['sequence'][sgRNA_interval_start:sgRNA_interval_end+1]
                        sgRNA_seq_this = refs[ref_name]['sequence'][s1inds[sgRNA_interval_start]:s1inds[sgRNA_interval_end]+1]
                        sgRNA_seq_orig_seq = refs[clone_ref_name]['sgRNA_orig_sequences'][idx]
                        fw_mismatches = CRISPRessoShared.get_mismatches(sgRNA_seq_this, sgRNA_seq_orig_seq, aln_matrix, args.needleman_wunsch_gap_open,args.needleman_wunsch_gap_extend)
                        rv_mismatches = CRISPRessoShared.get_mismatches(sgRNA_seq_this, CRISPRessoShared.reverse_complement(sgRNA_seq_orig_seq), aln_matrix, args.needleman_wunsch_gap_open,args.needleman_wunsch_gap_extend)
                        best_aln_mismatches = fw_mismatches
                        if len(rv_mismatches) < len(fw_mismatches):
                            best_aln_mismatches = rv_mismatches
                        this_sgRNA_mismatches.append(sorted(best_aln_mismatches))

                        ref_seq_length = len(refs[ref_name]['sequence'])
                        window_around_cut=max(1, args.plot_window_size)
                        st=max(0, sgRNA_cut_point_this-window_around_cut+1)
                        en=min(ref_seq_length-1, sgRNA_cut_point_this+window_around_cut+1)
                        this_sgRNA_plot_idxs.append(sorted(list(range(st, en))))

                    old_this_sgRNA_plot_idxs = []
                    for plot_idx_list in refs[clone_ref_name]['sgRNA_plot_idxs']:
                        old_this_sgRNA_plot_idxs.append([s1inds[x] for x in plot_idx_list])

                    this_include_idxs = []
                    start_val = -1
                    last_val = -1
                    for idx, val in enumerate(refs[clone_ref_name]['include_idxs']):
                        if start_val == -1:
                            start_val = val
                        elif val != last_val + 1:
                            this_include_idxs.extend(list(range(s1inds[start_val], s1inds[last_val]+1)))
                            start_val = val
                        last_val = val

                    if start_val != -1:
                        this_include_idxs.extend(list(range(s1inds[start_val], s1inds[last_val]+1)))

                    #subtract any indices in 'exclude_idxs' -- e.g. in case some of the cloned include_idxs were near the read ends (excluded)
                    this_exclude_idxs = sorted(list(set(refs[ref_name]['exclude_idxs'])))
                    this_include_idxs = sorted(list(set(np.setdiff1d(this_include_idxs, this_exclude_idxs))))

                    refs[ref_name]['gap_incentive'] = this_gap_incentive
                    refs[ref_name]['sgRNA_cut_points'] = this_cut_points
                    refs[ref_name]['sgRNA_intervals'] = this_sgRNA_intervals
                    refs[ref_name]['sgRNA_sequences'] = refs[clone_ref_name]['sgRNA_sequences']
                    refs[ref_name]['sgRNA_plot_cut_points'] = refs[clone_ref_name]['sgRNA_plot_cut_points']
                    refs[ref_name]['sgRNA_plot_idxs'] = this_sgRNA_plot_idxs
                    refs[ref_name]['sgRNA_mismatches'] = this_sgRNA_mismatches
                    refs[ref_name]['sgRNA_orig_sequences'] = refs[clone_ref_name]['sgRNA_orig_sequences']
                    refs[ref_name]['sgRNA_names'] = refs[clone_ref_name]['sgRNA_names']
                    refs[ref_name]['include_idxs'] = this_include_idxs
                    refs[ref_name]['contains_guide'] = refs[clone_ref_name]['contains_guide']

                #quantification window coordinates override other options
                if amplicon_quant_window_coordinates_arr[clone_ref_idx] != "":
                    this_include_idxs = []
                    these_coords = amplicon_quant_window_coordinates_arr[clone_ref_idx].split("_")
                    for coord in these_coords:
                        coordRE = re.match(r'^(\d+)-(\d+)$', coord)
                        if coordRE:
                            start = s1inds[int(coordRE.group(1))]
                            end = s1inds[int(coordRE.group(2)) + 1]
                            this_include_idxs.extend(range(start, end))
                        else:
                            raise NTException("Cannot parse analysis window coordinate '" + str(coord))
                    #subtract any indices in 'exclude_idxs' -- e.g. in case some of the cloned include_idxs were near the read ends (excluded)
                    this_exclude_idxs = sorted(list(set(refs[ref_name]['exclude_idxs'])))
                    this_include_idxs = sorted(list(set(np.setdiff1d(this_include_idxs, this_exclude_idxs))))
                    refs[ref_name]['include_idxs'] = this_include_idxs
                    refs[ref_name]['exclude_idxs'] = this_exclude_idxs

                if needs_exon_positions and clone_has_exons:
                    this_exon_positions = set()
                    this_splicing_positions = []
                    this_exon_intervals = []
                    this_exon_len_mods = []
                    this_seq_length = refs[ref_name]['sequence_length']
                    for (exon_interval_start, exon_interval_end) in refs[clone_ref_name]['exon_intervals']:
                        this_exon_start = s1inds[exon_interval_start]
                        this_exon_end = s1inds[exon_interval_end -1] + 1 # in case where exon is entire amplicon, because the exon_end is beyond the upper bound, we have to make this change
                        this_exon_intervals.append((this_exon_start, this_exon_end))
                        this_exon_len_mods.append(((this_exon_end - this_exon_start)-(exon_interval_end - exon_interval_start)))
                        this_exon_positions = this_exon_positions.union(set(range(this_exon_start, this_exon_end)))
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

        # if we're writing a bam, find the alignment location
        bam_header = '@HD\tVN:1.0\tSO:unknown\n'
        if args.bam_output:
            if args.bowtie2_index == '':
                bam_output_ref_fa = _jp('CRISPResso_output.fa')
                with open(bam_output_ref_fa, 'w') as fa_out:
                    for ref_name in ref_names:
                        refs[ref_name]['aln_genome'] = 'None'
                        refs[ref_name]['aln_chr'] = ref_name
                        refs[ref_name]['aln_start'] = 1
                        refs[ref_name]['aln_end'] = refs[ref_name]['sequence_length']
                        refs[ref_name]['aln_strand'] = '+'
                        bam_header += '@SQ\tSN:%s\tLN:%s\n'%(ref_name, refs[ref_name]['sequence_length'])
                        fa_out.write('>%s\n%s\n'%(ref_name,refs[ref_name]['sequence']))
            else:
                check_program('bowtie2')

                if (os.path.isfile(args.bowtie2_index+'.1.bt2l')):
                    CRISPRessoShared.check_file(args.bowtie2_index+'.1.bt2l')
                else:
                    CRISPRessoShared.check_file(args.bowtie2_index+'.1.bt2')

                filename_aligned_amplicons_sam = _jp('CRISPResso_amplicons_aligned.sam')
                filename_aligned_amplicons_sam_log = _jp('CRISPResso_amplicons_aligned.sam.log')
                filename_amplicon_seqs_fasta = _jp('CRISPResso_amplicons_to_align.fa')

                files_to_remove.append(filename_aligned_amplicons_sam)
                files_to_remove.append(filename_amplicon_seqs_fasta)

                with open(filename_amplicon_seqs_fasta, 'w') as fastas:
                    for ref_name in ref_names:
                        ref_seq = refs[ref_name]['sequence']
                        fastas.write('>%s\n%s\n'%(ref_name, ref_seq))

                aligner_command = 'bowtie2 -x %s -p %s -f -U %s 2> %s > %s ' %(args.bowtie2_index, args.n_processes, \
                                  filename_amplicon_seqs_fasta, filename_aligned_amplicons_sam_log, filename_aligned_amplicons_sam)
                bowtie_status = sb.call(aligner_command, shell=True)
                if bowtie_status:
                    raise CRISPRessoShared.BadParameterException('Bowtie2 failed to align amplicons to the genome, please check the output file and log here: %s.'% filename_aligned_amplicons_sam_log)

                bam_header = ''
                with open(filename_aligned_amplicons_sam) as aln:
                    for line in aln.readlines():
                        if line.startswith('@'):
                            bam_header += line
                            continue
                        line_els = line.split("\t")
                        if line_els[2] == "*":
                            raise CRISPRessoShared.BadParameterException('The amplicon [%s] is not mappable to the reference genome provided!' % line_els[0])

                        else:
                            aln_len = CRISPRessoShared.get_ref_length_from_cigar(line_els[5])
                            seq_start = int(line_els[3])
                            seq_stop = seq_start + aln_len
                            strand = "-" if (int(line_els[1]) & 0x10) else "+"
                            info('The amplicon [%s] was mapped to: %s:%d-%d ' % (line_els[0], line_els[2], seq_start, seq_stop))
                        ref_name = line_els[0]
                        refs[ref_name]['aln_genome'] = args.bowtie2_index
                        refs[ref_name]['aln_chr'] = line_els[2]
                        refs[ref_name]['aln_start'] = seq_start
                        refs[ref_name]['aln_end'] = seq_stop
                        refs[ref_name]['aln_strand'] = strand

        N_READS_INPUT = 0
        if args.fastq_r1:
            N_READS_INPUT = get_n_reads_fastq(args.fastq_r1)
        elif args.bam_input:
            N_READS_INPUT = get_n_reads_bam(args.bam_input, args.bam_chr_loc)

        if N_READS_INPUT == 0:
            raise CRISPRessoShared.BadParameterException('The input contains 0 reads.')

        args_string_arr = []
        clean_args_string_arr = []
        for arg in vars(args):
            val = str(getattr(args, arg))
            args_string_arr.append("%s: %s"%(str(arg), val))
            if os.sep in val:
                val = os.path.basename(val)
            clean_args_string_arr.append("%s: %s"%(str(arg), val))

        if args.write_cleaned_report:
            crispresso2_info['running_info']['args_string'] = '\n'.join(sorted(clean_args_string_arr))
        else:
            crispresso2_info['running_info']['args_string'] = '\n'.join(sorted(args_string_arr))

        crispresso2_info['running_info']['start_time'] = start_time
        crispresso2_info['running_info']['start_time_string'] = start_time_string

        if args.split_interleaved_input:
            if args.fastq_r2!='' or args.bam_input != '':
                raise CRISPRessoShared.BadParameterException('The option --split_interleaved_input is available only when a single fastq file is specified!')
            else:
                info('Splitting paired end single fastq file into two files...')
                args.fastq_r1, args.fastq_r2=split_interleaved_fastq(args.fastq_r1,
                    output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq', '')).replace('.gz', '')+'_splitted_r1.fastq.gz'),
                    output_filename_r2=_jp(os.path.basename(args.fastq_r1.replace('.fastq', '')).replace('.gz', '')+'_splitted_r2.fastq.gz'),)
                files_to_remove+=[args.fastq_r1, args.fastq_r2]
                N_READS_INPUT = N_READS_INPUT/2

                info('Done!')

        if args.min_average_read_quality>0 or args.min_single_bp_quality>0 or args.min_bp_quality_or_N>0:
            if args.bam_input != '':
                raise CRISPRessoShared.BadParameterException('The read filtering options are not available with bam input')
            info('Filtering reads with average bp quality < %d and single bp quality < %d and replacing bases with quality < %d with N ...' % (args.min_average_read_quality, args.min_single_bp_quality, args.min_bp_quality_or_N))
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
                output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq', '')).replace('.gz', '')+'_filtered.fastq.gz')
                output_filename_r2=_jp(os.path.basename(args.fastq_r2.replace('.fastq', '')).replace('.gz', '')+'_filtered.fastq.gz')

                from CRISPResso2 import filterFastqs
                filterFastqs.filterFastqs(fastq_r1=args.fastq_r1, fastq_r2=args.fastq_r2, fastq_r1_out=output_filename_r1, fastq_r2_out=output_filename_r2, min_bp_qual_in_read=min_single_bp_quality, min_av_read_qual=min_av_quality, min_bp_qual_or_N=min_bp_quality_or_N)

                args.fastq_r1 = output_filename_r1
                args.fastq_r2 = output_filename_r2
                files_to_remove += [output_filename_r1]
                files_to_remove += [output_filename_r2]

            else:
                output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq', '')).replace('.gz', '')+'_filtered.fastq.gz')

                from CRISPResso2 import filterFastqs
                filterFastqs.filterFastqs(fastq_r1=args.fastq_r1, fastq_r1_out=output_filename_r1, min_bp_qual_in_read=min_single_bp_quality, min_av_read_qual=min_av_quality, min_bp_qual_or_N=min_bp_quality_or_N)

                args.fastq_r1 = output_filename_r1
                files_to_remove += [output_filename_r1]


        #Trim and merge reads
        if args.bam_input != '' and args.trim_sequences:
            raise CRISPRessoShared.BadParameterException('Read trimming options are not available with bam input')
        elif args.fastq_r1 != '' and args.fastq_r2 == '': #single end reads
            if not args.trim_sequences: #no trimming or merging required
                output_forward_filename=args.fastq_r1
            else:
                output_forward_filename=_jp('reads.trimmed.fq.gz')
                #Trimming with trimmomatic
                cmd='%s SE -phred33 %s  %s %s >>%s 2>&1'\
                % (args.trimmomatic_command, args.fastq_r1,
                   output_forward_filename,
                   args.trimmomatic_options_string.replace('NexteraPE-PE.fa', 'TruSeq3-SE.fa'),
                   log_filename)
                #print cmd
                TRIMMOMATIC_STATUS=sb.call(cmd, shell=True)

                if TRIMMOMATIC_STATUS:
                        raise CRISPRessoShared.TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')
                crispresso2_info['trimmomatic_command'] = cmd

                files_to_remove += [output_forward_filename]

            processed_output_filename=output_forward_filename

        elif args.fastq_r1 != '' and args.fastq_r2 != '':#paired end reads
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
                        args.fastq_r1, args.fastq_r2, output_forward_paired_filename,
                        output_forward_unpaired_filename, output_reverse_paired_filename,
                        output_reverse_unpaired_filename, args.trimmomatic_options_string, log_filename)
                #print cmd
                TRIMMOMATIC_STATUS=sb.call(cmd, shell=True)
                if TRIMMOMATIC_STATUS:
                    raise CRISPRessoShared.TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')
                crispresso2_info['trimmomatic_command'] = cmd

                files_to_remove += [output_forward_paired_filename]
                files_to_remove += [output_reverse_paired_filename]

                info('Done!')

            #for paired-end reads, merge them
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
                max_overlap = max(10, min(avg_read_length, expected_max_overlap+indel_overlap_tolerance))
                #min overlap is either 4bp (as in crispresso1) or the expected amplicon length - indel tolerance
                min_overlap = max(4, expected_min_overlap-indel_overlap_tolerance)
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
            cmd='%s "%s" "%s" --min-overlap %d --max-overlap %d --allow-outies -z -d %s -o %s >>%s 2>&1' %\
            (args.flash_command,
                 output_forward_paired_filename,
                 output_reverse_paired_filename,
                 min_overlap,
                 max_overlap,
                 OUTPUT_DIRECTORY,
                 output_prefix,
                 log_filename)

            info('Running FLASH command: ' + cmd)
            crispresso2_info['flash_command'] = cmd
            FLASH_STATUS=sb.call(cmd, shell=True)
            if FLASH_STATUS:
                raise CRISPRessoShared.FlashException('Flash failed to run, please check the log file.')

            flash_hist_filename=_jp('out.hist')
            flash_histogram_filename=_jp('out.histogram')
            flash_not_combined_1_filename=_jp('out.notCombined_1.fastq.gz')
            flash_not_combined_2_filename=_jp('out.notCombined_2.fastq.gz')

            processed_output_filename=_jp('out.extendedFrags.fastq.gz')
            if os.path.isfile(processed_output_filename) is False:
                raise CRISPRessoShared.FlashException('Flash failed to produce merged reads file, please check the log file.')

            files_to_remove+=[processed_output_filename, flash_hist_filename, flash_histogram_filename,\
                    flash_not_combined_1_filename, flash_not_combined_2_filename, _jp('out.hist.innie'), _jp('out.histogram.innie'), _jp('out.histogram.outie'), _jp('out.hist.outie')]

            if (args.force_merge_pairs):
                 old_flashed_filename = processed_output_filename
                 new_merged_filename=_jp('out.forcemerged_uncombined.fastq.gz')
                 num_reads_force_merged = CRISPRessoShared.force_merge_pairs(flash_not_combined_1_filename, flash_not_combined_2_filename, new_merged_filename)
                 new_output_filename=_jp('out.forcemerged.fastq.gz')
                 merge_command = "cat %s %s > %s"%(processed_output_filename, new_merged_filename, new_output_filename)
                 MERGE_STATUS=sb.call(merge_command, shell=True)
                 if MERGE_STATUS:
                     raise FlashException('Force-merging read pairs failed to run, please check the log file.')
                 processed_output_filename = new_output_filename

                 files_to_remove+=[new_merged_filename]
                 files_to_remove+=[new_output_filename]
                 if args.debug:
                     info('Wrote force-merged reads to ' + new_merged_filename)

            info('Done!')


        #count reads
        N_READS_AFTER_PREPROCESSING = 0
        if args.bam_input:
            N_READS_AFTER_PREPROCESSING = N_READS_INPUT
        else:
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
        if args.bam_input:
            aln_stats = process_bam(args.bam_input, args.bam_chr_loc, crispresso2_info['bam_output'], variantCache, ref_names, refs, args)
        elif args.fastq_output:
            aln_stats = process_fastq_write_out(processed_output_filename, crispresso2_info['fastq_output'], variantCache, ref_names, refs, args)
        elif args.bam_output:
            bam_header += '@PG\tID:crispresso2\tPN:crispresso2\tVN:'+CRISPRessoShared.__version__+'\tCL:"'+crispresso_cmd_to_write+'"\n'
            aln_stats = process_single_fastq_write_bam_out(processed_output_filename, crispresso2_info['bam_output'], bam_header, variantCache, ref_names, refs, args)
        else:
            aln_stats = process_fastq(processed_output_filename, variantCache, ref_names, refs, args)

        info('Done!')

        if args.prime_editing_pegRNA_scaffold_seq != "":
            #introduce a new ref (that we didn't align to) called 'Scaffold Incorporated' -- copy it from the ref called 'prime-edited'
            new_ref = deepcopy(refs['Prime-edited'])
            new_ref['name'] = "Scaffold-incorporated"
            ref_names.append("Scaffold-incorporated")
            refs["Scaffold-incorporated"] = new_ref

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

        inserted_n_dicts = {} # dict of number of insertions for all reads: dict{len}->number of reads with that insertion length
        deleted_n_dicts = {}
        substituted_n_dicts = {}
        effective_len_dicts = {} # dict of effective lengths for all reads

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


            inserted_n_dicts                    [ref_name] = defaultdict(int)
            deleted_n_dicts                     [ref_name] = defaultdict(int)
            substituted_n_dicts                 [ref_name] = defaultdict(int)
            effective_len_dicts                 [ref_name] = defaultdict(int)

            hists_inframe                       [ref_name] = defaultdict(int)
            hists_inframe                       [ref_name][0] = 0
            hists_frameshift                    [ref_name] = defaultdict(int)
            hists_frameshift                    [ref_name][0] = 0
        #end initialize data structures for each ref
        def get_allele_row(reference_name, variant_count, aln_ref_names_str, aln_ref_scores_str, variant_payload, write_detailed_allele_table):
            """
            gets a row for storing allele information in the allele table
            parameters:
                reference_name: string Reference name to write
                variant_count: number of times this allele appears
                aln_ref_names_str: '&'-joined string of references this allele aligned to
                aln_ref_scores_str: '&'-joined string of references this allele aligned to
                variant_payload: payload object (dict) with keys containing information about the allele
                write_detailed_allele_table: bool for whether to write detailed row
            returns:
                row to put into allele table
            """
            if args.write_detailed_allele_table:
                allele_row = {'#Reads':variant_count,
                    'Aligned_Sequence': variant_payload['aln_seq'],
                    'Reference_Sequence':variant_payload['aln_ref'],
                    'n_inserted':variant_payload['insertion_n'],
                    'n_deleted':variant_payload['deletion_n'],
                    'n_mutated':variant_payload['substitution_n'],
                    'Reference_Name':reference_name,
                    'Read_Status':variant_payload['classification'],
                    'Aligned_Reference_Names':aln_ref_names_str,
                    'Aligned_Reference_Scores':aln_ref_scores_str,
                    'ref_positions':variant_payload['ref_positions'],
                    'all_insertion_positions': variant_payload['all_insertion_positions'], #arr with 1's where there are insertions (including those outside of include_idxs quantification window)
                    'all_insertion_left_positions': variant_payload['all_insertion_left_positions'], #arr with 1's to the left of where the insertion occurs
                    'insertion_positions': variant_payload['insertion_positions'], # arr with 1's where there are insertions (1bp before and 1bp after insertion) that overlap with include_idxs quantification window
                    'insertion_coordinates': variant_payload['insertion_coordinates'], # one entry per insertion, tuple of (start,end)
                    'insertion_sizes': variant_payload['insertion_sizes'],
                    'all_deletion_positions': variant_payload['all_deletion_positions'], #arr with 1's where there are insertions
                    'deletion_positions': variant_payload['deletion_positions'], #arr with 1's where there are insertions that overlap the include_idxs quantification window
                    'deletion_coordinates': variant_payload['deletion_coordinates'], # one entry per deletion
                    'deletion_sizes': variant_payload['deletion_sizes'], # correspond to entries in 'deletion_coordinates'
                    'all_substitution_positions': variant_payload['all_substitution_positions'],
                    'substitution_positions': variant_payload['substitution_positions'],
                    'substitution_values': variant_payload['substitution_values']
    			}
            else:
                allele_row = {'#Reads':variant_count,
                       'Aligned_Sequence': variant_payload['aln_seq'],
                       'Reference_Sequence':variant_payload['aln_ref'],
                       'n_inserted':variant_payload['insertion_n'],
                       'n_deleted':variant_payload['deletion_n'],
                       'n_mutated':variant_payload['substitution_n'],
                       'Reference_Name':reference_name,
                       'Read_Status':variant_payload['classification'],
                       'Aligned_Reference_Names':aln_ref_names_str,
                       'Aligned_Reference_Scores':aln_ref_scores_str,
                       'ref_positions':variant_payload['ref_positions']
    			}
            return allele_row

        #end get_allele_row() definition

        #take care of empty seqs
        cache_fastq_seq = ''
        variantCache[cache_fastq_seq]['count'] = 0

        ###iterate through variants
        for variant in variantCache:
            #skip variant if there were none observed
            variant_count = variantCache[variant]['count']
            if (variant_count == 0):
                continue

            #check to see if this sequence's reverse complement is in the variant
            rc_variant = CRISPRessoShared.reverse_complement(variant)
            if rc_variant in variantCache and variantCache[rc_variant]['count'] > 0:
                variant_count += variantCache[rc_variant]['count']
                variantCache[rc_variant]['count'] = 0
                variantCache[variant]['count'] = variant_count
            N_TOTAL += variant_count

            aln_ref_names = variantCache[variant]['aln_ref_names'] #list of references this seq aligned to
            aln_ref_names_str = '&'.join(aln_ref_names)
            aln_ref_scores = variantCache[variant]['aln_scores']
            aln_ref_scores_str = '&'.join([str(x) for x in aln_ref_scores])
            class_name = variantCache[variant]['class_name'] #for classifying read e.g. 'HDR_MODIFIED' for pie chart

            if class_name not in class_counts:
                    class_counts[class_name] = 0
            class_counts[class_name]+=variant_count

            #if class is AMBIGUOUS (set above if the args.expand_ambiguous_alignments param is false) don't add the modifications in this allele to the allele summaries
            if class_name == "AMBIGUOUS":
                variant_payload = variantCache[variant]["variant_"+aln_ref_names[0]]
                allele_row = get_allele_row('AMBIGUOUS_'+aln_ref_names[0], variant_count, aln_ref_names_str, aln_ref_scores_str, variant_payload, args.write_detailed_allele_table)
                alleles_list.append(allele_row)
                continue #for ambiguous reads, don't add indels to reference totals

            #iterate through payloads -- if a read aligned equally-well to two references, it could have more than one payload
            for ref_name in aln_ref_names:
                variant_payload = variantCache[variant]["variant_"+ref_name]
                if args.discard_indel_reads and (variant_payload['deletion_n'] > 0 or variant_payload['insertion_n'] > 0):
                    counts_discarded[ref_name] += variant_count
                    allele_row = get_allele_row('DISCARDED_'+aln_ref_names[0],variant_count,aln_ref_names_str,aln_ref_scores_str,variant_payload,args.write_detailed_allele_table)
                    alleles_list.append(allele_row)
                    continue

                counts_total[ref_name] += variant_count
                if variant_payload['classification'] == 'MODIFIED':
                    counts_modified[ref_name] += variant_count
                else:
                    counts_unmodified[ref_name] += variant_count

                allele_row = get_allele_row(ref_name, variant_count, aln_ref_names_str, aln_ref_scores_str, variant_payload, args.write_detailed_allele_table)
                alleles_list.append(allele_row)

                this_effective_len= refs[ref_name]['sequence_length'] #how long is this alignment (insertions increase length, deletions decrease length)

                this_has_insertions = False
                all_insertion_count_vectors[ref_name][variant_payload['all_insertion_positions']]+=variant_count
                all_insertion_left_count_vectors[ref_name][variant_payload['all_insertion_left_positions']]+=variant_count

                if not args.ignore_insertions:
                    inserted_n_dicts[ref_name][variant_payload['insertion_n']] += variant_count
                    insertion_count_vectors[ref_name][variant_payload['insertion_positions']]+=variant_count
                    this_effective_len = this_effective_len + variant_payload['insertion_n']
                    if variant_payload['insertion_n'] > 0:
                         counts_insertion[ref_name] += variant_count
                         this_has_insertions = True


                this_has_deletions = False
                all_deletion_count_vectors[ref_name][variant_payload['all_deletion_positions']]+=variant_count
                if not args.ignore_deletions:
                    deleted_n_dicts[ref_name][variant_payload['deletion_n']] += variant_count
                    deletion_count_vectors[ref_name][variant_payload['deletion_positions']]+=variant_count
                    this_effective_len = this_effective_len - variant_payload['deletion_n']
                    if variant_payload['deletion_n'] > 0:
                         counts_deletion[ref_name] += variant_count
                         this_has_deletions = True

                effective_len_dicts[ref_name][this_effective_len] += variant_count

                this_has_substitutions = False
                all_substitution_count_vectors[ref_name][variant_payload['all_substitution_positions']] += variant_count

                if not args.ignore_substitutions:
                    substituted_n_dicts[ref_name][variant_payload['substitution_n']] += variant_count
                    substitution_count_vectors[ref_name][variant_payload['substitution_positions']] += variant_count
                    if variant_payload['substitution_n'] > 0:
                         counts_substitution[ref_name] += variant_count
                         this_has_substitutions = True

                    nucs = ['A', 'T', 'C', 'G', 'N']
                    for nuc in nucs:
                        isNuc = [n == nuc for n in variant_payload['all_substitution_values']]
                        if(np.sum(isNuc) > 0):
                            locs = np.array(variant_payload['all_substitution_positions'])[isNuc]
                            all_substitution_base_vectors[ref_name + "_" + nuc ][locs] += variant_count


                if this_has_deletions:
                    if this_has_insertions:
                        if this_has_substitutions:
                            counts_insertion_and_deletion_and_substitution[ref_name] += variant_count
                        else:
                            counts_insertion_and_deletion[ref_name] += variant_count
                    else:
                        if this_has_substitutions:
                            counts_deletion_and_substitution[ref_name] += variant_count
                        else:
                            counts_only_deletion[ref_name] += variant_count
                else: #no deletions
                    if this_has_insertions:
                        if this_has_substitutions:
                            counts_insertion_and_substitution[ref_name] += variant_count
                        else:
                            counts_only_insertion[ref_name] += variant_count
                    else:
                        if this_has_substitutions:
                            counts_only_substitution[ref_name] += variant_count

                #set all_base_count_vectors
                aln_seq = variant_payload['aln_seq']
                ref_pos = variant_payload['ref_positions']
                for i in range(len(aln_seq)):
                    if ref_pos[i] < 0:
                        continue
                    nuc = aln_seq[i]
                    all_base_count_vectors[ref_name + "_" + nuc][ref_pos[i]] += variant_count

                exon_len_mods = refs[ref_name]['exon_len_mods'] #for each exon, how much length did this reference modify it?
                tot_exon_len_mod = sum(exon_len_mods) #for all exons, how much length was modified?
                if this_has_insertions or this_has_deletions or this_has_substitutions or tot_exon_len_mod != 0: #only count modified reads
                    exon_positions = refs[ref_name]['exon_positions']
                    splicing_positions = refs[ref_name]['splicing_positions']
                    insertion_coordinates = variant_payload['insertion_coordinates']
                    insertion_sizes = variant_payload['insertion_sizes']
                    all_insertion_positions = variant_payload['all_insertion_positions']
                    all_insertion_left_positions = variant_payload['all_insertion_left_positions']
                    insertion_positions = variant_payload['insertion_positions']
                    deletion_coordinates = variant_payload['deletion_coordinates']
                    deletion_sizes = variant_payload['deletion_sizes']
                    all_deletion_positions = variant_payload['all_deletion_positions']
                    deletion_positions = variant_payload['deletion_positions']
                    all_substitution_positions = variant_payload['all_substitution_positions']
                    substitution_positions = variant_payload['substitution_positions']

                    length_modified_positions_exons=[]
                    current_read_exons_modified = False
                    current_read_spliced_modified = False

                    for idx_ins, (ins_start, ins_end) in enumerate(insertion_coordinates):
                        insertion_length_vectors[ref_name][ins_start]+=(insertion_sizes[idx_ins]*variant_count)
                        insertion_length_vectors[ref_name][ins_end]+=(insertion_sizes[idx_ins]*variant_count)

                        if refs[ref_name]['contains_coding_seq']:
                            if set(exon_positions).intersection((ins_start, ins_end)): # check that we are inserting in one exon
                                current_read_exons_modified = True
                                set1 = set(exon_positions).intersection((ins_start, ins_end))
                                length_modified_positions_exons.append((insertion_sizes[idx_ins]))

                    for idx_del, (del_start, del_end) in enumerate(deletion_coordinates):
                        deletion_length_vectors[ref_name][list(range(del_start, del_end))] += (deletion_sizes[idx_del]*variant_count)

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
                            counts_splicing_sites_modified[ref_name] += variant_count

                        if tot_exon_len_mod != 0:
                            effective_length = sum(length_modified_positions_exons) + tot_exon_len_mod
                            if (effective_length % 3) == 0:
                                # indels have restored exon frame
                                counts_modified_non_frameshift[ref_name] += variant_count
                                hists_inframe[ref_name][effective_length] += variant_count
                            else:
                                counts_modified_frameshift[ref_name] += variant_count
                                hists_frameshift[ref_name][effective_length] += variant_count
                        # if modified check if frameshift
                        elif current_read_exons_modified:

                            if not length_modified_positions_exons:
                                # there are no indels
                                counts_modified_non_frameshift[ref_name] += variant_count
                                hists_inframe[ref_name][0] += variant_count
                            else:
                                effective_length = sum(length_modified_positions_exons)

                                if (effective_length % 3) == 0:
                                    counts_modified_non_frameshift[ref_name] += variant_count
                                    hists_inframe[ref_name][effective_length] += variant_count
                                else:
                                    counts_modified_frameshift[ref_name] += variant_count
                                    hists_frameshift[ref_name][effective_length] += variant_count

                        # the indels and subtitutions are outside the exon/s  so we don't care!
                        else:
                            counts_non_modified_non_frameshift[ref_name] += variant_count
                            insertion_count_vectors_noncoding[ref_name][insertion_positions] += variant_count
                            deletion_count_vectors_noncoding[ref_name][deletion_positions] += variant_count
                            substitution_count_vectors_noncoding[ref_name][substitution_positions] += variant_count
                            hists_inframe[ref_name][0] += variant_count
                #if unmodified but the tot_exon_len_mod is != 0
                elif tot_exon_len_mod != 0:
                    if (tot_exon_len_mod % 3) == 0:
                        # indels have restored exon frame
                        counts_modified_non_frameshift[ref_name] += variant_count
                        hists_inframe[ref_name][effective_length] += variant_count
                    else:
                        counts_modified_frameshift[ref_name] += variant_count
                        hists_frameshift[ref_name][effective_length] += variant_count
        #done iterating through variant cache objects



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

            all_indelsub_count_vectors[ref_name] = all_insertion_count_vectors[ref_name] + all_deletion_count_vectors[ref_name] + all_substitution_count_vectors[ref_name]
            indelsub_count_vectors[ref_name] = insertion_count_vectors[ref_name] + deletion_count_vectors[ref_name] + substitution_count_vectors[ref_name]

        # For HDR work, create a few more arrays, where all reads are aligned to ref1 (the first reference) and the indels are computed with regard to ref1
        if args.expected_hdr_amplicon_seq != "" or args.prime_editing_pegRNA_extension_seq != "":
            ref1_name = ref_names[0]
            ref1_len = refs[ref1_name]['sequence_length']

            ref1_all_insertion_count_vectors = {} #all insertions (including quantification window bases) with respect to ref1
            ref1_all_insertion_left_count_vectors = {} # 'all_insertion_left_positions' #arr with 1's to the left of where the insertion occurs
            ref1_all_deletion_count_vectors = {}
            ref1_all_substitution_count_vectors = {}
            ref1_all_indelsub_count_vectors = {}
            ref1_all_base_count_vectors = {}

            for ref_name in ref_names:
                ref1_all_insertion_count_vectors[ref_name] = np.zeros(ref1_len)
                ref1_all_insertion_left_count_vectors[ref_name] = np.zeros(ref1_len)
                ref1_all_deletion_count_vectors[ref_name] = np.zeros(ref1_len)
                ref1_all_substitution_count_vectors[ref_name] = np.zeros(ref1_len)
                ref1_all_indelsub_count_vectors[ref_name] = np.zeros(ref1_len)

                for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                    ref1_all_base_count_vectors[ref_name+"_"+nuc ] = np.zeros(ref1_len)

            #for ref1 we will add all other indels to indels that have already been found..
            ref1_all_insertion_count_vectors[ref_names[0]] = all_insertion_count_vectors[ref_names[0]].copy()
            ref1_all_insertion_left_count_vectors[ref_names[0]] = all_insertion_left_count_vectors[ref_names[0]].copy()
            ref1_all_deletion_count_vectors[ref_names[0]] = all_deletion_count_vectors[ref_names[0]].copy()
            ref1_all_substitution_count_vectors[ref_names[0]] = all_substitution_count_vectors[ref_names[0]].copy()
            ref1_all_indelsub_count_vectors[ref_names[0]] = all_indelsub_count_vectors[ref_names[0]].copy()

            for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                ref1_all_base_count_vectors[ref_names[0]+"_"+nuc ] = all_base_count_vectors[ref_names[0]+"_"+nuc].copy()

            #then go through and align and add other reads
            for variant in variantCache:
                #skip variant if there were none observed
                variant_count = variantCache[variant]['count']
                if (variant_count == 0):
                    continue

                aln_ref_names = variantCache[variant]['aln_ref_names'] #list of references this seq aligned to
                if len(aln_ref_names) == 1 and ref_names[0] == aln_ref_names[0]: #if this read was only aligned to ref1, skip it because we already included the indels when we initialized the ref1_all_deletion_count_vectors array
                    continue

                class_name = variantCache[variant]['class_name'] #for classifying read e.g. 'HDR_MODIFIED' for pie chart
                #if class is AMBIGUOUS (set above if the args.expand_ambiguous_alignments param is false) don't add the modifications in this allele to the allele summaries
                if class_name == "AMBIGUOUS":
                    continue #for ambiguous reads, don't add indels to reference totals

                #align this variant to ref1 sequence
                ref_aln_name, s1, s2, score = variantCache[variant]['ref_aln_details'][0]
                if args.use_legacy_insertion_quantification:
                    payload = CRISPRessoCOREResources.find_indels_substitutions_legacy(s1, s2, refs[ref1_name]['include_idxs'])
                else:
                    payload = CRISPRessoCOREResources.find_indels_substitutions(s1, s2, refs[ref1_name]['include_idxs'])

                #indels in this alignment against ref1 should be recorded for each ref it was originally assigned to, as well as for ref1
                #for example, if this read aligned to ref3, align this read to ref1, and add the resulting indels to ref1_all_insertion_count_vectors[ref3] as well as ref1_all_insertion_count_vectors[ref1]
                #   Thus, ref1_all_insertion_count_vectors[ref3] will show the position of indels of reads that aligned to ref3, but mapped onto ref1
                #   And ref1_alle_insertion_count_vectors[ref1] will show the position of indels of all reads, mapped onto ref1 (as copied above)
                for ref_name in aln_ref_names:
                    if ref_name == ref_names[0]:
                        continue
                    ref1_all_insertion_count_vectors[ref_name][payload['all_insertion_positions']]+=variant_count
                    ref1_all_insertion_left_count_vectors[ref_name][payload['all_insertion_left_positions']]+=variant_count
                    ref1_all_indelsub_count_vectors[ref_name][payload['all_insertion_positions']]+=variant_count

                    ref1_all_deletion_count_vectors[ref_name][payload['all_deletion_positions']]+=variant_count
                    ref1_all_indelsub_count_vectors[ref_name][payload['all_deletion_positions']]+=variant_count

                    ref1_all_substitution_count_vectors[ref_name][payload['all_substitution_positions']]+=variant_count
                    ref1_all_indelsub_count_vectors[ref_name][payload['all_substitution_positions']]+=variant_count

                    aln_seq = s1
                    ref_pos = payload['ref_positions']
                    for i in range(len(aln_seq)):
                        if ref_pos[i] < 0:
                            continue
                        nuc = aln_seq[i]
                        ref1_all_base_count_vectors[ref_name + "_" + nuc][ref_pos[i]] += variant_count

        info('Done!')

        #order class_counts
        decorated_class_counts = []
        for class_count_name in class_counts:
            thisRefInd = 100
            thisIsMod = 1
            for idx, ref_name in enumerate(ref_names):
                if class_count_name.startswith(ref_name):
                    thisRefInd = idx
                    break
            if "UNMODIFIED" in class_count_name:
                thisIsMod = 0
            decorated_class_counts.append((thisRefInd, thisIsMod, class_count_name))
        decorated_class_counts.sort()
        class_counts_order = [class_count_name for thisRefInd, thisIsMod, class_count_name in decorated_class_counts]

        if N_TOTAL == 0:
            raise CRISPRessoShared.NoReadsAlignedException('Error: No alignments were found')

        #create alleles table
        info('Calculating allele frequencies...')

        #set up allele table
        df_alleles = pd.DataFrame(alleles_list)
        #df_alleles['%Reads']=df_alleles['#Reads']/df_alleles['#Reads'].sum()*100 # sum of #reads will be >= N_TOTAL because an allele appears once for each reference it aligns to
        df_alleles['%Reads']=df_alleles['#Reads']/N_TOTAL*100
        df_alleles[['n_deleted', 'n_inserted', 'n_mutated']] = df_alleles[['n_deleted', 'n_inserted', 'n_mutated']].astype(int)

        df_alleles.sort_values(by='#Reads', ascending=False, inplace=True)

        def calculate_99_max(d):
            """
            Input is a dictionary of keys->counts
            After sorting keys, what is the largest key that contains 99% of counts
            Returns 0 if empty
            """
            curr_sum = 0
            total = sum(d.values())
            cutoff = total*.99
            for i in sorted(d.keys()):
                curr_sum += d[i]
                if curr_sum > cutoff:
                    return i
            return 0

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
                xmin, xmax=-min_cut, ref_len-max_cut
            else:
                min_cut=ref_len/2
                max_cut=ref_len/2
                xmin, xmax=-min_cut, +max_cut

            refs[ref_name]['min_cut'] = min_cut
            refs[ref_name]['max_cut'] = max_cut

            max_mut=max(15, max(substituted_n_dicts[ref_name].keys() or [0])) # the or bit is for when there are no keys
            max_ins=max(15, max(inserted_n_dicts[ref_name].keys() or [0]))
            max_del=max(15, max(deleted_n_dicts[ref_name].keys() or [0]))

            x_bins_mut = np.arange(max_mut+1)
            y_values_mut = np.array([substituted_n_dicts[ref_name][x] for x in x_bins_mut])
            x_bins_ins = np.arange(max_ins+1)
            y_values_ins = np.array([inserted_n_dicts[ref_name][x] for x in x_bins_ins])
            x_bins_del = np.arange(max_del+1)
            y_values_del = np.array([deleted_n_dicts[ref_name][x] for x in x_bins_del])

            refs[ref_name]['y_values_mut'] = y_values_mut
            refs[ref_name]['x_bins_mut'] = x_bins_mut
            refs[ref_name]['y_values_ins'] = y_values_ins
            refs[ref_name]['x_bins_ins'] = x_bins_ins
            refs[ref_name]['y_values_del'] = y_values_del
            refs[ref_name]['x_bins_del'] = x_bins_del

            min_eff_len = min(ref_len-15, min(effective_len_dicts[ref_name].keys() or [0]))
            max_eff_len = max(ref_len+15, max(effective_len_dicts[ref_name].keys() or [0]))
            hlengths = np.arange(min_eff_len, max_eff_len+1)

            hdensity = np.array([effective_len_dicts[ref_name][x] for x in hlengths])
            hlengths = hlengths-ref_len
            center_index=np.where(hlengths==0)[0][0]

            refs[ref_name]['hdensity'] = hdensity
            refs[ref_name]['hlengths'] = hlengths
            refs[ref_name]['center_index'] = center_index

            count_tot = counts_total[ref_name]
            if count_tot > 0:
                # normalize effect vectors
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
            ref_info_file = open(ref_info_file_name, 'w')
            refString = ( 'name' + "\t" +
                'sequence' + "\t" +
                'sequence_length' + "\t" +
                'min_aln_score' + "\t" +
                'gap_incentive' + "\t" +
                'sgRNA_cut_points' + "\t" +
                'sgRNA_plot_cut_points' + "\t" +
                'sgRNA_intervals' + "\t" +
                'sgRNA_plot_idxs' + "\t" +
                'sgRNA_mismatches' + "\t" +
                'sgRNA_sequences' + "\t" +
                'sgRNA_orig_sequences' + "\t" +
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
                    str(refs[ref_name]['sgRNA_plot_cut_points']) + "\t" +
                    str(refs[ref_name]['sgRNA_intervals']) + "\t" +
                    str(refs[ref_name]['sgRNA_sequences']) + "\t" +
                    str(refs[ref_name]['sgRNA_orig_sequences']) + "\t" +
                    str(refs[ref_name]['sgRNA_plot_idxs']) + "\t" +
                    str(refs[ref_name]['sgRNA_mismatches']) + "\t" +
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

        crispresso2_info['results']['ref_names'] = ref_names
        crispresso2_info['results']['refs'] = refs

        info('Saving processed data...')

        #write alleles table
        #crispresso1Cols = ["Aligned_Sequence","Reference_Sequence","NHEJ","UNMODIFIED","HDR","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        #df_alleles.loc[:,crispresso1Cols].to_csv(_jp('Alleles_frequency_table.txt'),sep='\t',header=True,index=None)
        #crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
#        crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads","Aligned_Reference_Names","Aligned_Reference_Scores"]
#        crispresso2Cols = ["Read_Sequence","Amplicon_Sequence","Amplicon_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        crispresso2Cols = ["Aligned_Sequence", "Reference_Sequence", "Reference_Name", "Read_Status", "n_deleted", "n_inserted", "n_mutated", "#Reads", "%Reads"]

        allele_frequency_table_filename = 'Alleles_frequency_table.txt'
        allele_frequency_table_fileLoc = _jp(allele_frequency_table_filename)

        allele_frequency_table_zip_filename = _jp('Alleles_frequency_table.zip')

        if args.dsODN == "":
            if args.write_detailed_allele_table:
                df_alleles.to_csv(allele_frequency_table_fileLoc, sep='\t', header=True, index=None)
            else:
                df_alleles.loc[:, crispresso2Cols].to_csv(allele_frequency_table_fileLoc, sep='\t', header=True, index=None)
        else:
            info('Writing file for alleles with dsODN')
            df_alleles["contains dsODN fw"] = df_alleles["Aligned_Sequence"].str.find(args.dsODN) > 0
            df_alleles["contains dsODN rv"] = df_alleles["Aligned_Sequence"].str.find(CRISPRessoShared.reverse_complement(args.dsODN)) > 0
            df_alleles["contains dsODN"] = df_alleles["contains dsODN fw"] | df_alleles["contains dsODN rv"]

            dsODN_cols = crispresso2Cols[:]
            dsODN_cols.append("contains dsODN")

            if len(args.dsODN) > 6:
                sub_dsODN = args.dsODN[3:-3]
                df_alleles["contains dsODN fragment fw"] = df_alleles["Aligned_Sequence"].str.find(sub_dsODN) > 0
                df_alleles["contains dsODN fragment rv"] = df_alleles["Aligned_Sequence"].str.find(CRISPRessoShared.reverse_complement(sub_dsODN)) > 0
                df_alleles["contains dsODN fragment"] = df_alleles["contains dsODN fragment fw"] | df_alleles["contains dsODN fragment rv"]
            dsODN_cols.append("contains dsODN fragment")

            if args.write_detailed_allele_table:
                df_alleles.to_csv(allele_frequency_table_fileLoc, sep='\t', header=True, index=None)
            else:
                df_alleles.loc[:, dsODN_cols].to_csv(allele_frequency_table_fileLoc, sep='\t', header=True, index=None)

        with zipfile.ZipFile(allele_frequency_table_zip_filename, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as myzip:
            myzip.write(allele_frequency_table_fileLoc, allele_frequency_table_filename)
        os.remove(allele_frequency_table_fileLoc)
        crispresso2_info['running_info']['allele_frequency_table_filename'] = os.path.basename(allele_frequency_table_filename) #filename is the name of the file in the zip
        crispresso2_info['running_info']['allele_frequency_table_zip_filename'] = os.path.basename(allele_frequency_table_zip_filename)



        if args.crispresso1_mode:
            with open(_jp('Quantification_of_editing_frequency.txt'), 'w+') as outfile:
                outfile.write("Quantification of editing frequency:\n")
                for ref_name in ref_names:
                    n_unmod = counts_unmodified[ref_name]
                    n_mod = counts_modified[ref_name]
                    n_discarded = counts_discarded[ref_name]

                    n_insertion = counts_insertion[ref_name]
                    n_deletion = counts_deletion[ref_name]
                    n_substitution = counts_substitution[ref_name]

                    outfile.write("%s: Unmodified: %d Modified: %d Discarded: %d\n" % (ref_name, n_unmod, n_mod, n_discarded))
                    outfile.write("(%d reads with insertions, %d reads with deletions, %d reads with substitutions)\n" % (n_insertion, n_deletion, n_substitution))

                outfile.write('Total Aligned:%d reads ' % N_TOTAL)

        quant_of_editing_freq_filename =_jp('CRISPResso_quantification_of_editing_frequency.txt')
        with open(quant_of_editing_freq_filename, 'w+') as outfile:
            outfile.write('Amplicon\tUnmodified%\tModified%\tReads_in_input\tReads_aligned_all_amplicons\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions\n')
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
                    unmod_pct = round(100*n_unmod/float(n_aligned), 8)
                    mod_pct = round(100*n_mod/float(n_aligned), 8)

                vals = [ref_name]
                vals.extend([str(x) for x in [unmod_pct, mod_pct, N_READS_INPUT, N_TOTAL, n_aligned, n_unmod, n_mod, n_discarded, n_insertion, n_deletion, n_substitution, n_only_insertion, n_only_deletion, n_only_substitution, n_insertion_and_deletion, n_insertion_and_substitution, n_deletion_and_substitution, n_insertion_and_deletion_and_substitution]])
                outfile.write("\t".join(vals) + "\n")

        crispresso2_info['running_info']['quant_of_editing_freq_filename'] = os.path.basename(quant_of_editing_freq_filename)


        #write statistics
        if args.crispresso1_mode:
            with open(_jp('Mapping_statistics.txt'), 'w+') as outfile:
                outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d\n' % (N_READS_INPUT, N_READS_AFTER_PREPROCESSING, N_TOTAL))

        mapping_stats_filename = _jp('CRISPResso_mapping_statistics.txt')
        with open(mapping_stats_filename, 'w+') as outfile:
            outfile.write('READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n')
            outfile.write("\t".join([str(x) for x in[N_READS_INPUT, N_READS_AFTER_PREPROCESSING, N_TOTAL, aln_stats['N_COMPUTED_ALN'], aln_stats['N_CACHED_ALN'], aln_stats['N_COMPUTED_NOTALN'], aln_stats['N_CACHED_NOTALN']]]) + "\n")
        crispresso2_info['running_info']['alignment_stats'] = aln_stats
        crispresso2_info['running_info']['mapping_stats_filename'] = os.path.basename(mapping_stats_filename)

        def save_vector_to_file(vector, filename):
            #np.savetxt(_jp('%s.txt' %name), np.vstack([(np.arange(len(vector))+1),vector]).T, fmt=['%d','%.18e'],delimiter='\t', newline='\n', header='amplicon position\teffect',footer='', comments='# ')
            np.savetxt(filename, np.vstack([(np.arange(len(vector))+1), vector]).T, fmt=['%d', '%.18e'], delimiter='\t', newline='\n', header='amplicon position\teffect', footer='', comments='# ')

        def save_count_vectors_to_file(vectors, vectorNames, refSeq, filename):
            outfile = open(filename, "w")
            outfile.write("Sequence\t"+"\t".join(list(refSeq))+"\n") #first row: reference sequence
            for vector, vectorName in zip(vectors, vectorNames):
                outfile.write(vectorName +"\t" + "\t".join([str(x) for x in vector]) + "\n") #next, vectors are printed
            outfile.close()

        crispresso2_info['results']['alignment_stats']['insertion_pct_vectors'] = insertion_pct_vectors
        crispresso2_info['results']['alignment_stats']['deletion_pct_vectors'] = insertion_pct_vectors
        crispresso2_info['results']['alignment_stats']['substitution_pct_vectors'] = insertion_pct_vectors
        crispresso2_info['results']['alignment_stats']['indelsub_pct_vectors'] = insertion_pct_vectors


        #set unique plot name to appear as prefix to files for each reference
        seen_ref_names = {}
        for ref_name in ref_names:
            #only show reference name in filenames if more than one reference
            ref_plot_name = ref_name
            if len(ref_names) == 1 and ref_names[0] == "Reference":
                ref_plot_name = ""
                seen_ref_names[ref_plot_name] = 1
                refs[ref_name]['ref_plot_name'] = ref_plot_name
                continue
            if len(ref_plot_name) > 21:
                ref_plot_name = ref_plot_name[0:21]
            #make sure it is unique
            orig_ref_plot_name = ref_plot_name
            ind = 2
            while(ref_plot_name in seen_ref_names):
                ref_plot_name = orig_ref_plot_name+"_"+ind
                ind+=1
            ref_plot_name += "."
            seen_ref_names[ref_plot_name] = 1
            refs[ref_name]['ref_plot_name'] = ref_plot_name

        for ref_name in ref_names:
            ref_plot_name = refs[ref_name]['ref_plot_name']

            #n_this_category = counts_total[ref_name]
            #if n_this_category < 1:
            #    continue

            if not args.suppress_plots:
                ins_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_insertion.txt')
                save_vector_to_file(insertion_pct_vectors[ref_name], ins_pct_vector_filename)
                crispresso2_info['results']['refs'][ref_name]['insertion_pct_vector_filename'] = os.path.basename(ins_pct_vector_filename)

                del_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_deletion.txt')
                save_vector_to_file(deletion_pct_vectors[ref_name], del_pct_vector_filename)
                crispresso2_info['results']['refs'][ref_name]['deletion_pct_vector_filename'] = os.path.basename(del_pct_vector_filename)

                sub_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_substitution.txt')
                save_vector_to_file(substitution_pct_vectors[ref_name], sub_pct_vector_filename)
                crispresso2_info['results']['refs'][ref_name]['substitution_pct_vector_filename'] = os.path.basename(sub_pct_vector_filename)

                indelsub_pct_vector_filename = _jp(ref_plot_name+'Effect_vector_combined.txt')
                save_vector_to_file(indelsub_pct_vectors[ref_name], indelsub_pct_vector_filename)
                crispresso2_info['results']['refs'][ref_name]['combined_pct_vector_filename'] = os.path.basename(indelsub_pct_vector_filename)

            #save mods in quantification window
            quant_window_mod_count_filename = _jp(ref_plot_name+'Quantification_window_modification_count_vectors.txt')
            save_count_vectors_to_file([insertion_count_vectors[ref_name],
                        deletion_count_vectors[ref_name],
                        substitution_count_vectors[ref_name],
                        indelsub_count_vectors[ref_name],
                        [counts_total[ref_name]]*refs[ref_name]['sequence_length']],
                        ['Insertions', 'Deletions', 'Substitutions', 'All_modifications', 'Total'],
                            refs[ref_name]['sequence'], quant_window_mod_count_filename)
            crispresso2_info['results']['refs'][ref_name]['quant_window_mod_count_filename'] = os.path.basename(quant_window_mod_count_filename)

            #save all mods
            mod_count_filename = _jp(ref_plot_name+'Modification_count_vectors.txt')
            save_count_vectors_to_file([all_insertion_count_vectors[ref_name],
                        all_insertion_left_count_vectors[ref_name],
                        all_deletion_count_vectors[ref_name],
                        all_substitution_count_vectors[ref_name],
                        all_indelsub_count_vectors[ref_name],
                        [counts_total[ref_name]]*refs[ref_name]['sequence_length']],
                        ['Insertions', 'Insertions_Left', 'Deletions', 'Substitutions', 'All_modifications', 'Total'],
                            refs[ref_name]['sequence'], mod_count_filename)
            crispresso2_info['results']['refs'][ref_name]['mod_count_filename'] = os.path.basename(mod_count_filename)
            crispresso2_info['results']['refs'][ref_name]['mod_count_filename_caption'] = "A tab-separated file showing the number of modifications for each position in the amplicon. " \
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
                with open(frameshift_analysis_filename, 'w+') as outfile:
                        outfile.write('Frameshift analysis:\n\tNoncoding mutation:%d reads\n\tIn-frame mutation:%d reads\n\tFrameshift mutation:%d reads\n' %(NON_MODIFIED_NON_FRAMESHIFT, MODIFIED_NON_FRAMESHIFT, MODIFIED_FRAMESHIFT))
                crispresso2_info['results']['refs'][ref_name]['frameshift_analysis_filename'] = os.path.basename(frameshift_analysis_filename)
                crispresso2_info['results']['refs'][ref_name]['frameshift_analysis_filename_caption'] = "A text file describing the number of noncoding, in-frame, and frameshift mutations. This report file is produced when the amplicon contains a coding sequence."

                splice_sites_analysis_filename = _jp(ref_plot_name+'Splice_sites_analysis.txt')
                with open(splice_sites_analysis_filename, 'w+') as outfile:
                        outfile.write('Splice sites analysis:\n\tUnmodified:%d reads\n\tPotential splice sites modified:%d reads\n' %(counts_total[ref_name]- SPLICING_SITES_MODIFIED, SPLICING_SITES_MODIFIED))
                crispresso2_info['results']['refs'][ref_name]['splice_sites_analysis_filename'] = os.path.basename(splice_sites_analysis_filename)
                crispresso2_info['results']['refs'][ref_name]['splice_sites_analysis_filename_caption'] = "A text file describing the number of splicing sites that are unmodified and modified. This file report is produced when the amplicon contains a coding sequence."

                ins_pct_vector_noncoding_filename = _jp(ref_plot_name+'Effect_vector_insertion_noncoding.txt')
                save_vector_to_file(insertion_pct_vectors_noncoding[ref_name], ins_pct_vector_noncoding_filename)
                crispresso2_info['results']['refs'][ref_name]['insertion_pct_vector_noncoding_filename'] = os.path.basename(ins_pct_vector_noncoding_filename)
                crispresso2_info['results']['refs'][ref_name]['insertion_pct_vector_noncoding_filename_caption'] = "A tab-separated text file with a one-row header that shows the percentage of reads with a noncoding insertion at each base in the " + ref_name + " sequence. " \
                    "The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding insertion at that location. This report file is produced when the amplicon contains a coding sequence."

                del_pct_vector_noncoding_filename = _jp(ref_plot_name+'Effect_vector_deletion_noncoding.txt')
                save_vector_to_file(deletion_pct_vectors_noncoding[ref_name], del_pct_vector_noncoding_filename)
                crispresso2_info['results']['refs'][ref_name]['deletion_pct_vector_noncoding_filename'] = os.path.basename(del_pct_vector_noncoding_filename)
                crispresso2_info['results']['refs'][ref_name]['deletion_pct_vector_noncoding_filename_caption'] = "A tab-separated text file with a one-row header that shows the percentage of reads with a noncoding deletion at each base in the " + ref_name + " sequence. " \
                    "The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding deletion at that location. This report file is produced when the amplicon contains a coding sequence."

                sub_pct_vector_noncoding_filename = _jp(ref_plot_name+'Effect_vector_substitution_noncoding.txt')
                save_vector_to_file(substitution_pct_vectors_noncoding[ref_name], sub_pct_vector_noncoding_filename)
                crispresso2_info['results']['refs'][ref_name]['substitution_pct_vector_noncoding_filename'] = os.path.basename(sub_pct_vector_noncoding_filename)
                crispresso2_info['results']['refs'][ref_name]['substitution_pct_vector_noncoding_filename_caption'] = "A tab-separated text file with a one-row header that shows the percentage of reads with a noncoding substitution at each base in the " + ref_name + " sequence. " \
                    "The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a nondcoding substitution at that location. This report file is produced when the amplicon contains a coding sequence."

            if args.dump:
                if refs[ref_name]['sgRNA_cut_points']:
                    with open(_jp(ref_plot_name + 'Cut_points.json'), 'w') as fh:
                        json.dump(refs[ref_name]['sgRNA_cut_points'], fh)

                if refs[ref_name]['sgRNA_intervals']:
                    with open(_jp(ref_plot_name + 'sgRNA_intervals.json'), 'w') as fh:
                        json.dump(refs[ref_name]['sgRNA_intervals'], fh)

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
                pd.DataFrame({'indel_size':pd.Series(hlengths, dtype='int'),'fq':hdensity})[['indel_size', 'fq']].to_csv(indel_histogram_file, index=None, sep='\t')
                crispresso2_info['results']['refs'][ref_name]['indel_histogram_filename'] = os.path.basename(indel_histogram_file)
                crispresso2_info['results']['refs'][ref_name]['indel_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the length of indels (both insertions and deletions) in the " + ref_name +" sequence in the quantification window. " \
                    "Indels outside of the quantification window are not included. The indel_size column shows the number of substitutions, and the fq column shows the number of reads having an indel of that length."


                insertion_histogram_file = _jp(ref_plot_name+'Insertion_histogram.txt')
                pd.DataFrame({'ins_size':x_bins_ins,'fq':y_values_ins})[["ins_size", "fq"]].to_csv(insertion_histogram_file, index=None, sep='\t')
                crispresso2_info['results']['refs'][ref_name]['insertion_histogram_filename'] = os.path.basename(insertion_histogram_file)
                crispresso2_info['results']['refs'][ref_name]['insertion_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of insertions in the " + ref_name +" sequence in the quantification window. " \
                    "Insertions outside of the quantification window are not included. The ins_size column shows the number of insertions, and the fq column shows the number of reads having that number of insertions."


                deletion_histogram_file = _jp(ref_plot_name+'Deletion_histogram.txt')
                pd.DataFrame({'del_size':-1*x_bins_del,'fq':y_values_del})[['del_size', 'fq']].to_csv(deletion_histogram_file, index=None, sep='\t')
                crispresso2_info['results']['refs'][ref_name]['deletion_histogram_filename'] = os.path.basename(deletion_histogram_file)
                crispresso2_info['results']['refs'][ref_name]['deletion_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of deletions in the " + ref_name +" sequence in the quantification window. " \
                    "Deletions outside of the quantification window are not included. The del_size column shows the number of deletions, and the fq column shows the number of reads having that number of deletions."


                substitution_histogram_file = _jp(ref_plot_name+'Substitution_histogram.txt')
                pd.DataFrame({'sub_count':x_bins_mut,'fq':y_values_mut})[['sub_count', 'fq']].to_csv(substitution_histogram_file, index=None, sep='\t')
                crispresso2_info['results']['refs'][ref_name]['substitution_histogram_filename'] = os.path.basename(substitution_histogram_file)
                crispresso2_info['results']['refs'][ref_name]['substitution_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of substitutions in the " + ref_name +" sequence in the quantification window. " \
                    "Substitutions outside of the quantification window are not included. The sub_size column shows the number of substitutions, and the fq column shows the number of reads having that number of substitutions."

            if args.dump:
                np.savez(_jp(ref_plot_name+'Effect_vector_insertion'), insertion_pct_vectors[ref_name])
                np.savez(_jp(ref_plot_name+'Effect_vector_deletion'), deletion_pct_vectors[ref_name])
                np.savez(_jp(ref_plot_name+'Effect_vector_substitution'), substitution_pct_vectors[ref_name])

                np.savez(_jp(ref_plot_name+'Effect_vector_combined'), indelsub_pct_vectors[ref_name])

                np.savez(_jp(ref_plot_name+'Position_dependent_vector_avg_insertion_size'), insertion_length_vectors[ref_name])
                np.savez(_jp(ref_plot_name+'Position_dependent_vector_avg_deletion_size'), deletion_length_vectors[ref_name])

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
        def get_plot_title_with_ref_name(plotTitle, ref_name):
            if n_refs > 1:
                return (plotTitle + ": " + ref_name)
            return plotTitle

        def count_alternate_alleles(sub_base_vectors, ref_name, ref_sequence, ref_total_aln_reads):
            #create vectors with all allele frequencies -- not just the substitution (the reference allele will not be 0)
            alph = ['A', 'C', 'G', 'T', 'N']

            #count the total number of times each substitution occurs
            count_sub_base_vectors = {}
            alt_nuc_counts = {}
            for a in alph:
                alt_nuc_counts[a] = {}
                count_sub_base_vectors[a] = list(sub_base_vectors[ref_name+"_"+a])
                for b in alph:
                    alt_nuc_counts[a][b] = 0

            for idx, c in enumerate(ref_sequence):
                tot_sub_at_idx = 0
                for a in alph:
                    sub = sub_base_vectors[ref_name+"_" + a][idx]
                    alt_nuc_counts[c][a] += sub
                    tot_sub_at_idx += sub

            #df_subs = pd.DataFrame([count_sub_base_vectors["A"],count_sub_base_vectors["C"],count_sub_base_vectors["G"],count_sub_base_vectors["T"],count_sub_base_vectors["N"]])
            df_subs = pd.DataFrame([count_sub_base_vectors[a] for a in alph])
            df_subs.index = alph
            df_subs.columns = list(ref_sequence)
            return (df_subs, alt_nuc_counts)

            ############

        plot_pool = ProcessPoolExecutor(n_processes)
        plot_results = []
        ###############################################################################################################################################
        ### FIGURE 1: Alignment
        if not args.suppress_plots:
            plot_1a_root = _jp("1a.Read_barplot")
            plot_1a_input = {
                'N_READS_INPUT': N_READS_INPUT,
                'N_READS_AFTER_PREPROCESSING': N_READS_AFTER_PREPROCESSING,
                'N_TOTAL': N_TOTAL,
                'plot_root': plot_1a_root,
                'save_png': save_png
            }
            if n_processes > 1:
                plot_results.append(plot_pool.submit(
                    CRISPRessoPlot.plot_read_barplot,
                    **plot_1a_input,
                ))
            else:
                CRISPRessoPlot.plot_read_barplot(**plot_1a_input)

            crispresso2_info['results']['general_plots']['plot_1a_root'] = os.path.basename(plot_1a_root)
            crispresso2_info['results']['general_plots']['plot_1a_caption'] = "Figure 1a: The number of reads in input fastqs, after preprocessing, and after alignment to amplicons."
            crispresso2_info['results']['general_plots']['plot_1a_data'] = [('Mapping statistics', os.path.basename(mapping_stats_filename))]

            plot_1b_root = _jp("1b.Alignment_pie_chart")
            plot_1c_root = _jp('1c.Alignment_barplot')
            plot_1bc_input = {
                'class_counts_order': class_counts_order,
                'class_counts': class_counts,
                'ref_names': ref_names,
                'expected_hdr_amplicon_seq': args.expected_hdr_amplicon_seq,
                'N_TOTAL': N_TOTAL,
                'piechart_plot_root': plot_1b_root,
                'barplot_plot_root': plot_1c_root,
                'save_png': save_png
            }
            crispresso2_info['results']['general_plots']['plot_1b_root'] = os.path.basename(plot_1b_root)
            crispresso2_info['results']['general_plots']['plot_1b_caption'] = "Figure 1b: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles."
            if args.expected_hdr_amplicon_seq != "":
                crispresso2_info['results']['general_plots']['plot_1b_caption'] = "Figure 1b: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles. NHEJ reads align more closely to the unmodified reference sequence, but have mutations present in the specified quantification window. HDR reads align to the HDR reference sequence and have no mutations in the specified quantification window. Imperfect HDR reads have mutations in the specified window. AMBIGUOUS reads align equally well to the unmodified and HDR reference sequences."
            crispresso2_info['results']['general_plots']['plot_1b_data'] = [('Quantification of editing', os.path.basename(quant_of_editing_freq_filename))]

            crispresso2_info['results']['general_plots']['plot_1c_root'] = os.path.basename(plot_1c_root)
            crispresso2_info['results']['general_plots']['plot_1c_caption'] = "Figure 1c: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles."
            crispresso2_info['results']['general_plots']['plot_1c_data'] = [('Quantification of editing', os.path.basename(quant_of_editing_freq_filename))]
            if n_processes > 1:
                plot_results.append(plot_pool.submit(
                    CRISPRessoPlot.plot_class_piechart_and_barplot,
                    **plot_1bc_input,
                ))
            else:
                CRISPRessoPlot.plot_class_piechart_and_barplot(
                    **plot_1bc_input,
                )
            # to test, run: plot_pool.apply_async(CRISPRessoPlot.plot_class_piechart_and_barplot, kwds=plot_1bc_input).get()


            #1d for dsODN
            if args.dsODN != "":
                #(1b) a piechart of classes
                n_contain_dsODN = df_alleles[df_alleles['contains dsODN'] == True]['#Reads'].sum()
                n_not_contain_dsODN = df_alleles[df_alleles['contains dsODN'] == False]['#Reads'].sum()
                labels = ['Contains dsODN\n('+str(n_contain_dsODN)+" reads)",
                    'No dsODN\n('+str(n_not_contain_dsODN) + ' reads)']
                sizes = [100*n_contain_dsODN/float(N_TOTAL),
                    100*n_not_contain_dsODN/float(N_TOTAL)]
                plot_root = _jp("1d.Detection_of_dsODN")
                plot_1d_input = {
                    'sizes': sizes,
                    'labels': labels,
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_class_dsODN_piechart,
                        **plot_1d_input,
                    ))
                else:
                    CRISPRessoPlot.plot_class_dsODN_piechart(**plot_1d_input)

                crispresso2_info['results']['general_plots']['plot_1d_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['general_plots']['plot_1d_caption'] = "Figure 1d: Frequency of detection of dsODN " + args.dsODN
                crispresso2_info['results']['general_plots']['plot_1d_data'] = [('Allele table', os.path.basename(allele_frequency_table_filename))]
        ###############################################################################################################################################

        #hold mod pct dfs for each amplicon/guide
        mod_pct_dfs = {}

        for ref_name in ref_names:
            ref_len = refs[ref_name]['sequence_length']
            ref_seq = refs[ref_name]['sequence']
            min_cut = refs[ref_name]['min_cut']
            max_cut = refs[ref_name]['max_cut']
            hdensity = refs[ref_name]['hdensity']
            hlengths = refs[ref_name]['hlengths']
            center_index = refs[ref_name]['center_index']
            include_idxs_list = refs[ref_name]['include_idxs']
            quantification_window_ref_seq = [list(ref_seq)[x] for x in include_idxs_list]
            sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
            sgRNA_orig_sequences = refs[ref_name]['sgRNA_orig_sequences']
            cut_points = refs[ref_name]['sgRNA_cut_points']
            plot_cut_points = refs[ref_name]['sgRNA_plot_cut_points']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
            sgRNA_names = refs[ref_name]['sgRNA_names']
            sgRNA_mismatches = refs[ref_name]['sgRNA_mismatches']
            tot_aln_reads = counts_total[ref_name]
            n_this_category = counts_total[ref_name]

            #only show reference name in filenames if more than one reference
            ref_plot_name = refs[ref_name]['ref_plot_name']

            if n_this_category < 1:
                continue

            #plot quilt for this amplicon  (if not crispresso1 mode)
            if not args.crispresso1_mode:
                ##nucleotide counts
                df_nuc_freq = pd.DataFrame([base_count_vectors[ref_name+"_A"], base_count_vectors[ref_name+"_C"], base_count_vectors[ref_name+"_G"], base_count_vectors[ref_name+"_T"], base_count_vectors[ref_name+"_N"], base_count_vectors[ref_name+'_-']])
                df_nuc_freq.index = ['A', 'C', 'G', 'T', 'N', '-']
                df_nuc_freq.columns = quantification_window_ref_seq
                #print table showing nuc frequencies (sum to total alleles) (in quantification window)
                quant_window_nuc_freq_filename = _jp(ref_plot_name + 'Quantification_window_nucleotide_frequency_table.txt')
                df_nuc_freq.to_csv(quant_window_nuc_freq_filename, sep='\t', header=True, index=True)
                crispresso2_info['results']['refs'][ref_name]['quant_window_nuc_freq_filename'] = os.path.basename(quant_window_nuc_freq_filename)

                df_nuc_pct = df_nuc_freq.divide(tot_aln_reads)
                quant_window_nuc_pct_filename = _jp(ref_plot_name + 'Quantification_window_nucleotide_percentage_table.txt')
                df_nuc_pct.to_csv(quant_window_nuc_pct_filename, sep='\t', header=True, index=True)
                crispresso2_info['results']['refs'][ref_name]['quant_window_nuc_pct_filename'] = os.path.basename(quant_window_nuc_pct_filename)

                df_nuc_freq_all = pd.DataFrame([all_base_count_vectors[ref_name+"_A"], all_base_count_vectors[ref_name+"_C"], all_base_count_vectors[ref_name+"_G"], all_base_count_vectors[ref_name+"_T"], all_base_count_vectors[ref_name+"_N"], all_base_count_vectors[ref_name+'_-']])
                df_nuc_freq_all.index = ['A', 'C', 'G', 'T', 'N', '-']
                df_nuc_freq_all.columns = list(ref_seq)
                #print table showing nuc frequencies (sum to total alleles) (in entire region)
                nuc_freq_filename = _jp(ref_plot_name + 'Nucleotide_frequency_table.txt')
                df_nuc_freq_all.to_csv(nuc_freq_filename, sep='\t', header=True, index=True)
                crispresso2_info['results']['refs'][ref_name]['nuc_freq_filename'] = os.path.basename(nuc_freq_filename)

                df_nuc_pct_all = df_nuc_freq_all.divide(tot_aln_reads)
                nuc_pct_filename = _jp(ref_plot_name + 'Nucleotide_percentage_table.txt')
                df_nuc_pct_all.to_csv(nuc_pct_filename, sep='\t', header=True, index=True)
                crispresso2_info['results']['refs'][ref_name]['nuc_pct_filename'] = os.path.basename(nuc_pct_filename)

                if args.base_editor_output:
                    #substitution frequencies
                    df_sub_freq, alt_nuc_counts = count_alternate_alleles(
                        sub_base_vectors = substitution_base_vectors,
                        ref_name = ref_name,
                        ref_sequence = quantification_window_ref_seq,
                        ref_total_aln_reads = tot_aln_reads
                        )

                    #print table showing sub frequencies
                    quant_window_sub_freq_filename =_jp(ref_plot_name + 'Quantification_window_substitution_frequency_table.txt')
                    df_sub_freq.to_csv(quant_window_sub_freq_filename, sep='\t', header=True, index=True)
                    crispresso2_info['results']['refs'][ref_name]['quant_window_sub_freq_filename'] = os.path.basename(quant_window_sub_freq_filename)


                    df_sub_freq_all, alt_nuc_counts_all = count_alternate_alleles(
                        sub_base_vectors = all_substitution_base_vectors,
                        ref_name = ref_name,
                        ref_sequence = ref_seq,
                        ref_total_aln_reads = tot_aln_reads
                        )

                    sub_freq_table_filename = _jp(ref_plot_name + 'Substitution_frequency_table.txt')
                    df_sub_freq_all.to_csv(sub_freq_table_filename, sep='\t', header=True, index=True)
                    crispresso2_info['results']['refs'][ref_name]['sub_freq_table_filename'] = os.path.basename(sub_freq_table_filename)

                if not args.suppress_plots:
                    mod_pcts = []
                    tot = float(counts_total[ref_name])
                    mod_pcts.append(np.concatenate((['Insertions'], np.array(all_insertion_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Insertions_Left'], np.array(all_insertion_left_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Deletions'], np.array(all_deletion_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Substitutions'], np.array(all_substitution_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['All_modifications'], np.array(all_indelsub_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate((['Total'], [counts_total[ref_name]]*refs[ref_name]['sequence_length'])))
                    colnames = ['Modification']+list(ref_seq)
                    modification_percentage_summary_df = pd.DataFrame(mod_pcts, columns=colnames).apply(pd.to_numeric, errors='ignore')

                    nuc_df_for_plot = df_nuc_pct_all.reset_index().rename(columns={'index':'Nucleotide'})
                    nuc_df_for_plot.insert(0, 'Batch', ref_name) #this function was designed for plottin batch... so just add a column in there to make it happy
                    mod_df_for_plot = modification_percentage_summary_df.copy()
                    mod_df_for_plot.insert(0, 'Batch', ref_name)

                    plot_root = _jp('2a.'+ref_plot_name + 'Nucleotide_percentage_quilt')
                    plot_2a_input = {
                        'nuc_pct_df': nuc_df_for_plot,
                        'mod_pct_df': mod_df_for_plot,
                        'fig_filename_root': plot_root,
                        'save_also_png': save_png,
                        'sgRNA_intervals': sgRNA_intervals,
                        'sgRNA_names': sgRNA_names,
                        'sgRNA_mismatches': sgRNA_mismatches,
                        'quantification_window_idxs': include_idxs_list,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_nucleotide_quilt,
                            **plot_2a_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_nucleotide_quilt(**plot_2a_input)
                    crispresso2_info['results']['refs'][ref_name]['plot_2a_root'] = os.path.basename(plot_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_2a_caption'] = "Figure 2a: Nucleotide distribution across amplicon. At each base in the reference amplicon, the percentage of each base as observed in sequencing reads is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the percentage of reads for which that base was deleted. Brown bars between bases show the percentage of reads having an insertion at that position."
                    crispresso2_info['results']['refs'][ref_name]['plot_2a_data'] = [('Nucleotide frequency table', os.path.basename(nuc_freq_filename))]

                    crispresso2_info['results']['refs'][ref_name]['plot_2b_roots'] = []
                    crispresso2_info['results']['refs'][ref_name]['plot_2b_captions'] = []
                    crispresso2_info['results']['refs'][ref_name]['plot_2b_datas'] = []
                    for i in range(len(cut_points)):
                        cut_point = cut_points[i]
                        sgRNA = sgRNA_orig_sequences[i]
                        sgRNA_name = sgRNA_names[i]

                        sgRNA_label = "sgRNA_"+sgRNA # for file names
                        sgRNA_legend = "sgRNA " + sgRNA # for legends
                        if sgRNA_name != "":
                            sgRNA_label = sgRNA_name
                            sgRNA_legend = sgRNA_name + " (" + sgRNA +")"
                        sgRNA_label = CRISPRessoShared.slugify(sgRNA_label)

                        #get nucleotide columns to print for this sgRNA
                        sel_cols = [0, 1]
                        plot_half_window = max(1, args.plot_window_size)
                        new_sel_cols_start = max(2, cut_point-plot_half_window+1)
                        new_sel_cols_end = min(ref_len, cut_point+plot_half_window+1)
                        sel_cols.extend(list(range(new_sel_cols_start+2, new_sel_cols_end+2)))
                        #get new intervals
                        new_sgRNA_intervals = []
                        #add annotations for each sgRNA (to be plotted on this sgRNA's plot)
                        for (int_start, int_end) in refs[ref_name]['sgRNA_intervals']:
                            new_sgRNA_intervals += [(int_start - new_sel_cols_start, int_end - new_sel_cols_start)]
                        new_include_idx = []
                        for x in include_idxs_list:
                            new_include_idx += [x - new_sel_cols_start]
                        plot_root = _jp('2b.'+ref_plot_name + 'Nucleotide_percentage_quilt_around_' + sgRNA_label)
                        plot_2b_input = {
                            'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
                            'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
                            'fig_filename_root': plot_root,
                            'save_also_png': save_png,
                            'sgRNA_intervals': new_sgRNA_intervals,
                            'sgRNA_names': sgRNA_names,
                            'sgRNA_mismatches': sgRNA_mismatches,
                            'quantification_window_idxs': new_include_idx,
                        }
                        if n_processes > 1:
                            plot_results.append(plot_pool.submit(
                                CRISPRessoPlot.plot_nucleotide_quilt, **plot_2b_input,
                            ))
                        else:
                            CRISPRessoPlot.plot_nucleotide_quilt(
                                **plot_2b_input,
                            )
                        crispresso2_info['results']['refs'][ref_name]['plot_2b_roots'].append(os.path.basename(plot_root))
                        crispresso2_info['results']['refs'][ref_name]['plot_2b_captions'].append('Figure 2b: Nucleotide distribution around the ' + sgRNA_legend + '.')
                        crispresso2_info['results']['refs'][ref_name]['plot_2b_datas'].append([('Nucleotide frequency in quantification window', os.path.basename(quant_window_nuc_freq_filename))])


            ###############################################################################################################################################
            #(3)visualize effective lengths of reads aligning to this amplicon

            if not args.suppress_plots:

                if counts_total[ref_name] < 1:
                    continue

                xmin = min(hlengths)
                xmax = max(hlengths)
                #get 99% cutoff
                if not args.plot_histogram_outliers:
                    sum_cutoff = .99 * hdensity.sum()
                    sum_so_far = 0
                    for idx, val in enumerate(hlengths):
                        sum_so_far += hdensity[idx]
                        if sum_so_far > sum_cutoff:
                            xmax = val
                            break
                    sum_so_far = 0
                    for idx, val in enumerate(hlengths[::-1]):
                        sum_so_far += hdensity[idx]
                        if sum_so_far > sum_cutoff:
                            xmin = val
                            break
                xmin = min(xmin, -15)
                xmax = max(xmax, 15)

                plot_root = _jp(
                    '3a.' + ref_plot_name + 'Indel_size_distribution',
                )
                plot_3a_input = {
                    'hdensity': hdensity,
                    'hlengths': hlengths,
                    'center_index': center_index,
                    'n_this_category': counts_total[ref_name],
                    'xmin': xmin,
                    'xmax': xmax,
                    'title': get_plot_title_with_ref_name(
                        'Indel size distribution', ref_name,
                    ),
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_indel_size_distribution,
                        **plot_3a_input,
                    ))
                else:
                    CRISPRessoPlot.plot_indel_size_distribution(**plot_3a_input)
                clipped_string = ""
                if xmax < max(hlengths):
                    clipped_string += " (Maximum " + str(int(max(hlengths))) + " not shown)"
                if xmin > min(hlengths):
                    clipped_string += " (Minimum " + str(int(min(hlengths))) + " not shown)"
                if clipped_string != "":
                    clipped_string = " Note that histograms are clipped to show 99% of the data. To show all data, run using the parameter '--plot_histogram_outliers'. " + clipped_string

                crispresso2_info['results']['refs'][ref_name]['plot_3a_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_name]['plot_3a_caption'] = "Figure 3a: Frequency distribution of alleles with indels (blue) and without indels (red)." + clipped_string
                crispresso2_info['results']['refs'][ref_name]['plot_3a_data'] = [('Indel histogram', os.path.basename(crispresso2_info['results']['refs'][ref_name]['indel_histogram_filename']))]
                ###############################################################################################################################################

                ###############################################################################################################################################
                #(3b) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);

                xmax_ins = max(x_bins_ins)
                if not args.plot_histogram_outliers:
                    sum_cutoff = 0.99 * hdensity.sum()
                    sum_so_far = 0
                    for idx, val in enumerate(y_values_ins):
                        sum_so_far += x_bins_ins[idx]
                        if sum_so_far > sum_cutoff:
                            xmax_ins = val
                            break
                xmax_ins = max(15, xmax_ins)

                clipped_string = ""
                if xmax_ins < max(x_bins_ins):
                    clipped_string += " (Insertion maximum " + str(int(max(x_bins_ins))) + " not shown)"
                xmax_del = max(x_bins_del)
                if not args.plot_histogram_outliers:
                    sum_cutoff = .99 * hdensity.sum()
                    sum_so_far = 0
                    for idx, val in enumerate(y_values_del):
                        sum_so_far += x_bins_del[idx]
                        if sum_so_far > sum_cutoff:
                            xmax_del = val
                            break
                xmax_del = max(15, xmax_del)

                if xmax_del < max(x_bins_del):
                    clipped_string += " (Deletion minimum -" + str(int(max(x_bins_del))) + " not shown)"
                xmax_mut = max(x_bins_mut)
                if not args.plot_histogram_outliers:
                    sum_cutoff = .99 * hdensity.sum()
                    sum_so_far = 0
                    for idx, val in enumerate(y_values_mut):
                        sum_so_far += x_bins_mut[idx]
                        if sum_so_far > sum_cutoff:
                            xmax_mut = val
                            break
                xmax_mut = max(15, xmax_mut)
                if xmax_mut < max(x_bins_mut):
                    clipped_string += " (Mutation maximum " + str(int(max(x_bins_mut))) + " not shown)"

                plot_root =  _jp(
                    '3b.' + ref_plot_name + 'Insertion_deletion_substitutions_size_hist',
                )
                plot_3b_input = {
                    'ref': refs[ref_name],
                    'counts_total': counts_total[ref_name],
                    'plot_path': plot_root,
                    'plot_titles': {
                        'ins': get_plot_title_with_ref_name(
                            'Insertions', ref_name,
                        ),
                        'del': get_plot_title_with_ref_name(
                            'Deletions', ref_name,
                        ),
                        'mut': get_plot_title_with_ref_name(
                            'Substitutions', ref_name,
                        ),
                    },
                    'xmax_del': xmax_del,
                    'xmax_ins': xmax_ins,
                    'xmax_mut': xmax_mut,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_frequency_deletions_insertions,
                        **plot_3b_input,
                    ))
                else:
                    CRISPRessoPlot.plot_frequency_deletions_insertions(
                        **plot_3b_input,
                    )

                if clipped_string != "":
                    clipped_string = " Note that histograms are clipped to show 99% of the data. To show all data, run using the parameter '--plot_histogram_outliers'. " + clipped_string

                crispresso2_info['results']['refs'][ref_name]['plot_3b_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_name]['plot_3b_caption'] = "Figure 3b: Left panel, frequency distribution of sequence modifications that increase read length with respect to the reference amplicon, classified as insertions (positive indel size). Middle panel, frequency distribution of sequence modifications that reduce read length with respect to the reference amplicon, classified as deletions (negative indel size). Right panel, frequency distribution of sequence modifications that do not alter read length with respect to the reference amplicon, which are classified as substitutions (number of substituted positions shown)." + clipped_string
                crispresso2_info['results']['refs'][ref_name]['plot_3b_data'] = [('Insertions frequency', crispresso2_info['results']['refs'][ref_name]['insertion_histogram_filename']), ('Deletions Frequency', crispresso2_info['results']['refs'][ref_name]['deletion_histogram_filename']), ('Substitutions Frequency', crispresso2_info['results']['refs'][ref_name]['substitution_histogram_filename'])]


                #(4) another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);
                #Indels location Plots

                n_this_category_modified = 0
                modifiedName = ref_name + "_MODIFIED"
                if modifiedName in class_counts:
                    n_this_category_modified = class_counts[modifiedName]

                y_max = max(all_indelsub_count_vectors[ref_name]) * 1.1
                plot_root = _jp(
                    '4a.' + ref_plot_name + 'Combined_insertion_deletion_substitution_locations',
                )
                plot_4a_input = {
                    'all_indelsub_count_vectors': all_indelsub_count_vectors[ref_name],
                    'include_idxs_list': include_idxs_list,
                    'cut_points': cut_points,
                    'plot_cut_points': plot_cut_points,
                    'sgRNA_intervals': sgRNA_intervals,
                    'n_total': N_TOTAL,
                    'n_this_category': n_this_category,
                    'ref_name': ref_name,
                    'num_refs': len(ref_names),
                    'ref_len': ref_len,
                    'y_max': y_max,
                    'plot_titles': {
                        'combined': get_plot_title_with_ref_name(
                            'Combined Insertions/Deletions/Substitutions', ref_name,
                        ),
                        'main': get_plot_title_with_ref_name(
                            'Mutation position distribution', ref_name,
                        ),
                    },
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_amplicon_modifications, **plot_4a_input,
                    ))
                else:
                    CRISPRessoPlot.plot_amplicon_modifications(**plot_4a_input)
                crispresso2_info['results']['refs'][ref_name]['plot_4a_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_name]['plot_4a_caption'] = "Figure 4a: Combined frequency of any modification across the amplicon. Modifications outside of the quantification window are also shown."
                crispresso2_info['results']['refs'][ref_name]['plot_4a_data'] = []

                plot_root = _jp(
                    '4b.' + ref_plot_name + 'Insertion_deletion_substitution_locations',
                )
                plot_4b_input = {
                    'include_idxs_list': include_idxs_list,
                    'all_insertion_count_vectors': all_insertion_count_vectors[ref_name],
                    'all_deletion_count_vectors': all_deletion_count_vectors[ref_name],
                    'all_substitution_count_vectors': all_substitution_count_vectors[ref_name],
                    'sgRNA_intervals': sgRNA_intervals,
                    'ref_len': ref_len,
                    'ref_name': ref_name,
                    'num_refs': len(ref_names),
                    'n_total': N_TOTAL,
                    'n_this_category': n_this_category,
                    'cut_points': cut_points,
                    'plot_cut_points': plot_cut_points,
                    'y_max': y_max,
                    'plot_title': get_plot_title_with_ref_name(
                        'Mutation position distribution', ref_name,
                    ),
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_modification_frequency,
                        **plot_4b_input,
                    ))
                else:
                    CRISPRessoPlot.plot_modification_frequency(**plot_4b_input)
                crispresso2_info['results']['refs'][ref_name]['plot_4b_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_name]['plot_4b_caption'] = "Figure 4b: Frequency of insertions (red), deletions (purple), and substitutions (green) across the entire amplicon, including modifications outside of the quantification window."
                crispresso2_info['results']['refs'][ref_name]['plot_4b_data'] = [('Modification frequency', os.path.basename(mod_count_filename))]

                plot_root = _jp(
                    '4c.' + ref_plot_name + 'Quantification_window_insertion_deletion_substitution_locations',
                )
                plot_4c_input = {
                    'insertion_count_vectors': insertion_count_vectors[ref_name],
                    'deletion_count_vectors': deletion_count_vectors[ref_name],
                    'substitution_count_vectors': substitution_count_vectors[ref_name],
                    'include_idxs_list': include_idxs_list,
                    'cut_points': cut_points,
                    'plot_cut_points': plot_cut_points,
                    'sgRNA_intervals': sgRNA_intervals,
                    'ref_len': ref_len,
                    'num_refs': len(ref_names),
                    'n_total': N_TOTAL,
                    'n_this_category': n_this_category,
                    'plot_title': get_plot_title_with_ref_name(
                        'Mutation position distribution', ref_name,
                    ),
                    'ref_name': ref_name,
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_quantification_window_locations,
                        **plot_4c_input,
                    ))
                else:
                    CRISPRessoPlot.plot_quantification_window_locations(
                        **plot_4c_input,
                    )
                crispresso2_info['results']['refs'][ref_name]['plot_4c_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_name]['plot_4c_caption'] = "Figure 4c: Frequency of insertions (red), deletions (purple), and substitutions (green) across the entire amplicon, considering only modifications that overlap with the quantification window."
                crispresso2_info['results']['refs'][ref_name]['plot_4c_data'] = [('Modification frequency in quantification window', os.path.basename(quant_window_mod_count_filename))]

                #Position dependent indels plot
                plot_root = _jp(
                    '4d.' + ref_plot_name + 'Position_dependent_average_indel_size',
                )
                plot_4d_input = {
                    'insertion_length_vectors': insertion_length_vectors[ref_name],
                    'deletion_length_vectors': deletion_length_vectors[ref_name],
                    'cut_points': cut_points,
                    'plot_cut_points': plot_cut_points,
                    'ref_len': ref_len,
                    'plot_titles': {
                        'ins': get_plot_title_with_ref_name(
                            'Position dependent insertion size', ref_name,
                        ),
                        'del': get_plot_title_with_ref_name(
                            'Position dependent deletion size', ref_name,
                        ),
                    },
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_position_dependent_indels,
                        **plot_4d_input,
                    ))
                else:
                    CRISPRessoPlot.plot_position_dependent_indels(
                        **plot_4d_input,
                    )
                crispresso2_info['results']['refs'][ref_name]['plot_4d_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_name]['plot_4d_caption'] = "Figure 4d: Position dependent insertion size(left) and deletion size (right), including only modifications that overlap with the quantification window."
                crispresso2_info['results']['refs'][ref_name]['plot_4d_data'] = []


                ###############################################################################################################################################
                #4e : for HDR, global modifications with respect to reference 1
                ###############################################################################################################################################

                if args.expected_hdr_amplicon_seq != "" and (ref_name == ref_names[0] or ref_name == "HDR"):
                    plot_4e_input = {
                        'ref1_all_insertion_count_vectors': ref1_all_insertion_count_vectors[ref_name],
                        'ref1_all_deletion_count_vectors': ref1_all_deletion_count_vectors[ref_name],
                        'ref1_all_substitution_count_vectors': ref1_all_substitution_count_vectors[ref_name],
                        'ref1': refs[ref_names[0]],
                        'include_idxs_list': include_idxs_list,
                        'n_total': N_TOTAL,
                        'ref_len': ref_len,
                        'ref_name': ref_names[0],
                        'save_also_png': save_png,
                    }
                    if ref_name == ref_names[0]:
                        plot_4e_input['plot_title'] = 'Mutation position distribution in all reads with reference to %s' % (ref_names[0])
                        plot_4e_input['plot_root'] = _jp(
                            '4e.' + ref_names[0] + '.Global_mutations_in_all_reads',
                        )
                        crispresso2_info['results']['refs'][ref_names[0]]['plot_4e_root'] = os.path.basename(plot_root)
                        crispresso2_info['results']['refs'][ref_names[0]]['plot_4e_caption'] = "Figure 4e: Positions of modifications in all reads when aligned to the reference sequence ("+ref_names[0]+"). Insertions: red, deletions: purple, substitutions: green. All modifications (including those outside the quantification window) are shown."
                        crispresso2_info['results']['refs'][ref_names[0]]['plot_4e_data'] = []
                    elif ref_name == "HDR":
                        plot_4e_input['plot_title'] = 'Mutation position distribution in %s reads with reference to %s'%(ref_name, ref_names[0])
                        plot_4e_input['plot_root'] = _jp(
                            '4f.' + ref_names[0] + '.Global_mutations_in_HDR_reads_with_reference_to_'+ref_names[0],
                        )
                        crispresso2_info['results']['refs'][ref_names[0]]['plot_4f_root'] = os.path.basename(plot_root)
                        crispresso2_info['results']['refs'][ref_names[0]]['plot_4f_caption'] = "Figure 4f: Positions of modifications in HDR reads with respect to the reference sequence ("+ref_names[0]+"). Insertions: red, deletions: purple, substitutions: green. All modifications (including those outside the quantification window) are shown."
                        crispresso2_info['results']['refs'][ref_names[0]]['plot_4f_data'] = []
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_global_modifications_reference,
                            **plot_4e_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_global_modifications_reference(
                            **plot_4e_input,
                        )

                ###############################################################################################################################################
                #4g : for HDR, nuc quilt comparison
                ###############################################################################################################################################

                if args.expected_hdr_amplicon_seq != "" and ref_name == ref_names[0]:
                    nuc_pcts = []
                    ref_names_for_hdr = [r for r in ref_names if counts_total[r] > 0]
                    for ref_name_for_hdr in ref_names_for_hdr:
                        tot = float(counts_total[ref_name_for_hdr])
                        for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                            nuc_pcts.append(np.concatenate(([ref_name_for_hdr, nuc], np.array(ref1_all_base_count_vectors[ref_name_for_hdr+"_"+nuc]).astype(np.float)/tot)))
                    colnames = ['Batch', 'Nucleotide']+list(refs[ref_names_for_hdr[0]]['sequence'])
                    hdr_nucleotide_percentage_summary_df = pd.DataFrame(nuc_pcts, columns=colnames).apply(pd.to_numeric, errors='ignore')

                    mod_pcts = []
                    for ref_name_for_hdr in ref_names_for_hdr:
                        tot = float(counts_total[ref_name_for_hdr])
                        mod_pcts.append(np.concatenate(([ref_name_for_hdr, 'Insertions'], np.array(ref1_all_insertion_count_vectors[ref_name_for_hdr]).astype(np.float)/tot)))
                        mod_pcts.append(np.concatenate(([ref_name_for_hdr, 'Insertions_Left'], np.array(ref1_all_insertion_left_count_vectors[ref_name_for_hdr]).astype(np.float)/tot)))
                        mod_pcts.append(np.concatenate(([ref_name_for_hdr, 'Deletions'], np.array(ref1_all_deletion_count_vectors[ref_name_for_hdr]).astype(np.float)/tot)))
                        mod_pcts.append(np.concatenate(([ref_name_for_hdr, 'Substitutions'], np.array(ref1_all_substitution_count_vectors[ref_name_for_hdr]).astype(np.float)/tot)))
                        mod_pcts.append(np.concatenate(([ref_name_for_hdr, 'All_modifications'], np.array(ref1_all_indelsub_count_vectors[ref_name_for_hdr]).astype(np.float)/tot)))
                        mod_pcts.append(np.concatenate(([ref_name_for_hdr, 'Total'], [counts_total[ref_name_for_hdr]]*refs[ref_names_for_hdr[0]]['sequence_length'])))
                    colnames = ['Batch', 'Modification']+list(refs[ref_names_for_hdr[0]]['sequence'])
                    hdr_modification_percentage_summary_df = pd.DataFrame(mod_pcts, columns=colnames).apply(pd.to_numeric, errors='ignore')
                    sgRNA_intervals = refs[ref_names_for_hdr[0]]['sgRNA_intervals']
                    sgRNA_names = refs[ref_names_for_hdr[0]]['sgRNA_names']
                    sgRNA_mismatches = refs[ref_names_for_hdr[0]]['sgRNA_mismatches']
#                    include_idxs_list = refs[ref_names_for_hdr[0]]['include_idxs']
                    include_idxs_list = [] # the quantification windows may be different between different amplicons

                    ref_plot_name = refs[ref_names_for_hdr[0]]['ref_plot_name']
                    mod_freq_filename = _jp(ref_plot_name + 'Reads_from_all_amplicons_modification_percent_table.txt')
                    hdr_modification_percentage_summary_df.rename(columns={'Batch':'Amplicon'}).to_csv(mod_freq_filename,sep='\t',header=True,index=False)
                    nuc_freq_filename = _jp(ref_plot_name + 'Reads_from_all_amplicons_nucleotide_percent_table.txt')
                    hdr_nucleotide_percentage_summary_df.rename(columns={'Batch':'Amplicon'}).to_csv(nuc_freq_filename,sep='\t',header=True,index=False)

                    plot_root = _jp('4g.HDR_nucleotide_percentage_quilt')
                    plot_4g_input = {
                        'nuc_pct_df': hdr_nucleotide_percentage_summary_df,
                        'mod_pct_df': hdr_modification_percentage_summary_df,
                        'fig_filename_root': plot_root,
                        'save_also_png': save_png,
                        'sgRNA_intervals': sgRNA_intervals,
                        'quantification_window_idxs': include_idxs_list,
                        'sgRNA_names': sgRNA_names,
                        'sgRNA_mismatches': sgRNA_mismatches,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_nucleotide_quilt,
                            **plot_4g_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_nucleotide_quilt(
                            **plot_4g_input,
                        )
                    crispresso2_info['results']['refs'][ref_names_for_hdr[0]]['plot_4g_root'] = os.path.basename(plot_root)
                    crispresso2_info['results']['refs'][ref_names_for_hdr[0]]['plot_4g_caption'] = "Figure 4g: Nucleotide distribution across all amplicons. At each base in the reference amplicon, the percentage of each base as observed in sequencing reads is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the percentage of reads for which that base was deleted. Brown bars between bases show the percentage of reads having an insertion at that position."
                    crispresso2_info['results']['refs'][ref_names_for_hdr[0]]['plot_4g_data'] = []
                    crispresso2_info['results']['refs'][ref_names_for_hdr[0]]['plot_4g_data'].append(('Nucleotide percentage table for all reads aligned to ' + ref_names[0],os.path.basename(nuc_freq_filename)))
                    crispresso2_info['results']['refs'][ref_names_for_hdr[0]]['plot_4g_data'].append(('Modification percentage table for all reads aligned to ' + ref_names[0],os.path.basename(mod_freq_filename)))
                    for ref_name_for_hdr in ref_names_for_hdr:
                        if 'nuc_freq_filename' in crispresso2_info['results']['refs'][ref_name_for_hdr]:
                            crispresso2_info['results']['refs'][ref_names_for_hdr[0]]['plot_4g_data'].append(('Nucleotide frequency table for ' + ref_name_for_hdr, os.path.basename(crispresso2_info['results']['refs'][ref_name_for_hdr]['nuc_freq_filename'])))


            ###############################################################################################################################################
            #(5, 6) frameshift analyses plots
            if (refs[ref_name]['contains_coding_seq']): #PERFORM FRAMESHIFT ANALYSIS
                #make frameshift plots
                ref_len = refs[ref_name]['sequence_length']
                cut_points = refs[ref_name]['sgRNA_cut_points']
                plot_cut_points = refs[ref_name]['sgRNA_plot_cut_points']
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
                    if MODIFIED_FRAMESHIFT + MODIFIED_NON_FRAMESHIFT + NON_MODIFIED_NON_FRAMESHIFT > 0:
                        plot_root = _jp(
                            '5.' + ref_plot_name + 'Frameshift_in-frame_mutations_pie_chart',
                        )
                        plot_5_input = {
                            'modified_frameshift': MODIFIED_FRAMESHIFT,
                            'modified_non_frameshift': MODIFIED_NON_FRAMESHIFT,
                            'non_modified_non_frameshift': NON_MODIFIED_NON_FRAMESHIFT,
                            'cut_points': cut_points,
                            'plot_cut_points': plot_cut_points,
                            'sgRNA_intervals': sgRNA_intervals,
                            'exon_intervals': exon_intervals,
                            'ref_len': ref_len,
                            'ref_name': ref_name,
                            'plot_root': plot_root,
                            'save_also_png': save_png,
                        }
                        if n_processes > 1:
                            plot_results.append(plot_pool.submit(
                                CRISPRessoPlot.plot_frameshift_analysis,
                                **plot_5_input,
                            ))
                        else:
                            CRISPRessoPlot.plot_frameshift_analysis(
                                **plot_5_input,
                            )
                        crispresso2_info['results']['refs'][ref_name]['plot_5_root'] = os.path.basename(plot_root)
                        crispresso2_info['results']['refs'][ref_name]['plot_5_caption'] = "Figure 5: Frameshift analysis of coding sequence reads affected by modifications (unmodified reads are excluded from this analysis)."
                        crispresso2_info['results']['refs'][ref_name]['plot_5_data'] = []

                     #profiles-----------------------------------------------------------------------------------
                    plot_root = _jp(
                        '6.' + ref_plot_name + 'Frameshift_in-frame_mutation_profiles',
                    )
                    plot_6_input = {
                        'hists_frameshift': hists_frameshift[ref_name],
                        'hists_inframe': hists_inframe[ref_name],
                        'plot_titles': {
                            'fs': get_plot_title_with_ref_name(
                                'Frameshift profile', ref_name,
                            ),
                            'if': get_plot_title_with_ref_name(
                                'In-frame profile', ref_name,
                            ),
                        },

                        'plot_root': plot_root,
                        'save_also_png': save_png,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_frameshift_frequency,
                            **plot_6_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_frameshift_frequency(
                            **plot_6_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_6_root'] = os.path.basename(plot_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_6_caption'] = "Figure 6: Frameshift and in-frame mutagenesis profiles indicating position affected by modification. The y axis shows the number of reads and percentage of all reads in that category (frameshifted (top) or in-frame (bottom)). %d reads with no length modifications are not shown."%hists_inframe[ref_name][0]
                    crispresso2_info['results']['refs'][ref_name]['plot_6_data'] = []
                    if 'indel_histogram_filename' in crispresso2_info['results']['refs'][ref_name]:
                        crispresso2_info['results']['refs'][ref_name]['plot_6_data'] = [('Indel histogram for ' + ref_name, os.path.basename(crispresso2_info['results']['refs'][ref_name]['indel_histogram_filename']))]

                    #-----------------------------------------------------------------------------------------------------------

                    #non coding
                    plot_root = _jp('7.'+ref_plot_name+'Insertion_deletion_substitution_locations_noncoding')
                    plot_7_input = {
                        'insertion_count_vectors_noncoding': insertion_count_vectors_noncoding[ref_name],
                        'deletion_count_vectors_noncoding': deletion_count_vectors_noncoding[ref_name],
                        'substitution_count_vectors_noncoding': substitution_count_vectors_noncoding[ref_name],
                        'include_idxs_list': include_idxs_list,
                        'cut_points': cut_points,
                        'plot_cut_points': plot_cut_points,
                        'ref_len': ref_len,
                        'sgRNA_intervals': sgRNA_intervals,
                        'plot_title': get_plot_title_with_ref_name(
                            'Noncoding mutation position distribution',
                            ref_name,
                        ),
                        'plot_root': plot_root,
                        'save_also_png': save_png,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_non_coding_mutations,
                            **plot_7_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_non_coding_mutations(
                            **plot_7_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_7_root'] = os.path.basename(plot_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_7_caption'] = "Figure 7: Reads with insertions (red), deletions (purple), and substitutions (green) mapped to reference amplicon position exclusively in noncoding region/s (that is, without mutations affecting coding sequences). The predicted cleavage site is indicated by a vertical dashed line. Only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are plotted."
                    crispresso2_info['results']['refs'][ref_name]['plot_7_data'] = []

                    plot_root = _jp('8.'+ref_plot_name+'Potential_splice_sites_pie_chart')
                    plot_8_input = {
                        'splicing_sites_modified': SPLICING_SITES_MODIFIED,
                        'count_total': count_total,
                        'plot_root': plot_root,
                        'save_also_png': save_png,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_potential_splice_sites,
                            **plot_8_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_potential_splice_sites(
                            **plot_8_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_8_root'] = os.path.basename(plot_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_8_caption'] = "Figure 8: Predicted impact on splice sites. Potential splice sites modified refers to reads in which the either of the two intronic positions adjacent to exon junctions are disrupted."
                    crispresso2_info['results']['refs'][ref_name]['plot_8_data'] = []

            #end contains coding seq

            ######PLOT
            if not args.crispresso1_mode and args.base_editor_output:
                if not args.suppress_plots:
                    fig_filename_root= _jp('10a.'+ref_plot_name+'Substitution_frequencies_at_each_bp')
                    plot_10a_input = {
                        'ref_len': ref_len,
                        'ref_seq': ref_seq,
                        'ref_name': ref_name,
                        'ref_count': tot_aln_reads,
                        'all_substitution_base_vectors': all_substitution_base_vectors,
                        'plot_title': get_plot_title_with_ref_name(
                            'Substitution frequency', ref_name,
                        ),
                        'fig_filename_root': fig_filename_root,
                        'save_also_png': save_png,
                        'quantification_window_idxs': include_idxs_list,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_subs_across_ref,
                            **plot_10a_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_subs_across_ref(
                            **plot_10a_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_10a_root'] = os.path.basename(fig_filename_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_10a_caption'] = "Figure 10a: Substitution frequencies across the amplicon."
                    if 'nuc_freq_filename' in crispresso2_info['results']['refs'][ref_name]:
                        nuc_freq_filename = crispresso2_info['results']['refs'][ref_name]['nuc_freq_filename']
                        crispresso2_info['results']['refs'][ref_name]['plot_10a_data'] = [('Nucleotide frequencies', os.path.basename(nuc_freq_filename))]

                    #plot all substitution rates in entire region
                    fig_filename_root = _jp('10b.'+ref_plot_name+'Substitution_frequency_barplot')
                    plot_10b_input = {
                        'alt_nuc_counts': alt_nuc_counts_all,
                        'plot_title': get_plot_title_with_ref_name(
                            'Substitution frequency\nin entire amplicon', ref_name,
                        ),
                        'fig_filename_root': fig_filename_root,
                        'save_also_png': save_png
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_sub_freqs,
                            **plot_10b_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_sub_freqs(
                            **plot_10b_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_10b_root'] = os.path.basename(fig_filename_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_10b_caption'] = "Figure 10b: Substitution frequencies across the amplicon."
                    crispresso2_info['results']['refs'][ref_name]['plot_10b_data'] = [('Nucleotide frequencies', os.path.basename(nuc_freq_filename))]

                    #plot all substitution rates in quantification_window
                    fig_filename_root = _jp('10c.'+ref_plot_name+'Substitution_frequency_barplot_in_quantification_window')
                    plot_10c_input = {
                        'alt_nuc_counts': alt_nuc_counts,
                        'plot_title': get_plot_title_with_ref_name('Substitution frequency\nin quantification window', ref_name),
                        'fig_filename_root': fig_filename_root,
                        'save_also_png': save_png
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_sub_freqs,
                            **plot_10c_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_sub_freqs(
                            **plot_10c_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_10c_root'] = os.path.basename(fig_filename_root)
                    crispresso2_info['results']['refs'][ref_name]['plot_10c_caption'] = "Figure 10c: Substitution frequencies in the quantification window"
                    crispresso2_info['results']['refs'][ref_name]['plot_10c_data'] = [('Nucleotide frequencies in quantification window', os.path.basename(quant_window_sub_freq_filename))]

            ##new plots alleles around cut_sites
            sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
            sgRNA_cut_points = refs[ref_name]['sgRNA_cut_points']
            sgRNA_plot_cut_points = refs[ref_name]['sgRNA_plot_cut_points']
            sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
            sgRNA_names = refs[ref_name]['sgRNA_names']
            sgRNA_plot_idxs = refs[ref_name]['sgRNA_plot_idxs']
            sgRNA_mismatches = refs[ref_name]['sgRNA_mismatches']

            crispresso2_info['results']['refs'][ref_name]['plot_9_roots'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_9_captions'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_9_datas'] = []
            crispresso2_info['results']['refs'][ref_name]['allele_frequency_files'] = []

            crispresso2_info['results']['refs'][ref_name]['plot_10d_roots'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10d_captions'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10d_datas'] = []

            crispresso2_info['results']['refs'][ref_name]['plot_10e_roots'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10e_captions'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10e_datas'] = []

            crispresso2_info['results']['refs'][ref_name]['plot_10f_roots'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10f_captions'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10f_datas'] = []

            crispresso2_info['results']['refs'][ref_name]['plot_10g_roots'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10g_captions'] = []
            crispresso2_info['results']['refs'][ref_name]['plot_10g_datas'] = []

            for ind, sgRNA_seq in enumerate(sgRNA_sequences):
                cut_point = sgRNA_cut_points[ind]
                plot_cut_point = sgRNA_plot_cut_points[ind]
                sgRNA_name = sgRNA_names[ind]
                plot_idxs = sgRNA_plot_idxs[ind]
                sgRNA = sgRNA_orig_sequences[ind]

                sgRNA_label = "sgRNA_"+sgRNA # for file names
                sgRNA_legend = "sgRNA " + sgRNA # for legends
                if sgRNA_name != "":
                    sgRNA_label = sgRNA_name
                    sgRNA_legend = sgRNA_name + " (" + sgRNA +")"
                sgRNA_label = CRISPRessoShared.slugify(sgRNA_label)

                plot_half_window = max(1, args.plot_window_size)
                df_alleles_around_cut=CRISPRessoShared.get_dataframe_around_cut(df_alleles.loc[df_alleles['Reference_Name'] == ref_name], cut_point, plot_half_window)
                count_total = counts_total[ref_name]
                if args.allele_plot_pcts_only_for_assigned_reference:
                    df_alleles_around_cut['%AllReads']=df_alleles_around_cut['%Reads']
                    df_alleles_around_cut['%Reads']=df_alleles_around_cut['#Reads']/count_total*100

                #write alleles table to file
                allele_filename = _jp(ref_plot_name+'Alleles_frequency_table_around_'+sgRNA_label+'.txt')
                df_alleles_around_cut.to_csv(allele_filename, sep='\t', header=True)
                crispresso2_info['results']['refs'][ref_name]['allele_frequency_files'].append(os.path.basename(allele_filename))

                ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_half_window+1:cut_point+plot_half_window+1]
                fig_filename_root = _jp('9.'+ref_plot_name+'Alleles_frequency_table_around_'+sgRNA_label)
                n_good = df_alleles_around_cut[df_alleles_around_cut['%Reads']>=args.min_frequency_alleles_around_cut_to_plot].shape[0]
                if not args.suppress_plots and n_good > 0:

                    df_to_plot = df_alleles_around_cut
                    if not args.expand_allele_plots_by_quantification:
                        df_to_plot = df_alleles_around_cut.groupby(['Aligned_Sequence', 'Reference_Sequence']).sum().reset_index().set_index('Aligned_Sequence')
                        df_to_plot.sort_values(by='%Reads', inplace=True, ascending=False)

                    new_sgRNA_intervals = []
                    #adjust coordinates of sgRNAs
                    new_sel_cols_start = cut_point - plot_half_window
                    for (int_start, int_end) in refs[ref_name]['sgRNA_intervals']:
                        new_sgRNA_intervals += [(int_start - new_sel_cols_start - 1, int_end - new_sel_cols_start - 1)]
                    plot_9_input = {
                        'reference_seq': ref_seq_around_cut,
                        'df_alleles': df_to_plot,
                        'fig_filename_root': fig_filename_root,
                        'MIN_FREQUENCY': args.min_frequency_alleles_around_cut_to_plot,
                        'MAX_N_ROWS': args.max_rows_alleles_around_cut_to_plot,
                        'SAVE_ALSO_PNG': save_png,
                        'plot_cut_point': plot_cut_point,
                        'sgRNA_intervals': new_sgRNA_intervals,
                        'sgRNA_names': sgRNA_names,
                        'sgRNA_mismatches': sgRNA_mismatches,
                        'annotate_wildtype_allele': args.annotate_wildtype_allele,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_alleles_table,
                            **plot_9_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_alleles_table(
                            **plot_9_input,
                        )
                    crispresso2_info['results']['refs'][ref_name]['plot_9_roots'].append(os.path.basename(fig_filename_root))
                    crispresso2_info['results']['refs'][ref_name]['plot_9_captions'].append("Figure 9: Visualization of the distribution of identified alleles around the cleavage site for the " + sgRNA_legend + ". Nucleotides are indicated by unique colors (A = green; C = red; G = yellow; T = purple). Substitutions are shown in bold font. Red rectangles highlight inserted sequences. Horizontal dashed lines indicate deleted sequences. The vertical dashed line indicates the predicted cleavage site.")
                    crispresso2_info['results']['refs'][ref_name]['plot_9_datas'].append([('Allele frequency table', os.path.basename(allele_filename))])

                if not args.crispresso1_mode and args.base_editor_output:
                    ###guide-specific base editor plots
                    plot_ref_seq = ref_seq_around_cut
                    plot_nuc_pcts = df_nuc_pct_all.iloc[:, plot_idxs]
                    plot_nuc_freqs = df_nuc_freq_all.iloc[:, plot_idxs]

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
                    just_sel_nuc_pcts = plot_nuc_pcts.iloc[:, from_nuc_indices].copy() #only nucleotides targeted by base editing
                    #just_sel_nuc_pcts.columns = [char + str(pos+1) for pos,char in enumerate(list(just_sel_nuc_pcts.columns.values))]
                    just_sel_nuc_pcts.columns = [args.conversion_nuc_from + str(pos+1) for pos in from_nuc_indices]
                    just_sel_nuc_freqs = plot_nuc_freqs.iloc[:, from_nuc_indices].copy()
                    just_sel_nuc_freqs.columns = [args.conversion_nuc_from + str(pos+1) for pos in from_nuc_indices]

                    quant_window_sel_nuc_pct_filename = _jp(ref_plot_name + 'Selected_nucleotide_percentage_table_around_'+sgRNA_label+'.txt')
                    just_sel_nuc_pcts.to_csv(quant_window_sel_nuc_pct_filename, sep='\t', header=True, index=True)
#                   not storing the name because it is unique to this sgRNA
#                    crispresso2_info['quant_window_sel_nuc_pct_filename'] = os.path.basename(quant_window_sel_nuc_pct_filename)

                    quant_window_sel_nuc_freq_filename = _jp(ref_plot_name + 'Selected_nucleotide_frequency_table_around_'+sgRNA_label+'.txt')
                    just_sel_nuc_freqs.to_csv(quant_window_sel_nuc_freq_filename, sep='\t', header=True, index=True)
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
                        fig_filename_root = _jp(
                            '10d.' + ref_plot_name + 'Log2_nucleotide_frequency_around_' + sgRNA_label,
                        )
                        plot_10d_input = {
                            'df_nuc_freq': plot_nuc_freqs,
                            'tot_aln_reads': tot_aln_reads,
                            'plot_title': get_plot_title_with_ref_name(
                                'Log2 Nucleotide Frequencies Around the ' + sgRNA_legend,
                                ref_name,
                            ),
                            'fig_filename_root': fig_filename_root,
                            'save_also_png': save_png,
                            'quantification_window_idxs': plot_quant_window_idxs,
                        }
                        if n_processes > 1:
                            plot_results.append(plot_pool.submit(
                                CRISPRessoPlot.plot_log_nuc_freqs,
                                **plot_10d_input,
                            ))
                        else:
                            CRISPRessoPlot.plot_log_nuc_freqs(
                                **plot_10d_input,
                            )
                        crispresso2_info['results']['refs'][ref_name]['plot_10d_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['results']['refs'][ref_name]['plot_10d_captions'].append("Figure 10d: Log2 nucleotide frequencies for each position in the plotting window around the " + sgRNA_legend + ". The quantification window is outlined by the dotted box.")
                        crispresso2_info['results']['refs'][ref_name]['plot_10d_datas'].append([])

                        fig_filename_root = _jp(
                            '10e.' + ref_plot_name + 'Selected_conversion_at_' + args.conversion_nuc_from + 's_around_' + sgRNA_label,
                        )
                        plot_10e_input = {
                            'df_subs': plot_nuc_pcts,
                            'ref_name': ref_name,
                            'ref_sequence': plot_ref_seq,
                            'plot_title': get_plot_title_with_ref_name(
                                'Substitution Frequencies at ' + args.conversion_nuc_from + 's around the ' + sgRNA_legend,
                                ref_name,
                            ),
                            'conversion_nuc_from': args.conversion_nuc_from,
                            'fig_filename_root': fig_filename_root,
                            'save_also_png': save_png,
                        }
                        if n_processes > 1:
                            plot_results.append(plot_pool.submit(
                                CRISPRessoPlot.plot_conversion_at_sel_nucs,
                                **plot_10e_input,
                            ))
                        else:
                            CRISPRessoPlot.plot_conversion_at_sel_nucs(
                                **plot_10e_input,
                            )
                        crispresso2_info['results']['refs'][ref_name]['plot_10e_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['results']['refs'][ref_name]['plot_10e_captions'].append("Figure 10e: Proportion of each base at each nucleotide targeted by base editors in the plotting window around the " + sgRNA_legend + ". The number of each target base is annotated on the reference sequence at the bottom of the plot.")
                        crispresso2_info['results']['refs'][ref_name]['plot_10e_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from + 's', os.path.basename(quant_window_sel_nuc_freq_filename))])

                        fig_filename_root = _jp('10f.'+ref_plot_name+'Selected_conversion_no_ref_at_'+args.conversion_nuc_from+'s_around_'+sgRNA_label)
                        plot_10f_input = {
                            'df_subs': plot_nuc_pcts,
                            'ref_name': ref_name,
                            'ref_sequence': plot_ref_seq,
                            'plot_title': get_plot_title_with_ref_name('Substitution Frequencies at '+args.conversion_nuc_from+'s around the ' + sgRNA_legend, ref_name),
                            'conversion_nuc_from': args.conversion_nuc_from,
                            'fig_filename_root': fig_filename_root,
                            'save_also_png': save_png,
                        }
                        if n_processes > 1:
                            plot_results.append(plot_pool.submit(
                                CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref,
                                **plot_10f_input,
                            ))
                        else:
                            CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref(
                                **plot_10f_input,
                            )
                        crispresso2_info['results']['refs'][ref_name]['plot_10f_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['results']['refs'][ref_name]['plot_10f_captions'].append("Figure 10f: Non-reference base proportions. For target nucleotides in the plotting window, this plot shows the proportion of non-reference (non-"+args.conversion_nuc_from + ") bases as a percentage of all non-reference sequences. The number of each target base is annotated on the reference sequence at the bottom of the plot.")
                        crispresso2_info['results']['refs'][ref_name]['plot_10f_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from + 's', os.path.basename(quant_window_sel_nuc_freq_filename))])

                        fig_filename_root = _jp('10g.'+ref_plot_name+'Selected_conversion_no_ref_scaled_at_'+args.conversion_nuc_from+'s_around_'+sgRNA_label)
                        plot_10g_input = {
                            'df_subs': plot_nuc_pcts,
                            'ref_name': ref_name,
                            'ref_sequence': plot_ref_seq,
                            'plot_title': get_plot_title_with_ref_name(
                                'Substitution Frequencies at ' + args.conversion_nuc_from + 's around the ' + sgRNA_legend,
                                ref_name,
                            ),
                            'conversion_nuc_from': args.conversion_nuc_from,
                            'fig_filename_root': fig_filename_root,
                            'save_also_png': save_png
                        }
                        if n_processes > 1:
                            plot_results.append(plot_pool.submit(
                                CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref_scaled,
                                **plot_10g_input,
                            ))
                        else:
                            CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref_scaled(
                                **plot_10g_input,
                            )
                        crispresso2_info['results']['refs'][ref_name]['plot_10g_roots'].append(os.path.basename(fig_filename_root))
                        crispresso2_info['results']['refs'][ref_name]['plot_10g_captions'].append("Figure 10g: Non-reference base counts. For target nucleotides in the plotting window, this plot shows the number of non-reference (non-" + args.conversion_nuc_from + ") bases. The number of each target base is annotated on the reference sequence at the bottom of the plot.")
                        crispresso2_info['results']['refs'][ref_name]['plot_10g_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from +'s', os.path.basename(quant_window_sel_nuc_freq_filename))])

            #END GUIDE SPECIFIC PLOTS

        #(5, 6) GLOBAL frameshift analyses plots
        if args.coding_seq:
            global_MODIFIED_FRAMESHIFT = 0
            global_MODIFIED_NON_FRAMESHIFT = 0
            global_NON_MODIFIED_NON_FRAMESHIFT = 0
            global_SPLICING_SITES_MODIFIED = 0

            global_hists_frameshift = defaultdict(lambda :0)
            global_hists_frameshift[0] = 0  # fill with at least the zero value (in case there are no others)
            global_hists_inframe = defaultdict(lambda :0)
            global_hists_inframe[0] = 0

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

                    for (exon_len, count) in hists_frameshift[ref_name].items():
                        global_hists_frameshift[exon_len] += count
                    for (exon_len, count) in hists_inframe[ref_name].items():
                        global_hists_inframe[exon_len] += count


                    global_count_total += counts_total[ref_name]
                    global_count_modified += counts_modified[ref_name]
                    global_count_unmodified += counts_unmodified[ref_name]

            if not args.suppress_plots:
                plot_root = _jp('5a.Global_frameshift_in-frame_mutations_pie_chart')
                plot_5a_input = {
                    'global_modified_frameshift': global_MODIFIED_FRAMESHIFT,
                    'global_modified_non_frameshift': global_MODIFIED_NON_FRAMESHIFT,
                    'global_non_modified_non_frameshift': global_NON_MODIFIED_NON_FRAMESHIFT,
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_global_frameshift_analysis,
                        **plot_5a_input,
                    ))
                else:
                    CRISPRessoPlot.plot_global_frameshift_analysis(
                        **plot_5a_input,
                    )
                crispresso2_info['results']['general_plots']['plot_5a_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['general_plots']['plot_5a_caption'] = "Figure 5a: Frameshift analysis of coding sequence reads affected by modifications for all reads. Unmodified reference reads are excluded from this plot, and all HDR reads are included in this plot."
                crispresso2_info['results']['general_plots']['plot_5a_data'] = []

                 #profiles-----------------------------------------------------------------------------------
                plot_root = _jp('6a.Global_frameshift_in-frame_mutation_profiles')
                plot_6a_input = {
                    'global_hists_frameshift': global_hists_frameshift,
                    'global_hists_inframe': global_hists_inframe,
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_global_frameshift_in_frame_mutations,
                        **plot_6a_input,
                    ))
                else:
                    CRISPRessoPlot.plot_global_frameshift_in_frame_mutations(
                        **plot_6a_input,
                    )

                crispresso2_info['results']['general_plots']['plot_6a_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['general_plots']['plot_6a_caption'] = "Figure 6a: Frameshift and in-frame mutagenesis profiles for all reads indicating position affected by modification. The y axis shows the number of reads and percentage of all reads in that category (frameshifted (top) or in-frame (bottom)). %d reads with no length modifications are not shown."%global_hists_inframe[0]
                crispresso2_info['results']['general_plots']['plot_6a_data'] = []
                for ref_name in ref_names:
                    if 'indel_histogram_filename' in crispresso2_info['results']['refs'][ref_name]:
                        crispresso2_info['results']['general_plots']['plot_6a_data'].append(('Indel histogram for ' + ref_name, os.path.basename(crispresso2_info['results']['refs'][ref_name]['indel_histogram_filename'])))


                #-----------------------------------------------------------------------------------------------------------
                plot_root = _jp('8a.Global_potential_splice_sites_pie_chart')
                plot_8a_input = {
                    'global_splicing_sites_modified': global_SPLICING_SITES_MODIFIED,
                    'global_count_total': global_count_total,
                    'plot_root': plot_root,
                    'save_also_png': save_png,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_impact_on_splice_sites,
                        **plot_8a_input,
                    ))
                else:
                    CRISPRessoPlot.plot_impact_on_splice_sites(
                        **plot_8a_input,
                    )
                crispresso2_info['results']['general_plots']['plot_8a_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['general_plots']['plot_8a_caption'] = "Figure 8a: Predicted impact on splice sites for all reads. Potential splice sites modified refers to reads in which the either of the two intronic positions adjacent to exon junctions are disrupted."
                crispresso2_info['results']['general_plots']['plot_8a_data'] = []

            #end global coding seq plots

        #prime editing plots
        if args.prime_editing_pegRNA_extension_seq != "":
            # count length of scaffold insertions
            scaffold_insertion_sizes_filename = ""
            if args.prime_editing_pegRNA_scaffold_seq != "":
                #first, define the sequence we are looking for (extension plus the first base(s) of the scaffold)
                scaffold_dna_seq = CRISPRessoShared.reverse_complement(args.prime_editing_pegRNA_scaffold_seq)
                pe_seq = refs['Prime-edited']['sequence']
                pe_scaffold_dna_info = get_pe_scaffold_search(pe_seq, args.prime_editing_pegRNA_extension_seq, args.prime_editing_pegRNA_scaffold_seq, args.prime_editing_pegRNA_scaffold_min_match_length)

                df_alleles_scaffold = df_alleles.loc[df_alleles['Reference_Name'] == 'Scaffold-incorporated']

                def get_scaffold_len(row, scaffold_start_loc, scaffold_seq):
                    pe_read_possible_scaffold_loc = row['ref_positions'].index(scaffold_start_loc - 1)+1
                    i = 0
                    num_match_scaffold = 0
                    num_gaps = 0
                    matches_scaffold_to_this_point = True
                    has_gaps_to_this_point = True
                    aln_seq_to_test = row['Aligned_Sequence'][pe_read_possible_scaffold_loc:].replace("-", "")
                    while i < len(scaffold_seq) and (matches_scaffold_to_this_point or has_gaps_to_this_point):
                        if matches_scaffold_to_this_point and aln_seq_to_test[i] == scaffold_seq[i] :
                            num_match_scaffold += 1
                        else:
                            matches_scaffold_to_this_point = False

                        if has_gaps_to_this_point and row['Reference_Sequence'][pe_read_possible_scaffold_loc + i] == '-':
                            num_gaps += 1
                        else:
                            has_gaps_to_this_point = False
                        i += 1
                    return row['Aligned_Sequence'], row['Reference_Sequence'], num_match_scaffold, num_gaps, row['#Reads'], row['%Reads']

                df_scaffold_insertion_sizes = pd.DataFrame(list(df_alleles_scaffold.apply(lambda row: get_scaffold_len(row, pe_scaffold_dna_info[0], scaffold_dna_seq), axis=1).values),
                    columns=['Aligned_Sequence', 'Reference_Sequence', 'Num_match_scaffold', 'Num_gaps', '#Reads', '%Reads'])

                scaffold_insertion_sizes_filename = _jp('Scaffold_insertion_sizes.txt')
                df_scaffold_insertion_sizes.to_csv(scaffold_insertion_sizes_filename, sep='\t', header=True, index=False)
                crispresso2_info['Scaffold_insertion_sizes_filename'] = scaffold_insertion_sizes_filename

            #if pegRNA extension seq and plotting, plot summary plots
            if not args.suppress_plots:
                nuc_pcts = []
                ref_names_for_pe = [r for r in ref_names if counts_total[r] > 0]
                for ref_name in ref_names_for_pe:
                    tot = float(counts_total[ref_name])
                    for nuc in ['A', 'C', 'G', 'T', 'N', '-']:
                        nuc_pcts.append(np.concatenate(([ref_name, nuc], np.array(ref1_all_base_count_vectors[ref_name+"_"+nuc]).astype(np.float)/tot)))
                colnames = ['Batch', 'Nucleotide']+list(refs[ref_names[0]]['sequence'])
                pe_nucleotide_percentage_summary_df = pd.DataFrame(nuc_pcts, columns=colnames).apply(pd.to_numeric,errors='ignore')

                mod_pcts = []
                for ref_name in ref_names_for_pe:
                    tot = float(counts_total[ref_name])
                    mod_pcts.append(np.concatenate(([ref_name, 'Insertions'], np.array(ref1_all_insertion_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate(([ref_name, 'Insertions_Left'], np.array(ref1_all_insertion_left_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate(([ref_name, 'Deletions'], np.array(ref1_all_deletion_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate(([ref_name, 'Substitutions'], np.array(ref1_all_substitution_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate(([ref_name, 'All_modifications'], np.array(ref1_all_indelsub_count_vectors[ref_name]).astype(np.float)/tot)))
                    mod_pcts.append(np.concatenate(([ref_name, 'Total'], [counts_total[ref_name]]*refs[ref_names_for_pe[0]]['sequence_length'])))
                colnames = ['Batch', 'Modification']+list(refs[ref_names_for_pe[0]]['sequence'])
                pe_modification_percentage_summary_df = pd.DataFrame(mod_pcts, columns=colnames).apply(pd.to_numeric,errors='ignore')
                sgRNA_intervals = refs[ref_names_for_pe[0]]['sgRNA_intervals']
                sgRNA_names = refs[ref_names_for_pe[0]]['sgRNA_names']
                sgRNA_mismatches = refs[ref_names_for_pe[0]]['sgRNA_mismatches']
                include_idxs_list = refs[ref_names_for_pe[0]]['include_idxs']

                plot_root = _jp('11a.Prime_editing_nucleotide_percentage_quilt')
                plot_11a_input = {
                    'nuc_pct_df': pe_nucleotide_percentage_summary_df,
                    'mod_pct_df': pe_modification_percentage_summary_df,
                    'fig_filename_root': plot_root,
                    'save_also_png': save_png,
                    'sgRNA_intervals': sgRNA_intervals,
                    'sgRNA_names': sgRNA_names,
                    'sgRNA_mismatches': sgRNA_mismatches,
                    'quantification_window_idxs': include_idxs_list,
                }
                if n_processes > 1:
                    plot_results.append(plot_pool.submit(
                        CRISPRessoPlot.plot_nucleotide_quilt,
                        **plot_11a_input,
                    ))
                else:
                    CRISPRessoPlot.plot_nucleotide_quilt(
                        **plot_11a_input,
                    )
                crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11a_root'] = os.path.basename(plot_root)
                crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11a_caption'] = "Figure 11a: Nucleotide distribution across all amplicons. At each base in the reference amplicon, the percentage of each base as observed in sequencing reads is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the percentage of reads for which that base was deleted. Brown bars between bases show the percentage of reads having an insertion at that position."
                crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11a_data'] = [('Nucleotide frequency table for ' + ref_name, os.path.basename(crispresso2_info['results']['refs'][ref_name]['nuc_freq_filename'])) for ref_name in ref_names_for_pe]

                crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11b_roots'] = []
                crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11b_captions'] = []
                crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11b_datas'] = []

                pe_sgRNA_sequences = refs[ref_names_for_pe[0]]['sgRNA_sequences']
                pe_sgRNA_orig_sequences = refs[ref_names_for_pe[0]]['sgRNA_orig_sequences']
                pe_sgRNA_cut_points = refs[ref_names_for_pe[0]]['sgRNA_cut_points']
                pe_sgRNA_plot_cut_points = refs[ref_names_for_pe[0]]['sgRNA_plot_cut_points']
                pe_sgRNA_intervals = refs[ref_names_for_pe[0]]['sgRNA_intervals']
                pe_sgRNA_names = refs[ref_names_for_pe[0]]['sgRNA_names']
                pe_sgRNA_plot_idxs = refs[ref_names_for_pe[0]]['sgRNA_plot_idxs']
                pe_sgRNA_mismatches = refs[ref_names_for_pe[0]]['sgRNA_mismatches']
                pe_ref_len = refs[ref_names_for_pe[0]]['sequence_length']
                pe_include_idxs_list = refs[ref_names_for_pe[0]]['include_idxs']

                for i in range(len(pe_sgRNA_cut_points)):
                    cut_point = pe_sgRNA_cut_points[i]
                    sgRNA = pe_sgRNA_orig_sequences[i]
                    sgRNA_name = pe_sgRNA_names[i]

                    sgRNA_label = "sgRNA_"+sgRNA # for file names
                    sgRNA_legend = "sgRNA " + sgRNA # for legends
                    if sgRNA_name != "":
                        sgRNA_label = sgRNA_name
                        sgRNA_legend = sgRNA_name + " (" + sgRNA +")"
                    sgRNA_label = CRISPRessoShared.slugify(sgRNA_label)

                    #get nucleotide columns to print for this sgRNA
                    sel_cols = [0, 1]
                    plot_half_window = max(1, args.plot_window_size)
                    new_sel_cols_start = max(2, cut_point-plot_half_window+1)
                    new_sel_cols_end = min(pe_ref_len, cut_point+plot_half_window+1)
                    sel_cols.extend(list(range(new_sel_cols_start+2, new_sel_cols_end+2))) #+2 because the first two columns are Batch and Nucleotide
                    #get new intervals
                    new_sgRNA_intervals = []
                    #add annotations for each sgRNA (to be plotted on this sgRNA's plot)
                    for (int_start, int_end) in refs[ref_names_for_pe[0]]['sgRNA_intervals']:
                        new_sgRNA_intervals += [(int_start - new_sel_cols_start, int_end - new_sel_cols_start)]
                    new_include_idx = []
                    for x in pe_include_idxs_list:
                        new_include_idx += [x - new_sel_cols_start]
                    plot_root = _jp('11b.Nucleotide_percentage_quilt_around_' + sgRNA_label)
                    plot_11b_input = {
                        'nuc_pct_df': pe_nucleotide_percentage_summary_df.iloc[:, sel_cols],
                        'mod_pct_df': pe_modification_percentage_summary_df.iloc[:, sel_cols],
                        'fig_filename_root': plot_root,
                        'save_also_png': save_png,
                        'sgRNA_intervals': new_sgRNA_intervals,
                        'sgRNA_names': sgRNA_names,
                        'sgRNA_mismatches': sgRNA_mismatches,
                        'quantification_window_idxs': new_include_idx,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_nucleotide_quilt,
                            **plot_11b_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_nucleotide_quilt(
                            **plot_11b_input,
                        )
                    crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11b_roots'].append(os.path.basename(plot_root))
                    crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11b_captions'].append('Figure 11b: Nucleotide distribution around the ' + sgRNA_legend + '.')
                    crispresso2_info['results']['refs'][ref_names_for_pe[0]]['plot_11b_datas'].append([('Nucleotide frequency in quantification window for ' + ref_name, os.path.basename(crispresso2_info['results']['refs'][ref_name]['quant_window_nuc_freq_filename'])) for ref_name in ref_names_for_pe])

                if args.prime_editing_pegRNA_scaffold_seq != "" and df_scaffold_insertion_sizes.shape[0] > 0 and df_scaffold_insertion_sizes['Num_match_scaffold'].max() > 0 and df_scaffold_insertion_sizes['Num_gaps'].max() > 0:
                    plot_root = _jp('11c.Prime_editing_scaffold_insertion_sizes')
                    plot_11c_input = {
                        'df_scaffold_insertion_sizes': df_scaffold_insertion_sizes,
                        'plot_root': plot_root,
                        'save_also_png': save_png,
                    }
                    if n_processes > 1:
                        plot_results.append(plot_pool.submit(
                            CRISPRessoPlot.plot_scaffold_indel_lengths,
                            **plot_11c_input,
                        ))
                    else:
                        CRISPRessoPlot.plot_scaffold_indel_lengths(
                            **plot_11c_input,
                        )
                    crispresso2_info['results']['general_plots']['plot_11c_root'] = os.path.basename(plot_root)
                    crispresso2_info['results']['general_plots']['plot_11c_caption'] = "Figure 11a: Scaffold insertion lengths and deletion lengths in reads that contain a scaffold insertion. 'Length matching scaffold' shows the number of basepairs immediately after the pegRNA extension sequence that exactly match the scaffold RNA sequence. 'Insertion length' shows the length of the insertion immediately after the pegRNA extension sequence (including bases that do not match the scaffold sequence)."
                    crispresso2_info['results']['general_plots']['plot_11c_data'] = [('Scaffold insertion alleles with insertion sizes', os.path.basename(scaffold_insertion_sizes_filename))]

        # join plotting pool
        wait(plot_results)
        plot_pool.shutdown()

        info('Done!')

        if not args.keep_intermediate:
            info('Removing Intermediate files...')

            for file_to_remove in files_to_remove:
                try:
                    if os.path.islink(file_to_remove):
                        os.unlink(file_to_remove)
                    else:
                        os.remove(file_to_remove)
                except Exception as e:
                    warn('Skipping removal of: %s: %s' %(file_to_remove, e))

        crispresso2_info['results']['alignment_stats']['counts_total'] = counts_total
        crispresso2_info['results']['alignment_stats']['counts_modified'] = counts_modified
        crispresso2_info['results']['alignment_stats']['counts_unmodified'] = counts_unmodified
        crispresso2_info['results']['alignment_stats']['counts_discarded'] = counts_discarded

        crispresso2_info['results']['alignment_stats']['counts_insertion'] = counts_insertion
        crispresso2_info['results']['alignment_stats']['counts_deletion'] = counts_deletion
        crispresso2_info['results']['alignment_stats']['counts_substitution'] = counts_substitution

        crispresso2_info['results']['alignment_stats']['counts_only_insertion'] = counts_only_insertion
        crispresso2_info['results']['alignment_stats']['counts_only_deletion'] = counts_only_deletion
        crispresso2_info['results']['alignment_stats']['counts_only_substitution'] = counts_only_substitution
        crispresso2_info['results']['alignment_stats']['counts_insertion_and_deletion'] = counts_insertion_and_deletion
        crispresso2_info['results']['alignment_stats']['counts_insertion_and_substitution'] = counts_insertion_and_substitution
        crispresso2_info['results']['alignment_stats']['counts_deletion_and_substitution'] = counts_deletion_and_substitution
        crispresso2_info['results']['alignment_stats']['counts_insertion_and_deletion_and_substitution'] = counts_insertion_and_deletion_and_substitution
        crispresso2_info['results']['alignment_stats']['counts_modified_frameshift'] = counts_modified_frameshift
        crispresso2_info['results']['alignment_stats']['counts_modified_non_frameshift'] = counts_modified_non_frameshift
        crispresso2_info['results']['alignment_stats']['counts_non_modified_non_frameshift'] = counts_non_modified_non_frameshift
        crispresso2_info['results']['alignment_stats']['counts_splicing_sites_modified'] = counts_splicing_sites_modified
        crispresso2_info['results']['alignment_stats']['class_counts'] = class_counts

        end_time =  datetime.now()
        end_time_string =  end_time.strftime('%Y-%m-%d %H:%M:%S')
        running_time = end_time - start_time
        running_time_string =  str(running_time)

        crispresso2_info['running_info']['end_time'] = end_time
        crispresso2_info['running_info']['end_time_string'] = end_time_string
        crispresso2_info['running_info']['running_time'] = running_time
        crispresso2_info['running_info']['running_time_string'] = running_time_string


        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_report(crispresso2_info, report_name, OUTPUT_DIRECTORY, _ROOT)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        CRISPRessoShared.write_crispresso_info(
            crispresso2_info_file,
            crispresso2_info,
        )

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
        sys.exit(9)
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
