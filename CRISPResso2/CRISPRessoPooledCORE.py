# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''


import os
import errno
import sys
from copy import deepcopy
from datetime import datetime
import subprocess as sb
import glob
import argparse
import unicodedata
import string
import re
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoReport
from CRISPResso2 import CRISPRessoPlot
import traceback

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

####Support functions###
def get_data(path):
        return os.path.join(_ROOT, 'data', path)

def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoPooled!' % library_name)
                sys.exit(1)


#the dependencies are bowtie2 and samtools
def which(program):
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


def check_samtools():

    cmd_path=which('samtools')
    if cmd_path:
        return True
    else:
        sys.stdout.write('\nCRISPRessoPooled requires samtools')
        sys.stdout.write('\n\nPlease install samtools and add it to your path following the instructions at: http://www.htslib.org/download/')
        return False

def check_bowtie2():

    cmd_path1=which('bowtie2')
    cmd_path2=which('bowtie2-inspect')

    if cmd_path1 and cmd_path2:
        return True
    else:
        sys.stdout.write('\nCRISPRessoPooled requires Bowtie2!')
        sys.stdout.write('\n\nPlease install Bowtie2 and add it to your path following the instructions at: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2')
        return False

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    pd.set_option('display.float_format', '{:20,.2f}'.format)
    pd.set_option('display.max_colwidth', -1)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.float_format')
    pd.reset_option('display.max_colwidth')

#this is overkilling to run for many sequences,
#but for few is fine and effective.
def get_align_sequence(seq,bowtie2_index):

    cmd='''bowtie2 -x  %s -c -U %s''' %(bowtie2_index,seq) + ''' |\
    grep -v '@' | awk '{OFS="\t"; bpstart=$4; split ($6,a,"[MIDNSHP]"); n=0;  bpend=bpstart;\
    for (i=1; i in a; i++){\
      n+=1+length(a[i]); \
      if (substr($6,n,1)=="S"){\
          bpstart-=a[i];\
          if (bpend==$4)\
            bpend=bpstart;\
      } else if( (substr($6,n,1)!="I")  && (substr($6,n,1)!="H") )\
          bpend+=a[i];\
    }if ( ($2 % 32)>=16) print $3,bpstart,bpend,"-",$1,$10,$11;else print $3,bpstart,bpend,"+",$1,$10,$11;}' '''
    p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
    return p.communicate()[0]


#get n_reads and region data from region fastq file (location is pulled from filename)
def summarize_region_fastq_chunk(input_arr):
    ret_val = []
    for input in input_arr:
#        print('doing region ' + str(input))
        region_fastq,uncompressed_reference = input.split(" ")
        #region format: REGION_chr8_1077_1198.fastq.gz
        #But if the chr has underscores, it could look like this:
        #    REGION_chr8_KI270812v1_alt_1077_1198.fastq.gz
        region_info = os.path.basename(region_fastq).replace('.fastq.gz','').replace('.fastq','').split('_')
        chr_string = "_".join(region_info[1:len(region_info)-2]) #in case there are underscores
        region_string='%s:%s-%d' % (chr_string,region_info[-2],int(region_info[-1])-1)
        p = sb.Popen("samtools faidx %s %s | grep -v ^\> | tr -d '\n'" %(uncompressed_reference,region_string), shell=True,stdout=sb.PIPE)
        seq = p.communicate()[0]
        p = sb.Popen(('z' if region_fastq.endswith('.gz') else '' ) +"cat < %s | wc -l" % region_fastq, shell=True,stdout=sb.PIPE)
        n_reads = int(float(p.communicate()[0])/4.0)
        ret_val.append([chr_string] + region_info[-2:]+[region_fastq,n_reads,seq])
    return ret_val

#                                                              Consumes  Consumes
# Op  BAM Description                                             query  reference
# M   0   alignment match (can be a sequence match or mismatch)   yes   yes
# I   1   insertion to the reference                              yes   no
# D   2   deletion from the reference                             no    yes
# N   3   skipped region from the reference                       no    yes
# S   4   soft clipping (clipped sequences present in SEQ)        yes   no
# H   5   hard clipping (clipped sequences NOT present in SEQ)    no    no
# P   6   padding (silent deletion from padded reference)         no    no
# =   7   sequence match                                          yes   yes
# X   8   sequence mismatch                                       yes   yes
def get_read_length_from_cigar(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    """
    read_consuming_ops = ("M", "I", "S", "=", "X")
    result = 0
    ops = re.findall(r'(\d+)(\w)', cigar_string)
    for c in ops:
        length,op = c
        if op in read_consuming_ops:
            result += int(length)
    return result

def get_ref_length_from_cigar(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    reference sequence.
    """
    read_consuming_ops = ("M", "D", "N", "=", "X")
    result = 0
    ops = re.findall(r'(\d+)(\w)', cigar_string)
    for c in ops:
        length,op = c
        if op in read_consuming_ops:
            result += int(length)
    return result

def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     n_reads = int(float(p.communicate()[0])/4.0)
     return n_reads

def get_n_aligned_bam(bam_filename):
     p = sb.Popen("samtools view -F 0x904 -c %s" % bam_filename , shell=True,stdout=sb.PIPE)
     return int(p.communicate()[0])

#get a clean name that we can use for a filename
validFilenameChars = "+-_.() %s%s" % (string.ascii_letters, string.digits)

def clean_filename(filename):
    cleanedFilename = unicodedata.normalize('NFKD', unicode(filename)).encode('ASCII', 'ignore')
    return ''.join(c for c in cleanedFilename if c in validFilenameChars)

def find_overlapping_genes(row,df_genes):
    df_genes_overlapping=df_genes.ix[(df_genes.chrom==row.chr_id) &
                                     (df_genes.txStart<=row.bpend) &
                                     (row.bpstart<=df_genes.txEnd)]
    genes_overlapping=[]

    for idx_g,row_g in df_genes_overlapping.iterrows():
        if 'name' in row_g.keys() and 'name2' in row_g.keys():
            genes_overlapping.append( '%s (%s)' % (row_g.name2,row_g['name']))
        elif '#name' in row_g.keys() and 'name2' in row_g.keys():
            genes_overlapping.append( '%s (%s)' % (row_g.name2,row_g['#name']))
        elif '#name' in row_g.keys():
            genes_overlapping.append( '%s' % (row_g['#name']))
        elif 'name' in row_g.keys():
            genes_overlapping.append( '%s' % (row_g['name']))
        else:
            genes_overlapping.append( '%s' % (row_g[0]))



    row['gene_overlapping']=','.join(genes_overlapping)

    return row


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


def main():
    try:
        start_time =  datetime.now()
        start_time_string =  start_time.strftime('%Y-%m-%d %H:%M:%S')

        description = ['~~~CRISPRessoPooled~~~','-Analysis of CRISPR/Cas9 outcomes from POOLED deep sequencing data-']
        pooled_string = r'''
 _______________________
| __  __  __     __ __  |
||__)/  \/  \|  |_ |  \ |
||   \__/\__/|__|__|__/ |
|_______________________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description,pooled_string))

        parser = CRISPRessoShared.getCRISPRessoArgParser(parserTitle = 'CRISPRessoPooled Parameters',requiredParams={'fastq_r1':True})
        parser.add_argument('-f','--amplicons_file', type=str,  help='Amplicons description file. This file is a tab-delimited text file with up to 5 columns (2 required):\
        \nAMPLICON_NAME:  an identifier for the amplicon (must be unique)\nAMPLICON_SEQUENCE:  amplicon sequence used in the experiment\n\
        \nsgRNA_SEQUENCE (OPTIONAL):  sgRNA sequence used for this amplicon without the PAM sequence. Multiple guides can be given separated by commas and not spaces. If not available enter NA.\
        \nEXPECTED_AMPLICON_AFTER_HDR (OPTIONAL): expected amplicon sequence in case of HDR. If not available enter NA.\
        \nCODING_SEQUENCE (OPTIONAL): Subsequence(s) of the amplicon corresponding to coding sequences. If more than one separate them by commas and not spaces. If not available enter NA.', default='')

        #tool specific optional
        parser.add_argument('--gene_annotations', type=str, help='Gene Annotation Table from UCSC Genome Browser Tables (http://genome.ucsc.edu/cgi-bin/hgTables?command=start), \
        please select as table "knownGene", as output format "all fields from selected table" and as file returned "gzip compressed"', default='')
        parser.add_argument('-p','--n_processes',type=str, help='Specify the number of processes to use for analysis.\
        Please use with caution since increasing this parameter will significantly increase the memory required to run CRISPResso. Can be set to \'max\'.',default='1')
        parser.add_argument('-x','--bowtie2_index', type=str, help='Basename of Bowtie2 index for the reference genome', default='')
        parser.add_argument('--bowtie2_options_string', type=str, help='Override options for the Bowtie2 alignment command',default=' -k 1 --end-to-end -N 0 --np 0 ')
        parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=1000)
        parser.add_argument('--skip_failed',  help='Continue with pooled analysis even if one sample fails',action='store_true')
        parser.add_argument('--skip_reporting_problematic_regions',help='Skip reporting of problematic regions. By default, when both amplicons (-f) and genome (-x) are provided, problematic reads that align to the genome but to positions other than where the amplicons align are reported as problematic',action='store_true')
        parser.add_argument('--crispresso_command', help='CRISPResso command to call',default='CRISPResso')

        args = parser.parse_args()

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = set(['fastq_r1','fastq_r2','amplicon_seq','amplicon_name','output_folder','name'])
        crispresso_options_for_pooled = list(crispresso_options-options_to_ignore)

        files_to_remove = []


        info('Checking dependencies...')

        if check_samtools() and check_bowtie2():
            info('All the required dependencies are present!')
        else:
            sys.exit(1)

        #check files
        CRISPRessoShared.check_file(args.fastq_r1)
        if args.fastq_r2:
            CRISPRessoShared.check_file(args.fastq_r2)

        if args.bowtie2_index:
            if (os.path.isfile(args.bowtie2_index+'.1.bt2l')):
                CRISPRessoShared.check_file(args.bowtie2_index+'.1.bt2l')
            else:
                CRISPRessoShared.check_file(args.bowtie2_index+'.1.bt2')


        if args.amplicons_file:
            CRISPRessoShared.check_file(args.amplicons_file)

        if args.gene_annotations:
            CRISPRessoShared.check_file(args.gene_annotations)

        if args.amplicons_file and not args.bowtie2_index:
            RUNNING_MODE='ONLY_AMPLICONS'
            info('Only the Amplicon description file was provided. The analysis will be perfomed using only the provided amplicons sequences.')

        elif args.bowtie2_index and not args.amplicons_file:
            RUNNING_MODE='ONLY_GENOME'
            info('Only the bowtie2 reference genome index file was provided. The analysis will be perfomed using only genomic regions where enough reads align.')
        elif args.bowtie2_index and args.amplicons_file:
            RUNNING_MODE='AMPLICONS_AND_GENOME'
            info('Amplicon description file and bowtie2 reference genome index files provided. The analysis will be perfomed using the reads that are aligned ony to the amplicons provided and not to other genomic regions.')
        else:
            error('Please provide the amplicons description file (-f or --amplicons_file option) or the bowtie2 reference genome index file (-x or --bowtie2_index option) or both.')
            sys.exit(1)

        n_processes = 1
        if args.n_processes == "max":
            n_processes = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes = int(args.n_processes)


        ####TRIMMING AND MERGING
        get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')

        if not args.name:
                 if args.fastq_r2!='':
                         database_id='%s_%s' % (get_name_from_fasta(args.fastq_r1),get_name_from_fasta(args.fastq_r2))
                 else:
                         database_id='%s' % get_name_from_fasta(args.fastq_r1)

        else:
                 database_id=args.name



        OUTPUT_DIRECTORY='CRISPRessoPooled_on_%s' % database_id

        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory

        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoPooled_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        crispresso2_info_file = os.path.join(OUTPUT_DIRECTORY,'CRISPResso2Pooled_info.pickle')
        crispresso2_info = {} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['version'] = CRISPRessoShared.__version__
        crispresso2_info['args'] = deepcopy(args)

        crispresso2_info['log_filename'] = os.path.basename(log_filename)
        crispresso2_info['finished_steps'] = {}

        #keep track of args to see if it is possible to skip computation steps on rerun
        can_finish_incomplete_run = False
        if args.no_rerun:
            if os.path.exists(crispresso2_info_file):
                previous_run_data = cp.load(open(crispresso2_info_file,'rb'))
                if previous_run_data['version'] == CRISPRessoShared.__version__:
                    args_are_same = True
                    for arg in vars(args):
                        if arg is "no_rerun" or arg is "debug" or arg is "n_processes":
                            continue
                        if arg not in vars(previous_run_data['args']):
                            info('Comparing current run to previous run: old run had argument ' + str(arg) + ' \nRerunning.')
                            args_are_same = False
                        elif str(getattr(previous_run_data['args'],arg)) != str(getattr(args,arg)):
                            info('Comparing current run to previous run:\n\told argument ' + str(arg) + ' = ' + str(getattr(previous_run_data['args'],arg)) + '\n\tnew argument: ' + str(arg) + ' = ' + str(getattr(args,arg)) + '\nRerunning.')
                            args_are_same = False

                    if args_are_same:
                        if 'end_time_string' in previous_run_data:
                            info('Analysis already completed on %s!'%previous_run_data['end_time_string'])
                            sys.exit(0)
                        else:
                            can_finish_incomplete_run = True
                            #add previous run info to this run
                            if 'finished_steps' in previous_run_data:
                                for key in previous_run_data['finished_steps'].keys():
                                    crispresso2_info['finished_steps'][key] = previous_run_data['finished_steps'][key]
                                    if args.debug:
                                        info('finished: ' + key)
                else:
                    info('The no_rerun flag is set, but this analysis will be rerun because the existing run was performed using an old version of CRISPResso (' + str(previous_run_data['version']) + ').')

        #write this file early on so we can check the params if we have to rerun
        with open(crispresso2_info_file,"wb") as info_file:
            cp.dump(crispresso2_info, info_file )

        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        info('Processing input')

        #read filtering (for quality) is done at the individual crispresso run
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
                 cmd='%s SE -phred33 %s %s %s >>%s 2>&1'\
                 % (args.trimmomatic_command,args.fastq_r1,
                    output_forward_filename,
                    args.trimmomatic_options_string,
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
                 cmd='%s PE -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
                 % (args.trimmomatic_command,
                         args.fastq_r1,args.fastq_r2,output_forward_paired_filename,
                         output_forward_unpaired_filename,output_reverse_paired_filename,
                         output_reverse_unpaired_filename,args.trimmomatic_options_string,log_filename)
                 #print cmd
                 TRIMMOMATIC_STATUS=sb.call(cmd,shell=True)
                 if TRIMMOMATIC_STATUS:
                         raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

                 info('Done!')


             max_overlap_string = ""
             min_overlap_string = ""
             if args.max_paired_end_reads_overlap:
                 max_overlap_string = "--max-overlap " + str(args.max_paired_end_reads_overlap)
             if args.min_paired_end_reads_overlap:
                 min_overlap_string = args.min_paired_end_reads_overlap
             #Merging with Flash
             info('Merging paired sequences with Flash...')
             cmd=args.flash_command+' --allow-outies %s %s %s %s -z -d %s >>%s 2>&1' %\
             (output_forward_paired_filename,
              output_reverse_paired_filename,
              max_overlap_string,
              max_overlap_string,
              OUTPUT_DIRECTORY,log_filename)

             FLASH_STATUS=sb.call(cmd,shell=True)
             if FLASH_STATUS:
                 raise FlashException('Flash failed to run, please check the log file.')

             flash_hist_filename=_jp('out.hist')
             flash_histogram_filename=_jp('out.histogram')
             flash_not_combined_1_filename=_jp('out.notCombined_1.fastq.gz')
             flash_not_combined_2_filename=_jp('out.notCombined_2.fastq.gz')

             processed_output_filename=_jp('out.extendedFrags.fastq.gz')

             if args.force_merge_pairs:
                 old_flashed_filename = processed_output_filename
                 new_merged_filename=_jp('out.forcemerged_uncombined.fastq.gz')
                 num_reads_force_merged = CRISPRessoShared.force_merge_pairs(flash_not_combined_1_filename,flash_not_combined_2_filename,new_merged_filename)
                 new_output_filename=_jp('out.forcemerged.fastq.gz')
                 merge_command = "cat %s %s > %s"%(processed_output_filename,new_merged_filename,new_output_filename)
                 MERGE_STATUS=sb.call(merge_command,shell=True)
                 if MERGE_STATUS:
                     raise FlashException('Force-merging read pairs failed to run, please check the log file.')
                 processed_output_filename = new_output_filename

             info('Done!')


        if can_finish_incomplete_run and 'count_input_reads' in crispresso2_info['finished_steps']:
            (N_READS_INPUT,N_READS_AFTER_PREPROCESSING) = crispresso2_info['finished_steps']['count_input_reads']
        #count reads
        else:
            N_READS_INPUT=get_n_reads_fastq(args.fastq_r1)
            N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)
            crispresso2_info['finished_steps']['count_input_reads'] = (N_READS_INPUT,N_READS_AFTER_PREPROCESSING)
            with open(crispresso2_info_file,"wb") as info_file:
                cp.dump(crispresso2_info, info_file)

        #load gene annotation
        if args.gene_annotations:
            info('Loading gene coordinates from annotation file: %s...' % args.gene_annotations)
            try:
                df_genes=pd.read_csv(args.gene_annotations,compression='gzip',sep="\t")
                df_genes.txEnd=df_genes.txEnd.astype(int)
                df_genes.txStart=df_genes.txStart.astype(int)
                df_genes.head()
            except:
               info('Failed to load the gene annotations file.')


        if RUNNING_MODE=='ONLY_AMPLICONS' or  RUNNING_MODE=='AMPLICONS_AND_GENOME':

            #load and validate template file
            df_template=pd.read_csv(args.amplicons_file,names=[
                    'Name','Amplicon_Sequence','sgRNA',
                    'Expected_HDR','Coding_sequence'],comment='#',sep='\t',dtype={'Name':str})

            if str(df_template.iloc[0,1]).lower() == "amplicon_sequence":
                df_template.drop(0,axis=0,inplace=True)
                info('Detected header in amplicon file.')


            #remove empty amplicons/lines
            df_template.dropna(subset=['Amplicon_Sequence'],inplace=True)
            df_template.dropna(subset=['Name'],inplace=True)

            df_template.Amplicon_Sequence=df_template.Amplicon_Sequence.apply(CRISPRessoShared.capitalize_sequence)
            df_template.Expected_HDR=df_template.Expected_HDR.apply(CRISPRessoShared.capitalize_sequence)
            df_template.sgRNA=df_template.sgRNA.apply(CRISPRessoShared.capitalize_sequence)
            df_template.Coding_sequence=df_template.Coding_sequence.apply(CRISPRessoShared.capitalize_sequence)

            if not len(df_template.Amplicon_Sequence.unique())==df_template.shape[0]:
                duplicated_entries = df_template.Amplicon_Sequence[df_template.Amplicon_Sequence.duplicated()]
                raise Exception('The amplicon sequences must be distinct! (Duplicated entries: ' + str(duplicated_entries.values) + ')')

            if not len(df_template.Name.unique())==df_template.shape[0]:
                duplicated_entries = df_template.Name[df_template.Name.duplicated()]
                raise Exception('The amplicon names must be distinct! (Duplicated names: ' + str(duplicated_entries.values) + ')')

            df_template=df_template.set_index('Name')
            df_template.index=df_template.index.to_series().str.replace(' ','_')

            for idx,row in df_template.iterrows():

                wrong_nt=CRISPRessoShared.find_wrong_nt(row.Amplicon_Sequence)
                if wrong_nt:
                     raise NTException('The amplicon sequence %s contains wrong characters:%s' % (idx,' '.join(wrong_nt)))

                if not pd.isnull(row.sgRNA):

                    cut_points=[]

                    for current_guide_seq in row.sgRNA.strip().upper().split(','):

                        wrong_nt=CRISPRessoShared.find_wrong_nt(current_guide_seq)
                        if wrong_nt:
                            raise NTException('The sgRNA sequence %s contains wrong characters:%s'  % (current_guide_seq, ' '.join(wrong_nt)))

                        offset_fw=args.quantification_window_center+len(current_guide_seq)-1
                        offset_rc=(-args.quantification_window_center)-1
                        cut_points+=[m.start() + offset_fw for \
                                    m in re.finditer(current_guide_seq,  row.Amplicon_Sequence)]+[m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse_complement(current_guide_seq),  row.Amplicon_Sequence)]

                    if not cut_points:
                        warn('\nThe guide sequence/s provided: %s is(are) not present in the amplicon sequence:%s! \nNOTE: The guide will be ignored for the analysis. Please check your input!' % (row.sgRNA,row.Amplicon_Sequence))
                        df_template.ix[idx,'sgRNA']=''



        if RUNNING_MODE=='ONLY_AMPLICONS':
            #create a fasta file with all the amplicons
            amplicon_fa_filename=_jp('AMPLICONS.fa')
            fastq_gz_amplicon_filenames=[]
            with open(amplicon_fa_filename,'w+') as outfile:
                for idx,row in df_template.iterrows():
                    if row['Amplicon_Sequence']:
                        outfile.write('>%s\n%s\n' %(clean_filename('AMPL_'+idx),row['Amplicon_Sequence']))

                        #create place-holder fastq files
                        fastq_gz_amplicon_filenames.append(_jp('%s.fastq.gz' % clean_filename('AMPL_'+idx)))
                        open(fastq_gz_amplicon_filenames[-1], 'w+').close()

            df_template['Demultiplexed_fastq.gz_filename']=fastq_gz_amplicon_filenames
            info('Creating a custom index file with all the amplicons...')
            custom_index_filename=_jp('CUSTOM_BOWTIE2_INDEX')
            sb.call('bowtie2-build %s %s >>%s 2>&1' %(amplicon_fa_filename,custom_index_filename,log_filename), shell=True)


            #align the file to the amplicons (MODE 1)
            info('Align reads to the amplicons...')
            bam_filename_amplicons= _jp('CRISPResso_AMPLICONS_ALIGNED.bam')
            aligner_command= 'bowtie2 -x %s -p %s %s -U %s 2>>%s | samtools view -bS - > %s' %(custom_index_filename,n_processes,args.bowtie2_options_string,processed_output_filename,log_filename,bam_filename_amplicons)


            info('Alignment command: ' + aligner_command)
            sb.call(aligner_command,shell=True)

            N_READS_ALIGNED=get_n_aligned_bam(bam_filename_amplicons)

            s1=r"samtools view -F 4 %s 2>>%s | grep -v ^'@'" % (bam_filename_amplicons,log_filename)
            s2=r'''|awk '{ gzip_filename=sprintf("gzip >> OUTPUTPATH%s.fastq.gz",$3);\
            print "@"$1"\n"$10"\n+\n"$11  | gzip_filename;}' '''

            cmd=s1+s2.replace('OUTPUTPATH',_jp(''))
            sb.call(cmd,shell=True)

            info('Demultiplex reads and run CRISPResso on each amplicon...')
            n_reads_aligned_amplicons=[]
            crispresso_cmds = []
            for idx,row in df_template.iterrows():
                info('\n Processing:%s' %idx)
                n_reads_aligned_amplicons.append(get_n_reads_fastq(row['Demultiplexed_fastq.gz_filename']))
                crispresso_cmd= args.crispresso_command + ' -r1 %s -a %s -o %s --name %s' % (row['Demultiplexed_fastq.gz_filename'],row['Amplicon_Sequence'],OUTPUT_DIRECTORY,idx)

                if n_reads_aligned_amplicons[-1]>args.min_reads_to_use_region:
                    if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                        crispresso_cmd+=' -g %s' % row['sgRNA']

                    if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                        crispresso_cmd+=' -e %s' % row['Expected_HDR']

                    if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                        crispresso_cmd+=' -c %s' % row['Coding_sequence']

                    crispresso_cmd=CRISPRessoShared.propagate_crispresso_options(crispresso_cmd,crispresso_options_for_pooled,args)
                    crispresso_cmds.append(crispresso_cmd)

                else:
                    warn('Skipping amplicon [%s] because no reads align to it\n'% idx)

            CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,n_processes,'amplicon',args.skip_failed)

            df_template['n_reads']=n_reads_aligned_amplicons
            df_template['n_reads_aligned_%']=df_template['n_reads']/float(N_READS_ALIGNED)*100
            df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_AMPLICONS.txt'),sep='\t')



        if RUNNING_MODE=='AMPLICONS_AND_GENOME':
            info('Mapping amplicons to the reference genome...')

            filename_amplicon_aligned_locations = _jp('CRISPResso_amplicon_aligned_locations.csv')
            filename_aligned_amplicons_sam = _jp('CRISPResso_amplicons_aligned.sam')
            filename_aligned_amplicons_sam_log = _jp('CRISPResso_amplicons_aligned.sam.log')
            filename_amplicon_seqs_fasta = _jp('CRISPResso_amplicons_to_align.fa')

            if can_finish_incomplete_run and 'mapping_amplicons_to_reference_genome' in crispresso2_info['finished_steps']:
                info('Reading previously-computed alignment of amplicons to genome')
                additional_columns_df = pd.read_csv(filename_amplicon_aligned_locations,sep="\t")
                additional_columns_df.set_index('Name',inplace=True)
            else:
                #write amplicons as fastq for alignment
                with open(filename_amplicon_seqs_fasta,'w') as fastas:
                    for idx,row in df_template.iterrows():
                        fastas.write('>%s\n%s\n'%(idx,row.Amplicon_Sequence))

                aligner_command= 'bowtie2 -x %s -p %s %s -f -U %s --no-hd --no-sq 2> %s > %s ' %(args.bowtie2_index,n_processes,args.bowtie2_options_string, \
                    filename_amplicon_seqs_fasta,filename_aligned_amplicons_sam_log,filename_aligned_amplicons_sam)
                bowtie_status=sb.call(aligner_command,shell=True)
                if bowtie_status:
                        raise Bowtie2Exception('Bowtie2 failed to align amplicons to the genome, please check the output file.')

                additional_columns = []
                with open (filename_aligned_amplicons_sam) as aln:
                    for line in aln.readlines():
                        line_els = line.split("\t")
                        if line_els[2] == "*":
                            info('The amplicon [%s] is not mappable to the reference genome provided!' % idx )
                            additional_columns.append([line_els[0],'NOT_ALIGNED',0,-1,'+',''])
                        else:
                            aln_len = get_ref_length_from_cigar(line_els[5])
                            seq_start = int(line_els[3])
                            seq_stop = seq_start + aln_len
                            strand = "-" if (int(line_els[1]) & 0x10) else "+"
                            additional_columns.append([line_els[0],line_els[2],seq_start,seq_stop,strand,line_els[9]])
                            info('The amplicon [%s] was mapped to: %s:%d-%d ' % (line_els[0],line_els[2],seq_start,seq_stop))
                additional_columns_df = pd.DataFrame(additional_columns,columns=['Name','chr_id','bpstart','bpend','strand','Reference_Sequence']).set_index('Name')
                additional_columns_df.to_csv(filename_amplicon_aligned_locations,sep="\t",index_label='Name')

                crispresso2_info['finished_steps']['mapping_amplicons_to_reference_genome'] = True
                with open(crispresso2_info_file,"wb") as info_file:
                    cp.dump(crispresso2_info, info_file)

            files_to_remove.append(filename_amplicon_seqs_fasta)
            files_to_remove.append(filename_aligned_amplicons_sam)

            df_template=df_template.join(additional_columns_df)

            df_template.bpstart=df_template.bpstart.astype(int)
            df_template.bpend=df_template.bpend.astype(int)

            #Check reference is the same otherwise throw a warning
            for idx,row in df_template.iterrows():
                if row.Amplicon_Sequence != row.Reference_Sequence and row.Amplicon_Sequence != CRISPRessoShared.reverse_complement(row.Reference_Sequence):
                    warn('The amplicon sequence %s provided:\n%s\n\nis different from the reference sequence(both strands):\n\n%s\n\n%s\n' %(row.name,row.Amplicon_Sequence,row.Amplicon_Sequence,CRISPRessoShared.reverse_complement(row.Amplicon_Sequence)))


        if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':

            ###HERE we recreate the uncompressed genome file if not available###

            #check you have all the files for the genome and create a fa idx for samtools

            uncompressed_reference=args.bowtie2_index+'.fa'

            #if not os.path.exists(GENOME_LOCAL_FOLDER):
            #    os.mkdir(GENOME_LOCAL_FOLDER)

            if os.path.exists(uncompressed_reference):
                info('The uncompressed reference fasta file for %s is already present! Skipping generation.' % args.bowtie2_index)
            else:
                #uncompressed_reference=os.path.join(GENOME_LOCAL_FOLDER,'UNCOMPRESSED_REFERENCE_FROM_'+args.bowtie2_index.replace('/','_')+'.fa')
                info('Extracting uncompressed reference from the provided bowtie2 index since it is not available... Please be patient!')

                cmd_to_uncompress='bowtie2-inspect %s > %s 2>>%s' % (args.bowtie2_index,uncompressed_reference,log_filename)
                sb.call(cmd_to_uncompress,shell=True)

                info('Indexing fasta file with samtools...')
                #!samtools faidx {uncompressed_reference}
                sb.call('samtools faidx %s 2>>%s ' % (uncompressed_reference,log_filename),shell=True)


        #align in unbiased way the reads to the genome
        if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
            bam_filename_genome = _jp('%s_GENOME_ALIGNED.bam' % database_id)

            if can_finish_incomplete_run and 'n_reads_aligned_genome' in crispresso2_info['finished_steps']:
                info('Using previously-computed alignment of reads to genome')
                N_READS_ALIGNED = crispresso2_info['finished_steps']['n_reads_aligned_genome']
            else:
                info('Aligning reads to the provided genome index...')
                aligner_command= 'bowtie2 -x %s -p %s %s -U %s 2>>%s| samtools view -bS - | samtools sort -@ %d - -o %s' %(args.bowtie2_index,n_processes,
                    args.bowtie2_options_string,processed_output_filename,log_filename,n_processes,bam_filename_genome)
                info('aligning with command: ' + aligner_command)
                sb.call(aligner_command,shell=True)

                sb.call('samtools index %s' % bam_filename_genome,shell=True)

                N_READS_ALIGNED=get_n_aligned_bam(bam_filename_genome)

                #save progress up to this point
                crispresso2_info['finished_steps']['n_reads_aligned_genome'] = N_READS_ALIGNED
                with open(crispresso2_info_file,"wb") as info_file:
                    cp.dump(crispresso2_info, info_file )


            MAPPED_REGIONS=_jp('MAPPED_REGIONS/')
            if not os.path.exists(MAPPED_REGIONS):
                os.mkdir(MAPPED_REGIONS)

            if can_finish_incomplete_run and 'genome_demultiplexing' in crispresso2_info['finished_steps']:
                info('Using previously-computed demultiplexing of genomic reads')
            else:
                #REDISCOVER LOCATIONS and DEMULTIPLEX READS

                cmd = "samtools view -H %s" % bam_filename_genome
                p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
                chr_lines = p.communicate()[0].split("\n")
                chrs = []
                for chr_line in chr_lines:
                    m = re.match(r'@SQ\s+SN:(\S+)\s+LN:',chr_line)
                    if m:
                        chrs.append(m.group(1))

                s1=r'''samtools view -F 0x0004 %s __CHR__ 2>>%s |''' % (bam_filename_genome,log_filename)+\
                r'''awk '{OFS="\t"; bpstart=$4;  bpend=bpstart; split ($6,a,"[MIDNSHP]"); n=0;\
                for (i=1; i in a; i++){\
                    n+=1+length(a[i]);\
                    if (substr($6,n,1)=="S"){\
                        if (bpend==$4)\
                            bpstart-=a[i];\
                        else
                            bpend+=a[i];
                        }\
                    else if( (substr($6,n,1)!="I")  && (substr($6,n,1)!="H") )\
                            bpend+=a[i];\
                    }\
                    if ( ($2 % 32)>=16)\
                        print $3,bpstart,bpend,"-",$1,$10,$11;\
                    else\
                        print $3,bpstart,bpend,"+",$1,$10,$11;}' | '''

                s2=r'''  sort -k1,1 -k2,2n  | awk \
                'BEGIN{chr_id="NA";bpstart=-1;bpend=-1; fastq_filename="NA"}\
                { if ( (chr_id!=$1) || (bpstart!=$2) || (bpend!=$3) )\
                    {\
                    if (fastq_filename!="NA") {close(fastq_filename); system("gzip -f "fastq_filename)}\
                    chr_id=$1; bpstart=$2; bpend=$3;\
                    fastq_filename=sprintf("__OUTPUTPATH__REGION_%s_%s_%s.fastq",$1,$2,$3);\
                    }\
                print "@"$5"\n"$6"\n+\n"$7 >> fastq_filename;\
                }' '''
                s3=" && gzip -f __OUTPUTPATH__/REGION___CHR__*.fastq ";
                cmd=(s1+s2+s3).replace('__OUTPUTPATH__',MAPPED_REGIONS)

                chr_commands = []
                for chr in chrs:
                    chr_cmd=cmd.replace('__CHR__',chr)
                    chr_commands.append(chr_cmd)

                info('Demultiplexing reads by location...')
                CRISPRessoMultiProcessing.run_parallel_commands(chr_commands,n_processes=n_processes,descriptor='Demultiplexing reads by location',continue_on_fail=args.skip_failed)

                #todo: sometimes no reads align -- we should alert the user and display an error here

                #gzip the missing ones
                #sb.call('gzip -qf %s/*.fastq' % MAPPED_REGIONS,shell=True)

                crispresso2_info['finished_steps']['genome_demultiplexing'] = True
                with open(crispresso2_info_file,"wb") as info_file:
                    cp.dump(crispresso2_info, info_file)

        '''
        The most common use case, where many different target sites are pooled into a single
        high-throughput sequencing library for quantification, is not directly addressed by this implementation.
        Potential users of CRISPResso would need to write their own code to generate separate input files for processing.
        Importantly, this preprocessing code would need to remove any PCR amplification artifacts
        (such as amplification of sequences from a gene and a highly similar pseudogene )
        which may confound the interpretation of results.
        This can be done by mapping of input sequences to a reference genome and removing
        those that do not map to the expected genomic location, but is non-trivial for an end-user to implement.
        '''


        if RUNNING_MODE=='AMPLICONS_AND_GENOME':
            files_to_match=glob.glob(os.path.join(MAPPED_REGIONS,'REGION*'))
            n_reads_aligned_genome=[]
            fastq_region_filenames=[]

            if can_finish_incomplete_run and 'crispresso_amplicons_and_genome' in crispresso2_info['finished_steps']:
                info('Using previously-computed crispresso runs')
                (n_reads_aligned_genome,fastq_region_filenames,files_to_match) = crispresso2_info['finished_steps']['crispresso_amplicons_and_genome'];
            else:
                crispresso_cmds = []
                for idx,row in df_template.iterrows():

                    info('Processing amplicon: %s' % idx )

                    #check if we have reads
                    fastq_filename_region=os.path.join(MAPPED_REGIONS,'REGION_%s_%s_%s.fastq.gz' % (row['chr_id'],row['bpstart'],row['bpend']))

                    if os.path.exists(fastq_filename_region):

                        N_READS=get_n_reads_fastq(fastq_filename_region)
                        n_reads_aligned_genome.append(N_READS)
                        fastq_region_filenames.append(fastq_filename_region)
                        if fastq_filename_region in files_to_match:
                            files_to_match.remove(fastq_filename_region)
                        #else:
                             #info('Warning: Fastq filename ' + fastq_filename_region + ' is not in ' + str(files_to_match))
                             #debug here??
                        if N_READS>=args.min_reads_to_use_region:
                            info('\nThe amplicon [%s] has enough reads (%d) mapped to it! Running CRISPResso!\n' % (idx,N_READS))

                            crispresso_cmd= args.crispresso_command + ' -r1 %s -a %s -o %s --name %s' % (fastq_filename_region,row['Amplicon_Sequence'],OUTPUT_DIRECTORY,idx)

                            if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                                crispresso_cmd+=' -g %s' % row['sgRNA']

                            if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                                crispresso_cmd+=' -e %s' % row['Expected_HDR']

                            if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                                crispresso_cmd+=' -c %s' % row['Coding_sequence']

                            crispresso_cmd=CRISPRessoShared.propagate_crispresso_options(crispresso_cmd,crispresso_options_for_pooled,args)
                            info('Running CRISPResso:%s' % crispresso_cmd)
                            crispresso_cmds.append(crispresso_cmd)

                        else:
                             warn('The amplicon [%s] has not enough reads (%d) mapped to it! Skipping the execution of CRISPResso!' % (idx,N_READS))
                    else:
                        fastq_region_filenames.append('')
                        n_reads_aligned_genome.append(0)
                        warn("The amplicon %s doesn't have any reads mapped to it!\n Please check your amplicon sequence." %  idx)

                CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,n_processes,'amplicon',args.skip_failed)

                crispresso2_info['finished_steps']['crispresso_amplicons_and_genome'] = (n_reads_aligned_genome,fastq_region_filenames,files_to_match)
                with open(crispresso2_info_file,"wb") as info_file:
                    cp.dump(crispresso2_info, info_file)

            df_template['Amplicon_Specific_fastq.gz_filename']=fastq_region_filenames
            df_template['n_reads']=n_reads_aligned_genome
            df_template['n_reads_aligned_%']=df_template['n_reads']/float(N_READS_ALIGNED)*100

            if args.gene_annotations:
                df_template=df_template.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)

            df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_GENOME_AND_AMPLICONS.txt'),sep='\t')

            #write another file with the not amplicon regions

            if args.skip_reporting_problematic_regions:
                df_regions=pd.DataFrame(columns=['chr_id','bpstart','bpend','fastq_file','n_reads','Reference_sequence'])
            else:
                filename_problematic_regions = _jp('REPORTS_READS_ALIGNED_TO_GENOME_NOT_MATCHING_AMPLICONS.txt')
                if can_finish_incomplete_run and 'reporting_problematic_regions' in crispresso2_info['finished_steps']:
                    info('Skipping previously-computed reporting of problematic regions')
                    df_regions = pd.read_csv(filename_problematic_regions,sep='\t')
                else:
                    info('Reporting problematic regions...')
                    summarize_region_fastq_input = [f+" "+uncompressed_reference for f in files_to_match] #pass both params to parallel function
                    coordinates = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(summarize_region_fastq_input,summarize_region_fastq_chunk,n_processes=n_processes)
                    df_regions=pd.DataFrame(coordinates,columns=['chr_id','bpstart','bpend','fastq_file','n_reads','Reference_sequence'])
                    df_regions.dropna(inplace=True) #remove regions in chrUn

                    df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
                    df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])
                    df_regions['n_reads'] = pd.to_numeric(df_regions['n_reads'])

                    df_regions.bpstart=df_regions.bpstart.astype(int)
                    df_regions.bpend=df_regions.bpend.astype(int)

                    df_regions['n_reads_aligned_%']=df_regions['n_reads']/float(N_READS_ALIGNED)*100

                    if args.gene_annotations:
                        info('Checking overlapping genes...')
                        df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)

                    if np.sum(np.array(map(int,pd.__version__.split('.')))*(100,10,1))< 170:
                        df_regions.sort('n_reads',ascending=False,inplace=True)
                    else:
                        df_regions.sort_values(by='n_reads',ascending=False,inplace=True)

                    df_regions.fillna('NA').to_csv(filename_problematic_regions,sep='\t',index=None)

                    crispresso2_info['finished_steps']['reporting_problematic_regions'] = True
                    with open(crispresso2_info_file,"wb") as info_file:
                        cp.dump(crispresso2_info, info_file)


        if RUNNING_MODE=='ONLY_GENOME' :
            #Load regions and build REFERENCE TABLES
            filename_reads_aligned_to_genome_only = _jp('REPORT_READS_ALIGNED_TO_GENOME_ONLY.txt')
            if can_finish_incomplete_run and 'demultiplexing_genome_only_regions' in crispresso2_info['finished_steps']:
                info('Using previously-computed extraction of aligned regions')
                df_regions = pd.read_csv(filename_reads_aligned_to_genome_only,sep="\t")
            else:
                info('Parsing the demultiplexed files and extracting locations and reference sequences...')
                files_to_match = glob.glob(os.path.join(MAPPED_REGIONS,'REGION*.fastq.gz'))
                summarize_region_fastq_input = [f+" "+uncompressed_reference for f in files_to_match] #pass both params to parallel function
                coordinates = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(summarize_region_fastq_input,summarize_region_fastq_chunk,n_processes=n_processes)
                df_regions=pd.DataFrame(coordinates,columns=['chr_id','bpstart','bpend','fastq_file','n_reads','sequence'])

                df_regions.dropna(inplace=True) #remove regions in chrUn

                df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
                df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])
                df_regions['n_reads'] = pd.to_numeric(df_regions['n_reads'])

                df_regions.bpstart=df_regions.bpstart.astype(int)
                df_regions.bpend=df_regions.bpend.astype(int)

                df_regions['n_reads_aligned_%']=df_regions['n_reads']/float(N_READS_ALIGNED)*100

                if args.gene_annotations:
                    info('Checking overlapping genes...')
                    df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)

                if np.sum(np.array(map(int,pd.__version__.split('.')))*(100,10,1))< 170:
                    df_regions.sort('n_reads',ascending=False,inplace=True)
                else:
                    df_regions.sort_values(by='n_reads',ascending=False,inplace=True)

                df_regions.fillna('NA').to_csv(filename_reads_aligned_to_genome_only,sep='\t',index=None)

                crispresso2_info['finished_steps']['demultiplexing_genome_only_regions'] = True
                with open(crispresso2_info_file,"wb") as info_file:
                    cp.dump(crispresso2_info, info_file)


            #run CRISPResso
            #demultiplex reads in the amplicons and call crispresso!
            if can_finish_incomplete_run and 'crispresso_genome_only' in crispresso2_info['finished_steps']:
                info('Using previously-computed crispresso runs')
            else:
                info('Running CRISPResso on the regions discovered...')
                crispresso_cmds = []
                for idx,row in df_regions.iterrows():

                    if row.n_reads > args.min_reads_to_use_region:
                        info('\nRunning CRISPResso on: %s-%d-%d...'%(row.chr_id,row.bpstart,row.bpend ))
                        crispresso_cmd= args.crispresso_command + ' -r1 %s -a %s -o %s' %(row.fastq_file,row.sequence,OUTPUT_DIRECTORY)
                        crispresso_cmd=CRISPRessoShared.propagate_crispresso_options(crispresso_cmd,crispresso_options_for_pooled,args)
                        crispresso_cmds.append(crispresso_cmd)
                    else:
                        info('Skipping region: %s-%d-%d , not enough reads (%d)' %(row.chr_id,row.bpstart,row.bpend, row.n_reads))
                CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,n_processes,'region',args.skip_failed)

                crispresso2_info['finished_steps']['crispresso_genome_only'] = True
                with open(crispresso2_info_file,"wb") as info_file:
                    cp.dump(crispresso2_info, info_file)

        #write alignment statistics
        with open(_jp('MAPPING_STATISTICS.txt'),'w+') as outfile:
            outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d' % (N_READS_INPUT,N_READS_AFTER_PREPROCESSING,N_READS_ALIGNED))

        quantification_summary=[]

        if RUNNING_MODE=='ONLY_AMPLICONS' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
            df_final_data=df_template
        else:
            df_final_data=df_regions

        all_region_names = []
        all_region_read_counts = {}
        good_region_names = []
        good_region_folders = {}
        header = 'Name\tUnmodified%\tModified%\tReads_total\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions'
        header_els = header.split("\t")
        header_el_count = len(header_els)
        empty_line_els = [np.nan]*(header_el_count-1)
        n_reads_index = header_els.index('Reads_total') - 1
        for idx,row in df_final_data.iterrows():
                run_name = idx
                if RUNNING_MODE=='ONLY_AMPLICONS' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
                    run_name=idx
                else:
                    run_name='REGION_%s_%d_%d' %(row.chr_id,row.bpstart,row.bpend )
                folder_name = 'CRISPResso_on_%s'%run_name

                all_region_names.append(run_name)
                all_region_read_counts[run_name] = row.n_reads

                run_file = os.path.join(_jp(folder_name),'CRISPResso2_info.pickle')
                if not os.path.exists(run_file):
                    warn('Skipping the folder %s: not enough reads, incomplete, or empty folder.'% folder_name)
                    this_els = empty_line_els[:]
                    this_els[n_reads_index] = row.n_reads
                    to_add = [run_name]
                    to_add.extend(this_els)
                    quantification_summary.append(to_add)
                else:
                    run_data = cp.load(open(run_file,'rb'))
                    ref_name = run_data['ref_names'][0] #only expect one amplicon sequence
                    n_tot = row.n_reads
                    n_aligned = run_data['counts_total'][ref_name]
                    n_unmod = run_data['counts_unmodified'][ref_name]
                    n_mod = run_data['counts_modified'][ref_name]
                    n_discarded = run_data['counts_discarded'][ref_name]

                    n_insertion = run_data['counts_insertion'][ref_name]
                    n_deletion = run_data['counts_deletion'][ref_name]
                    n_substitution = run_data['counts_substitution'][ref_name]
                    n_only_insertion = run_data['counts_only_insertion'][ref_name]
                    n_only_deletion = run_data['counts_only_deletion'][ref_name]
                    n_only_substitution = run_data['counts_only_substitution'][ref_name]
                    n_insertion_and_deletion = run_data['counts_insertion_and_deletion'][ref_name]
                    n_insertion_and_substitution = run_data['counts_insertion_and_substitution'][ref_name]
                    n_deletion_and_substitution = run_data['counts_deletion_and_substitution'][ref_name]
                    n_insertion_and_deletion_and_substitution = run_data['counts_insertion_and_deletion_and_substitution'][ref_name]

                    unmod_pct = np.nan
                    mod_pct = np.nan
                    if n_aligned > 0:
                        unmod_pct = 100*n_unmod/float(n_aligned)
                        mod_pct = 100*n_mod/float(n_aligned)


                    vals = [run_name]
                    vals.extend([round(unmod_pct,8),round(mod_pct,8),n_aligned,n_tot,n_unmod,n_mod,n_discarded,n_insertion,n_deletion,n_substitution,n_only_insertion,n_only_deletion,n_only_substitution,n_insertion_and_deletion,n_insertion_and_substitution,n_deletion_and_substitution,n_insertion_and_deletion_and_substitution])
                    quantification_summary.append(vals)

                    good_region_names.append(run_name)
                    good_region_folders[idx] = folder_name


        samples_quantification_summary_filename = _jp('SAMPLES_QUANTIFICATION_SUMMARY.txt')

        df_summary_quantification=pd.DataFrame(quantification_summary,columns=header_els)
        if args.crispresso1_mode:
            crispresso1_columns=['Name','Unmodified%','Modified%','Reads_aligned','Reads_total']
            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename,sep='\t',index=None,columns=crispresso1_columns)
        else:

            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename,sep='\t',index=None)

        crispresso2_info['samples_quantification_summary_filename'] = os.path.basename(samples_quantification_summary_filename)
        crispresso2_info['final_data'] = df_final_data
        crispresso2_info['all_region_names'] = all_region_names
        crispresso2_info['all_region_read_counts'] = all_region_read_counts
        crispresso2_info['good_region_names'] = good_region_names
        crispresso2_info['good_region_folders'] = good_region_folders
        crispresso2_info['running_mode'] = RUNNING_MODE

        crispresso2_info['summary_plot_names'] = []
        crispresso2_info['summary_plot_titles'] = {}
        crispresso2_info['summary_plot_labels'] = {}
        crispresso2_info['summary_plot_datas'] = {}


        df_summary_quantification.set_index('Name')

        save_png = True
        if args.suppress_report:
            save_png = False

        if not args.suppress_plots:
            plot_root = _jp("CRISPRessoPooled_reads_summary")

            CRISPRessoPlot.plot_reads_total(plot_root,df_summary_quantification,save_png,args.min_reads_to_use_region)
            plot_name = os.path.basename(plot_root)
            crispresso2_info['summary_plot_root'] = plot_name
            crispresso2_info['summary_plot_names'].append(plot_name)
            crispresso2_info['summary_plot_titles'][plot_name] = 'CRISPRessoPooled Read Allocation Summary'
            crispresso2_info['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads allocated to each amplicon. The vertical line shows the cutoff for analysis, set using the --min_reads_to_use_region parameter.'
            crispresso2_info['summary_plot_datas'][plot_name] = [('CRISPRessoPooled summary',os.path.basename(samples_quantification_summary_filename))]

            plot_root = _jp("CRISPRessoPooled_modification_summary")
            CRISPRessoPlot.plot_unmod_mod_pcts(plot_root,df_summary_quantification,save_png,args.min_reads_to_use_region)
            plot_name = os.path.basename(plot_root)
            crispresso2_info['summary_plot_root'] = plot_name
            crispresso2_info['summary_plot_names'].append(plot_name)
            crispresso2_info['summary_plot_titles'][plot_name] = 'CRISPRessoPooled Modification Summary'
            crispresso2_info['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads aligned to each amplicon, divided into the reads that are modified and unmodified. The vertical line shows the cutoff for analysis, set using the --min_reads_to_use_region parameter.'
            crispresso2_info['summary_plot_datas'][plot_name] = [('CRISPRessoPooled summary',os.path.basename(samples_quantification_summary_filename))]




        #if many reads weren't aligned, print those out for the user
        if RUNNING_MODE != 'ONLY_GENOME':
            #N_READS_INPUT=get_n_reads_fastq(args.fastq_r1)
            #N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)
    		tot_reads_aligned = df_summary_quantification['Reads_aligned'].fillna(0).sum()
    		tot_reads = df_summary_quantification['Reads_total'].sum()

    		if RUNNING_MODE=='AMPLICONS_AND_GENOME':
    			this_bam_filename = bam_filename_genome
    		if RUNNING_MODE=='ONLY_AMPLICONS':
    			this_bam_filename = bam_filename_amplicons
    		#if less than 1/2 of reads aligned, find most common unaligned reads and advise the user
    		if N_READS_INPUT > 0 and tot_reads/float(N_READS_INPUT) < 0.5:
    			warn('Less than half (%d/%d) of reads aligned to amplicons. Finding most frequent unaligned reads.'%(tot_reads,N_READS_INPUT))
    			###
    			###this results in the unpretty messages being printed:
    			### sort: write failed: standard output: Broken pipe
    			### sort: write error
    			###
    			#cmd = "samtools view -f 4 %s | awk '{print $10}' | sort | uniq -c | sort -nr | head -n 10"%this_bam_filename
    			import signal
    			def default_sigpipe():
    				    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    			cmd = "samtools view -f 4 %s | head -n 10000 | awk '{print $10}' | sort | uniq -c | sort -nr | head -n 10 | awk '{print $2}'"%this_bam_filename
#    			print("command is: "+cmd)
#    		    p = sb.Popen(cmd, shell=True,stdout=sb.PIPE)
		    	p = sb.Popen(cmd, shell=True,stdout=sb.PIPE,preexec_fn=default_sigpipe)
    			top_unaligned = p.communicate()[0]
    			top_unaligned_filename=_jp('CRISPRessoPooled_TOP_UNALIGNED.txt')

    			with open(top_unaligned_filename,'w') as outfile:
    				outfile.write(top_unaligned)
    			warn('Perhaps one or more of the given amplicon sequences were incomplete or incorrect. Below is a list of the most frequent unaligned reads (in the first 10000 unaligned reads). Check this list to see if an amplicon is among these reads.\n%s'%top_unaligned)


        #cleaning up
        if not args.keep_intermediate:
             info('Removing Intermediate files...')

             if args.fastq_r2!='':
                 files_to_remove+=[processed_output_filename,flash_hist_filename,flash_histogram_filename,\
                              flash_not_combined_1_filename,flash_not_combined_2_filename]
                 if args.force_merge_pairs:
                    files_to_remove.append(new_merged_filename)
                    files_to_remove.append(old_flashed_filename)
             else:
                 files_to_remove+=[processed_output_filename]

             if args.trim_sequences and args.fastq_r2!='':
                 files_to_remove+=[output_forward_paired_filename,output_reverse_paired_filename,\
                                                   output_forward_unpaired_filename,output_reverse_unpaired_filename]

             if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
                     files_to_remove+=[bam_filename_genome]

             if RUNNING_MODE=='ONLY_AMPLICONS':
                files_to_remove+=[bam_filename_amplicons,amplicon_fa_filename]
                for bowtie2_file in glob.glob(_jp('CUSTOM_BOWTIE2_INDEX.*')):
                    files_to_remove.append(bowtie2_file)

             for file_to_remove in files_to_remove:
                 try:
                         if os.path.islink(file_to_remove):
                             #print 'LINK',file_to_remove
                             os.unlink(file_to_remove)
                         else:
                             os.remove(file_to_remove)
                 except:
                         warn('Skipping:%s' %file_to_remove)

        if not args.suppress_report and not args.suppress_plots:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Pooled_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_pooled_report_from_folder(report_name,crispresso2_info,OUTPUT_DIRECTORY,_ROOT)
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

        with open(crispresso2_info_file,"wb") as info_file:
            cp.dump(crispresso2_info, info_file )

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
