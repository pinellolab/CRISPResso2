# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''
import difflib
import os
import sys
from copy import deepcopy
from datetime import datetime
import subprocess as sb
import glob
import gzip
import unicodedata
import string
import re
import zipfile
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2.CRISPRessoReports import CRISPRessoReport
from CRISPResso2 import CRISPRessoPlot
import traceback

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

#get n_reads and region data from region fastq file (location is pulled from filename)
def summarize_region_fastq_chunk(input_arr):
    ret_val = []
    for input in input_arr:
#        print('doing region ' + str(input))
        region_fastq, uncompressed_reference = input.split(" ")
        #region format: REGION_chr8_1077_1198.fastq.gz
        #But if the chr has underscores, it could look like this:
        #    REGION_chr8_KI270812v1_alt_1077_1198.fastq.gz
        region_info = os.path.basename(region_fastq).replace('.fastq.gz', '').replace('.fastq', '').split('_')
        chr_string = "_".join(region_info[1:len(region_info)-2]) #in case there are underscores
        region_string='%s:%s-%d' % (chr_string, region_info[-2], int(region_info[-1])-1)
        p = sb.Popen("samtools faidx %s %s | grep -v ^\> | tr -d '\n'" %(uncompressed_reference, region_string), shell=True, stdout=sb.PIPE)
        seq = p.communicate()[0].decode('utf-8')
        p = sb.Popen(('z' if region_fastq.endswith('.gz') else '' ) +"cat < %s | wc -l" % region_fastq, shell=True, stdout=sb.PIPE)
        n_reads = int(float(p.communicate()[0])/4.0)
        ret_val.append([chr_string] + region_info[-2:]+[region_fastq, n_reads, seq])
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
        length, op = c
        if op in read_consuming_ops:
            result += int(length)
    return result

def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename, shell=True, stdout=sb.PIPE)
     n_reads = int(float(p.communicate()[0])/4.0)
     return n_reads

def get_n_reads_bam(bam_filename):
    p = sb.Popen("samtools view -c %s" % bam_filename, shell=True, stdout=sb.PIPE)
    return int(p.communicate()[0])

def get_n_aligned_bam(bam_filename):
    p = sb.Popen("samtools view -F 0x904 -c %s" % bam_filename, shell=True, stdout=sb.PIPE)
    return int(p.communicate()[0])

def get_n_aligned_bam_region(bam_filename, chr_name, chr_start, chr_end):
    p = sb.Popen("samtools view -F 0x904 -c %s %s:%d-%d" %(bam_filename, chr_name, chr_start, chr_end), shell=True, stdout=sb.PIPE)
    return int(p.communicate()[0])

def find_overlapping_genes(row, df_genes):
    df_genes_overlapping=df_genes.loc[(df_genes.chrom==row.chr_id) &
                                     (df_genes.txStart<=row.bpend) &
                                     (row.bpstart<=df_genes.txEnd)]
    genes_overlapping=[]

    for idx_g, row_g in df_genes_overlapping.iterrows():
        if 'name' in row_g.keys() and 'name2' in row_g.keys():
            genes_overlapping.append( '%s (%s)' % (row_g.name2, row_g['name']))
        elif '#name' in row_g.keys() and 'name2' in row_g.keys():
            genes_overlapping.append( '%s (%s)' % (row_g.name2, row_g['#name']))
        elif '#name' in row_g.keys():
            genes_overlapping.append( '%s' % (row_g['#name']))
        elif 'name' in row_g.keys():
            genes_overlapping.append( '%s' % (row_g['name']))
        else:
            genes_overlapping.append( '%s' % (row_g[0]))



    row['gene_overlapping']=','.join(genes_overlapping)

    return row


def normalize_name(name, fastq_r1, fastq_r2, aligned_pooled_bam):
    """Normalize the name according to the inputs and clean it.

    Parameters
    ----------
    name : str
        The name optionally given by the user.
    fastq_r1 : str
        The path to the first fastq file.
    fastq_r2 : str
        The path to the second fastq file.
    aligned_pooled_bam : str
        The path to the aligned bam file.

    Returns
    -------
    str
        The normalized name.
    """
    get_name_from_fastq = lambda x: os.path.basename(x).replace('.fastq', '').replace('.gz', '').replace('.fq', '')
    get_name_from_bam = lambda x: os.path.basename(x).replace('.bam', '')

    if not name:
        if aligned_pooled_bam is not None:
            return get_name_from_bam(aligned_pooled_bam)
        elif fastq_r2:
            return '%s_%s' % (get_name_from_fastq(fastq_r1), get_name_from_fastq(fastq_r2))
        else:
            return '%s' % get_name_from_fastq(fastq_r1)
    else:
        clean_name = CRISPRessoShared.slugify(name)
        if name != clean_name:
            warn(
                'The specified name {0} contained invalid characters and was changed to: {1}'.format(
                    name, clean_name,
                ),
            )
        return clean_name


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

        description = ['~~~CRISPRessoPooled~~~', '-Analysis of CRISPR/Cas9 outcomes from POOLED deep sequencing data-']
        pooled_string = r'''
 _______________________
| __  __  __     __ __  |
||__)/  \/  \|  |_ |  \ |
||   \__/\__/|__|__|__/ |
|_______________________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description, pooled_string))

        # if no args are given, print a simplified help message
        if len(sys.argv) == 1:
            print(CRISPRessoShared.format_cl_text('usage: CRISPRessoPooled [-r1 FASTQ_R1] [-r2 FASTQ_R2] [-f AMPLICONS_FILE] [-x GENOME_ROOT] [-n NAME]\n' + \
                'commonly-used arguments:\n' + \
                '-h, --help            show the full list of arguments\n' + \
                '-v, --version         show program\'s version number and exit\n' + \
                '-r1 FASTQ_R1          Input fastq file R1 (default: None)\n' + \
                '-r2 FASTQ_R2          Input fastq file R2 (default: None)\n' + \
                '-f AMPLICONS_FILE     Tab-separated file containing information for each amplicon, including at least columns for the amplicon_name and amplicon_seq (default: None)\n' + \
                '-x GENOME_ROOT        Folder that contains the bowtie2-indexed genome for optional unbiased alignment of reads (default: None, reads are only aligned to provided amplicon sequences)\n' + \
                '-n NAME, --name NAME  Name for the analysis (default: name based on input file name)\n'
            ))
            sys.exit()

        parser = CRISPRessoShared.getCRISPRessoArgParser("Pooled", parser_title = 'CRISPRessoPooled Parameters')
        # parser.add_argument('-f', '--amplicons_file', type=str,  help='Amplicons description file. This file is a tab-delimited text file with up to 14 columns (2 required):\
        # \namplicon_name:  an identifier for the amplicon (must be unique).\
        # \namplicon_seq:  amplicon sequence used in the experiment.\
        # \nguide_seq (OPTIONAL):  sgRNA sequence used for this amplicon without the PAM sequence. Multiple guides can be given separated by commas and not spaces.\
        # \nexpected_hdr_amplicon_seq (OPTIONAL): expected amplicon sequence in case of HDR.\
        # \ncoding_seq (OPTIONAL): Subsequence(s) of the amplicon corresponding to coding sequences. If more than one separate them by commas and not spaces.\
        # \nprime_editing_pegRNA_spacer_seq (OPTIONAL): pegRNA spacer sgRNA sequence used in prime editing. The spacer should not include the PAM sequence. The sequence should be given in the RNA 5\'->3\' order, so for Cas9, the PAM would be on the right side of the given sequence.\
        # \nprime_editing_nicking_guide_seq (OPTIONAL): Nicking sgRNA sequence used in prime editing. The sgRNA should not include the PAM sequence. The sequence should be given in the RNA 5\'->3\' order, so for Cas9, the PAM would be on the right side of the sequence.\
        # \nprime_editing_pegRNA_extension_seq (OPTIONAL): Extension sequence used in prime editing. The sequence should be given in the RNA 5\'->3\' order, such that the sequence starts with the RT template including the edit, followed by the Primer-binding site (PBS).\
        # \nprime_editing_pegRNA_scaffold_seq (OPTIONAL): If given, reads containing any of this scaffold sequence before extension sequence (provided by --prime_editing_extension_seq) will be classified as \'Scaffold-incorporated\'. The sequence should be given in the 5\'->3\' order such that the RT template directly follows this sequence. A common value ends with \'GGCACCGAGUCGGUGC\'.\
        # \nprime_editing_pegRNA_scaffold_min_match_length (OPTIONAL): Minimum number of bases matching scaffold sequence for the read to be counted as \'Scaffold-incorporated\'. If the scaffold sequence matches the reference sequence at the incorporation site, the minimum number of bases to match will be minimally increased (beyond this parameter) to disambiguate between prime-edited and scaffold-incorporated sequences.\
        # \nprime_editing_override_prime_edited_ref_seq (OPTIONAL): If given, this sequence will be used as the prime-edited reference sequence. This may be useful if the prime-edited reference sequence has large indels or the algorithm cannot otherwise infer the correct reference sequence.\
        # \nquantification_window_coordinates (OPTIONAL): Bp positions in the amplicon sequence specifying the quantification window. This parameter overrides values of the "--quantification_window_center", "-- cleavage_offset", "--window_around_sgrna" or "-- window_around_sgrna" values. Any indels/substitutions outside this window are excluded. Indexes are 0-based, meaning that the first nucleotide is position 0. Ranges are separated by the dash sign like "start-stop", and multiple ranges can be separated by the underscore (\_). A value of 0 disables this filter. (can be comma-separated list of values, corresponding to amplicon sequences given in --amplicon_seq e.g. 5-10,5-10_20-30 would specify the 5th-10th bp in the first reference and the 5th-10th and 20th-30th bp in the second reference) (default: None)\
        # \nquantification_window_size (OPTIONAL): Defines the size (in bp) of the quantification window extending from the position specified by the "--cleavage_offset" or "--quantification_window_center" parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp.\
        # \nquantification_window_center (OPTIONAL): Center of quantification window to use within respect to the 3\' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17.', default='')

        # #tool specific optional
        # parser.add_argument('--gene_annotations', type=str, help='Gene Annotation Table from UCSC Genome Browser Tables (http://genome.ucsc.edu/cgi-bin/hgTables?command=start), \
        # please select as table "knownGene", as output format "all fields from selected table" and as file returned "gzip compressed"', default='')
        # # rationale for setting the default scores:
        # # --end-to-end - no clipping, match bonus -ma is set to 0
        # # -N 0 number of mismatches allowed in seed alignment
        # # --np 0 where read (or ref have ambiguous character (N)) penalty is 0
        # # -mp 3,2 mismatch penalty - set max mismatch to -3 to coincide with the gap extension penalty (2 is the default min mismatch penalty)
        # # --score-min L,-5,-3*(1-H) For a given homology score, we allow up to (1-H) mismatches (-3) or gap extensions (-3) and one gap open (-5). This score translates to -5 + -3(1-H)L where L is the sequence length
        # parser.add_argument('--bowtie2_options_string', type=str, help='Override options for the Bowtie2 alignment command. By default, this is " --end-to-end -N 0 --np 0 -mp 3,2 --score-min L,-5,-3(1-H)" where H is the default homology score.', default='')
        # parser.add_argument('--use_legacy_bowtie2_options_string', help='Use legacy (more stringent) Bowtie2 alignment parameters: " -k 1 --end-to-end -N 0 --np 0 ".', action='store_true')
        # parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=1000)
        # parser.add_argument('--skip_failed',  help='Continue with pooled analysis even if one sample fails', action='store_true')
        # parser.add_argument('--skip_reporting_problematic_regions', help='Skip reporting of problematic regions. By default, when both amplicons (-f) and genome (-x) are provided, problematic reads that align to the genome but to positions other than where the amplicons align are reported as problematic', action='store_true')
        # parser.add_argument('--crispresso_command', help='CRISPResso command to call', default='CRISPResso')
        # parser.add_argument('--compile_postrun_references', help='If set, a file will be produced which compiles the reference sequences of frequent amplicons.', action='store_true')
        # parser.add_argument('--compile_postrun_reference_allele_cutoff', type=float, help='Only alleles with at least this percentage frequency in the population will be reported in the postrun analysis. This parameter is given as a percent, so 30 is 30%%.', default=30)
        # parser.add_argument('--alternate_alleles', type=str, help='Path to tab-separated file with alternate allele sequences for pooled experiments. This file has the columns "region_name","reference_seqs", and "reference_names" and gives the reference sequences of alternate alleles that will be passed to CRISPResso for each individual region for allelic analysis. Multiple reference alleles and reference names for a given region name are separated by commas (no spaces).', default='')
        # parser.add_argument('--limit_open_files_for_demux', help='If set, only one file will be opened during demultiplexing of read alignment locations. This will be slightly slower as the reads must be sorted, but may be necessary if the number of amplicons is greater than the number of files that can be opened due to OS constraints.', action='store_true')
        # parser.add_argument('--aligned_pooled_bam', type=str, help='Path to aligned input for CRISPRessoPooled processing. If this parameter is specified, the alignments in the given bam will be used to demultiplex reads. If this parameter is not set (default), input reads provided by --fastq_r1 (and optionally --fastq_r2) will be aligned to the reference genome using bowtie2. If the input bam is given, the corresponding reference fasta must also be given to extract reference genomic sequences via the parameter --bowtie2_index. Note that if the aligned reads are paired-end sequenced, they should already be merged into 1 read (e.g. via Flash) before alignment.', default=None)
        # parser.add_argument('--demultiplex_only_at_amplicons', help='If set, and an amplicon file (--amplicons_file) and reference sequence (--bowtie2_index) are provided, reads overlapping alignment positions of amplicons will be demultiplexed and assigned to that amplicon. If this flag is not set, the entire genome will be demultiplexed and reads with the same start and stop coordinates as an amplicon will be assigned to that amplicon.', action='store_true')

        args = parser.parse_args()

        CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = {'fastq_r1', 'fastq_r2', 'amplicon_seq', 'amplicon_name', 'output_folder', 'name', 'zip_output', 'split_interleaved_input'}
        crispresso_options_for_pooled = list(crispresso_options-options_to_ignore)

        files_to_remove = []

        OUTPUT_DIRECTORY = 'CRISPRessoPooled_on_{0}'.format(normalize_name(args.name, args.fastq_r1, args.fastq_r2, args.aligned_pooled_bam))

        if args.output_folder:
            OUTPUT_DIRECTORY = os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory

        try:
            info('Creating Folder %s' % OUTPUT_DIRECTORY)
            os.makedirs(OUTPUT_DIRECTORY)
            info('Done!')
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename = _jp('CRISPRessoPooled_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        logger.addHandler(CRISPRessoShared.StatusHandler(_jp('CRISPRessoPooled_status.txt')))

        if args.zip_output and not args.place_report_in_output_folder:
            logger.warn('Invalid arguement combination: If zip_output is True then place_report_in_output_folder must also be True. Setting place_report_in_output_folder to True.')
            args.place_report_in_output_folder = True

        info('Checking dependencies...')

        if check_samtools() and check_bowtie2():
            info('All the required dependencies are present!', {'percent_complete': 0})
        else:
            sys.exit(1)

        #check files
        if args.aligned_pooled_bam is not None:
            CRISPRessoShared.check_file(args.aligned_pooled_bam)
            if args.fastq_r1 != '':
                raise CRISPRessoShared.BadParameterException('Arguments for input fastq cannot be provided when input bam is also provided. Please provide either input reads (--fastq_r1) or input bam alignemnts (--aligned_pooled_bam), but not both.')
            if args.bowtie2_index == '':
                raise CRISPRessoShared.BadParameterException('The bowtie2 index must be provided when the aligned pooled bam is given in order to extract reference sequences for alignment. Please provide the bowtie2 index file with the parameter -x or --bowtie2_index.')
            if args.trim_sequences:
                raise CRISPRessoShared.BadParameterException('Cannot trim input sequences if input bam is provided.')


        elif args.fastq_r1:
            if args.fastq_r1 == '':
                raise CRISPRessoShared.BadParameterException('')
            CRISPRessoShared.check_file(args.fastq_r1)
            CRISPRessoShared.assert_fastq_format(args.fastq_r1)
            if args.fastq_r2:
                CRISPRessoShared.check_file(args.fastq_r2)
                CRISPRessoShared.assert_fastq_format(args.fastq_r2)
        else:
            parser.print_help()
            raise CRISPRessoShared.BadParameterException('Please provide input data for pooled analysis e.g. using the --fastq_r1 parameter.')

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
            info('Amplicon description file and bowtie2 reference genome index files provided. The analysis will be perfomed using the reads that are aligned only to the amplicons provided and not to other genomic regions.')
        else:
            error('Please provide the amplicons description file (-f or --amplicons_file option) or the bowtie2 reference genome index file (-x or --bowtie2_index option) or both.')
            sys.exit(1)

        bowtie2_options_string = args.bowtie2_options_string
        if args.bowtie2_options_string == "":
            if args.use_legacy_bowtie2_options_string:
                bowtie2_options_string = '-k 1 --end-to-end -N 0 --np 0'
            else:
                homology_param = -3 * (1-(args.default_min_aln_score/100.0))
                bowtie2_options_string = " --end-to-end -N 0 --np 0 --mp 3,2 --score-min L,-5," + str(homology_param) + " "



        if args.alternate_alleles:
            CRISPRessoShared.check_file(args.alternate_alleles)

        # for computation performed in CRISPRessoPooled (e.g. bowtie alignment, etc) use n_processes_for_pooled
        n_processes_for_pooled = 1
        if args.n_processes == "max":
            n_processes_for_pooled = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes_for_pooled = int(args.n_processes)

        # here, we set args.n_processes as 1 because this value is propagated to sub-CRISPResso runs (not for usage in CRISPRessoWGS)
        args.n_processes = 1

        ####TRIMMING AND MERGING

        crispresso2_info_file = os.path.join(
            OUTPUT_DIRECTORY, 'CRISPResso2Pooled_info.json',
        )
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)
        crispresso2_info['running_info']['finished_steps'] = {}

        #keep track of args to see if it is possible to skip computation steps on rerun
        can_finish_incomplete_run = False
        if args.no_rerun:
            if os.path.exists(crispresso2_info_file):
                previous_run_data = CRISPRessoShared.load_crispresso_info(crispresso_info_file_path=crispresso2_info_file)
                if previous_run_data['running_info']['version'] == CRISPRessoShared.__version__:
                    args_are_same = True
                    for arg in vars(args):
                        if arg == "no_rerun" or arg == "debug" or arg == "n_processes":
                            continue
                        if arg not in vars(previous_run_data['running_info']['args']):
                            info('Comparing current run to previous run: old run had argument ' + str(arg) + ' \nRerunning.')
                            args_are_same = False
                        elif str(getattr(previous_run_data['running_info']['args'], arg)) != str(getattr(args, arg)):
                            info('Comparing current run to previous run:\n\told argument ' + str(arg) + ' = ' + str(getattr(previous_run_data['running_info']['args'], arg)) + '\n\tnew argument: ' + str(arg) + ' = ' + str(getattr(args, arg)) + '\nRerunning.')
                            args_are_same = False

                    if args_are_same:
                        if 'end_time_string' in previous_run_data:
                            info('Analysis already completed on %s!'%previous_run_data['running_info']['end_time_string'])
                            sys.exit(0)
                        else:
                            can_finish_incomplete_run = True
                            #add previous run info to this run
                            if 'finished_steps' in previous_run_data['running_info']:
                                for key in previous_run_data['running_info']['finished_steps'].keys():
                                    crispresso2_info['running_info']['finished_steps'][key] = previous_run_data['running_info']['finished_steps'][key]
                                    if args.debug:
                                        info('finished: ' + key)
                else:
                    info('The no_rerun flag is set, but this analysis will be rerun because the existing run was performed using an old version of CRISPResso (' + str(previous_run_data['running_info']['version']) + ').')


        crispresso_cmd_to_write = ' '.join(sys.argv)
        if args.write_cleaned_report:
            cmd_copy = sys.argv[:]
            cmd_copy[0] = 'CRISPRessoPooled'
            for i in range(len(cmd_copy)):
                if os.sep in cmd_copy[i]:
                    cmd_copy[i] = os.path.basename(cmd_copy[i])

            crispresso_cmd_to_write = ' '.join(cmd_copy) #clean command doesn't show the absolute path to the executable or other files
        crispresso2_info['running_info']['command_used'] = crispresso_cmd_to_write

        #write this file early on so we can check the params if we have to rerun
        CRISPRessoShared.write_crispresso_info(
            crispresso2_info_file, crispresso2_info,
        )

        with open(log_filename, 'w+') as outfile:
            outfile.write('CRISPResso version %s\n[Command used]:\n%s\n\n[Execution log]:\n' %(CRISPRessoShared.__version__, crispresso_cmd_to_write))

        info('Processing input', {'percent_complete': 5})

        if args.split_interleaved_input:
            info('Splitting paired end single fastq file into two files...')
            args.fastq_r1, args.fastq_r2 = CRISPRessoShared.split_interleaved_fastq(
                args.fastq_r1,
                output_filename_r1=_jp('{0}_splitted_r1.fastq.gz'.format(os.path.basename(args.fastq_r1).replace('.fastq', '').replace('.gz', ''))),
                output_filename_r2=_jp('{0}_splitted_r2.fastq.gz'.format(os.path.basename(args.fastq_r1).replace('.fastq', '').replace('.gz', ''))),
            )
            files_to_remove += [args.fastq_r1, args.fastq_r2]

        # perform read trimming if necessary
        if args.aligned_pooled_bam is not None:
            # don't trim reads in aligned bams
            pass
        # read filtering (for quality) is done at the individual crispresso run
        elif args.fastq_r2 == '':  # single end reads
            # check if we need to trim
            if not args.trim_sequences:
                # create a symbolic link
                symlink_filename = _jp(os.path.basename(args.fastq_r1))
                CRISPRessoShared.force_symlink(os.path.abspath(args.fastq_r1), symlink_filename)
                output_forward_filename = symlink_filename
            else:
                output_forward_filename = _jp('reads.trimmed.fq.gz')
                # Trimming with trimmomatic
                info('Trimming sequences with Trimmomatic...', {'percent_complete': 7})
                cmd = '%s SE -phred33 %s %s %s >>%s 2>&1'\
                    % (args.trimmomatic_command, args.fastq_r1,
                        output_forward_filename,
                        args.trimmomatic_options_string,
                        log_filename)
                # print cmd
                TRIMMOMATIC_STATUS = sb.call(cmd, shell=True)

                if TRIMMOMATIC_STATUS:
                    raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

            processed_output_filename = output_forward_filename

        else:  # paired end reads case
            if not args.trim_sequences:
                output_forward_paired_filename = args.fastq_r1
                output_reverse_paired_filename = args.fastq_r2
            else:
                info('Trimming sequences with Trimmomatic...', {'percent_complete': 7})
                output_forward_paired_filename = _jp('output_forward_paired.fq.gz')
                output_forward_unpaired_filename = _jp('output_forward_unpaired.fq.gz')
                output_reverse_paired_filename = _jp('output_reverse_paired.fq.gz')
                output_reverse_unpaired_filename = _jp('output_reverse_unpaired.fq.gz')

                # Trimming with trimmomatic
                cmd = '%s PE -phred33 %s  %s %s  %s  %s  %s %s >>%s 2>&1'\
                    % (args.trimmomatic_command,
                        args.fastq_r1, args.fastq_r2, output_forward_paired_filename,
                        output_forward_unpaired_filename, output_reverse_paired_filename,
                        output_reverse_unpaired_filename, args.trimmomatic_options_string, log_filename)
                # print cmd
                TRIMMOMATIC_STATUS = sb.call(cmd, shell=True)
                if TRIMMOMATIC_STATUS:
                    raise TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

                info('Done!')

            max_overlap_string = ""
            min_overlap_string = ""
            if args.max_paired_end_reads_overlap:
                max_overlap_string = "--max-overlap " + str(args.max_paired_end_reads_overlap)
            if args.min_paired_end_reads_overlap:
                min_overlap_string = "--min-overlap " + str(args.min_paired_end_reads_overlap)
            # Merging with Flash
            info('Merging paired sequences with Flash...', {'percent_complete': 10})
            cmd = args.flash_command+' --allow-outies %s %s %s %s -z -d %s >>%s 2>&1' %\
                (output_forward_paired_filename,
                 output_reverse_paired_filename,
                 max_overlap_string,
                 min_overlap_string,
                 OUTPUT_DIRECTORY, log_filename)

            if args.debug:
                info('Flash command: %s'%cmd)

            FLASH_STATUS = sb.call(cmd, shell=True)
            if FLASH_STATUS:
                raise FlashException('Flash failed to run, please check the log file.')

            flash_hist_filename = _jp('out.hist')
            flash_histogram_filename = _jp('out.histogram')
            flash_not_combined_1_filename = _jp('out.notCombined_1.fastq.gz')
            flash_not_combined_2_filename = _jp('out.notCombined_2.fastq.gz')

            processed_output_filename = _jp('out.extendedFrags.fastq.gz')

            if args.force_merge_pairs:
                old_flashed_filename = processed_output_filename
                new_merged_filename = _jp('out.forcemerged_uncombined.fastq.gz')
                num_reads_force_merged = CRISPRessoShared.force_merge_pairs(flash_not_combined_1_filename, flash_not_combined_2_filename, new_merged_filename)
                new_output_filename = _jp('out.forcemerged.fastq.gz')
                merge_command = "cat %s %s > %s"%(processed_output_filename, new_merged_filename, new_output_filename)
                MERGE_STATUS = sb.call(merge_command, shell=True)
                if MERGE_STATUS:
                    raise FlashException('Force-merging read pairs failed to run, please check the log file.')
                processed_output_filename = new_output_filename

            info('Done!')

        if can_finish_incomplete_run and 'count_input_reads' in crispresso2_info['running_info']['finished_steps']:
            (N_READS_INPUT, N_READS_AFTER_PREPROCESSING) = crispresso2_info['running_info']['finished_steps']['count_input_reads']
        # count reads
        else:
            if args.aligned_pooled_bam is not None:
                N_READS_INPUT = get_n_reads_bam(args.aligned_pooled_bam)
                N_READS_AFTER_PREPROCESSING = N_READS_INPUT
            else:
                N_READS_INPUT = get_n_reads_fastq(args.fastq_r1)
                if args.split_interleaved_input:
                    N_READS_INPUT /= 2
                N_READS_AFTER_PREPROCESSING = get_n_reads_fastq(processed_output_filename)

            crispresso2_info['running_info']['finished_steps']['count_input_reads'] = (N_READS_INPUT, N_READS_AFTER_PREPROCESSING)
            CRISPRessoShared.write_crispresso_info(
                crispresso2_info_file, crispresso2_info
            )

        # load gene annotation
        if args.gene_annotations:
            info('Loading gene coordinates from annotation file: %s...' % args.gene_annotations)
            try:
                df_genes = pd.read_csv(args.gene_annotations, compression='gzip', sep="\t")
                df_genes.txEnd = df_genes.txEnd.astype(int)
                df_genes.txStart = df_genes.txStart.astype(int)
                df_genes.head()
            except Exception:
                info('Failed to load the gene annotations file.')

        # possible column names accepted in amplicon input file
        default_input_amplicon_headers = ['amplicon_name', 'amplicon_seq', 'guide_seq', 'expected_hdr_amplicon_seq', 'coding_seq',
                     'prime_editing_pegRNA_spacer_seq', 'prime_editing_nicking_guide_seq',
                     'prime_editing_pegRNA_extension_seq', 'prime_editing_pegRNA_scaffold_seq',
                     'prime_editing_pegRNA_scaffold_min_match_length', 'prime_editing_override_prime_edited_ref_seq',
                     'quantification_window_coordinates', 'quantification_window_size', 'quantification_window_center']

        if RUNNING_MODE == 'ONLY_AMPLICONS' or RUNNING_MODE == 'AMPLICONS_AND_GENOME':
            # open amplicons file
            with open(args.amplicons_file, 'r') as amplicons_fin:

                head_line = amplicons_fin.readline().strip()
                while head_line[0] == "#":  # read past comments
                    head_line = amplicons_fin.readline()
                header_els = head_line.split('\t')

            head_lookup = CRISPRessoShared.get_crispresso_options_lookup("Core")  # dict of qwc -> quantification_window_coordinates

            # add legacy CRISPRessoPooled headers to the head_lookup
            # lowercase input header names for matching - they'll get fixed in the matching to default_input_amplicon_headers
            head_lookup['sgrna'] = 'guide_seq'
            head_lookup['expected_hdr'] = 'expected_hdr_amplicon_seq'
            head_lookup['name'] = 'amplicon_name'
            head_lookup['sgrna_sequence'] = 'guide_seq'
            head_lookup['expected_amplicon_after_hdr'] = 'expected_hdr_amplicon_seq'
            head_lookup['coding_sequence'] = 'coding_seq'

            lowercase_default_amplicon_headers = {h.lower(): h for h in default_input_amplicon_headers}

            headers = []
            unmatched_headers = []
            has_header = False
            for head in header_els:
                # Header based on header provided
                # Look up long name (e.g. qwc -> quantification_window_coordinates)
                #   as well as custom legacy headers above
                long_head = head
                if head in head_lookup:
                    long_head = head_lookup[head]
                if head.lower() in head_lookup:
                    long_head = head_lookup[head.lower()]
                long_head = long_head.lower()

                match = difflib.get_close_matches(long_head, lowercase_default_amplicon_headers, n=1)
                if not match:
                    unmatched_headers.append(head)
                else:
                    has_header = True
                    headers.append(lowercase_default_amplicon_headers[match[0]])
                    if args.debug:
                        info(f'Matching header {head} with {lowercase_default_amplicon_headers[match[0]]}.')

            if len(headers) > 5 and not has_header:
                raise CRISPRessoShared.BadParameterException('Incorrect number of columns provided without header.')
            elif has_header and len(unmatched_headers) > 0:
                raise CRISPRessoShared.BadParameterException('Unable to match headers: ' + str(unmatched_headers))

            if not has_header:
                # Default header
                headers = []
                for i in range(len(header_els)):
                    headers.append(default_input_amplicon_headers[i])

            if args.debug:
                info(f'Header variable names in order: {headers}')

            # load and validate template file
            df_template = pd.read_csv(args.amplicons_file, names=headers, comment='#', sep='\t', dtype={'amplicon_name':str})

            if has_header:
                df_template.drop(0, axis=0, inplace=True)
                info('Detected header in amplicon file.')


            #remove empty amplicons/lines
            df_template.dropna(subset=['amplicon_seq'], inplace=True)
            df_template.dropna(subset=['amplicon_name'], inplace=True)

            df_template.amplicon_seq=df_template.amplicon_seq.apply(CRISPRessoShared.capitalize_sequence)
            if 'expected_hdr_amplicon_seq' in df_template.columns:
                df_template.expected_hdr_amplicon_seq=df_template.expected_hdr_amplicon_seq.apply(CRISPRessoShared.capitalize_sequence)
            if 'guide_seq' in df_template.columns:
                df_template.guide_seq=df_template.guide_seq.apply(CRISPRessoShared.capitalize_sequence)
            if 'coding_seq' in df_template.columns:
                df_template.coding_seq=df_template.coding_seq.apply(CRISPRessoShared.capitalize_sequence)

            if not len(df_template.amplicon_seq.unique())==df_template.shape[0]:
                duplicated_entries = df_template.amplicon_seq[df_template.amplicon_seq.duplicated()]
                raise Exception('The amplicon sequences must be distinct! (Duplicated entries: ' + str(duplicated_entries.values) + ')')

            if not len(df_template.amplicon_name.unique())==df_template.shape[0]:
                duplicated_entries = df_template.amplicon_name[df_template.Name.duplicated()]
                raise Exception('The amplicon names must be distinct! (Duplicated names: ' + str(duplicated_entries.values) + ')')

            df_template=df_template.set_index('amplicon_name')
            df_template.index=df_template.index.to_series().str.replace(' ', '_')

            for idx, row in df_template.iterrows():
                wrong_nt=CRISPRessoShared.find_wrong_nt(row.amplicon_seq)
                if wrong_nt:
                    raise NTException('The amplicon sequence %s contains wrong characters:%s' % (idx, ' '.join(wrong_nt)))

                if 'guide_seq' in df_template.columns and not pd.isnull(row.guide_seq):
                    cut_points = []
                    guides = row.guide_seq.strip().upper().split(',')
                    guide_qw_centers = CRISPRessoShared.set_guide_array(args.quantification_window_center, guides, 'guide quantification center')
                    for idx, current_guide_seq in enumerate(guides):

                        wrong_nt = CRISPRessoShared.find_wrong_nt(current_guide_seq)
                        if wrong_nt:
                            raise NTException('The sgRNA sequence %s contains wrong characters:%s'  % (current_guide_seq, ' '.join(wrong_nt)))

                        offset_fw=guide_qw_centers[idx]+len(current_guide_seq)-1
                        offset_rc=(-guide_qw_centers[idx])-1
                        cut_points+=[m.start() + offset_fw for \
                                    m in re.finditer(current_guide_seq,  row.amplicon_seq)]+[m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse_complement(current_guide_seq),  row.amplicon_seq)]

                    if not cut_points:
                        warn('\nThe guide sequence/s provided: %s is(are) not present in the amplicon sequence:%s! \nNOTE: The guide will be ignored for the analysis. Please check your input!' % (row.guide_seq, row.amplicon_seq))
                        df_template.iloc[idx,df_template.columns.get_loc('guide_seq')] = ''



        if RUNNING_MODE=='ONLY_AMPLICONS':
            #create a fasta file with all the amplicons
            amplicon_fa_filename=_jp('AMPLICONS.fa')
            fastq_gz_amplicon_filenames=[]
            with open(amplicon_fa_filename, 'w+') as outfile:
                for idx, row in df_template.iterrows():
                    if row['amplicon_seq']:
                        outfile.write('>%s\n%s\n' %(CRISPRessoShared.clean_filename('AMPL_'+idx), row['amplicon_seq']))

                        #create place-holder fastq files
                        fastq_gz_amplicon_filenames.append(_jp('%s.fastq.gz' % CRISPRessoShared.clean_filename('AMPL_'+idx)))
                        open(fastq_gz_amplicon_filenames[-1], 'w+').close()

            df_template['Demultiplexed_fastq.gz_filename']=fastq_gz_amplicon_filenames
            info('Creating a custom index file with all the amplicons...')
            custom_index_filename=_jp('CUSTOM_BOWTIE2_INDEX')
            sb.call('bowtie2-build %s %s >>%s 2>&1' %(amplicon_fa_filename, custom_index_filename, log_filename), shell=True)


            #align the file to the amplicons (MODE 1)
            info('Align reads to the amplicons...')
            bam_filename_amplicons= _jp('CRISPResso_AMPLICONS_ALIGNED.bam')
            aligner_command= 'bowtie2 -x %s -p %s %s -U %s 2>>%s | samtools view -bS - > %s' %(custom_index_filename, n_processes_for_pooled, bowtie2_options_string, processed_output_filename, log_filename, bam_filename_amplicons)


            info('Alignment command: ' + aligner_command, {'percent_complete': 15})
            sb.call(aligner_command, shell=True)

            N_READS_ALIGNED = get_n_aligned_bam(bam_filename_amplicons)

            if args.limit_open_files_for_demux:
                bam_iter = CRISPRessoShared.get_command_output(
                    '(samtools sort {bam_file} | samtools view -F 4) 2>> {log_file}'.format(
                        bam_file=bam_filename_amplicons,
                        log_file=log_filename,
                    ),
                )
                curr_file, curr_chr = None, None
                for bam_line in bam_iter:
                    bam_line_els = bam_line.split('\t')
                    if len(bam_line_els) < 9:
                        if args.debug:
                            info('ERROR got unexpected line from bam: {0} with els: {1}'.format(
                                bam_line, str(bam_line_els),
                            ))
                        continue
                    line_chr = bam_line_els[2]

                    # at the first line open new file, or at next amplicon
                    # close previous file and open new one
                    if curr_chr != line_chr:
                        if curr_file is not None:
                            curr_file.close()
                        curr_file = gzip.open(
                            _jp('{0}.fastq.gz'.format(line_chr)),
                            'wt',
                        )
                    curr_file.write('@{read_name}\n{seq}\n+\n{qual}\n'.format(
                        read_name=bam_line_els[0],
                        seq=bam_line_els[9],
                        qual=bam_line_els[10],
                    ))
                    curr_chr = line_chr
                if curr_file is not None:
                    curr_file.close()
            else:
                s1 = r"samtools view -F 4 %s 2>>%s | grep -v ^'@'" % (bam_filename_amplicons,log_filename)
                s2 = r'''|awk '{ gzip_filename=sprintf("gzip >> OUTPUTPATH%s.fastq.gz",$3);\
                print "@"$1"\n"$10"\n+\n"$11  | gzip_filename;}' '''

                cmd = s1+s2.replace('OUTPUTPATH', _jp(''))
                sb.call(cmd, shell=True)

            alternate_alleles = {}
            if args.alternate_alleles:
                with open(args.alternate_alleles, 'r') as alt_in:
                    alt_head_els = alt_in.readline().lower().rstrip().split("\t")
                    region_name_ind = 0
                    allele_seq_ind = 1
                    allele_name_ind = 2
                    if 'region_name' in alt_head_els:
                        region_name_ind = alt_head_els.index('region_name')
                    else:
                        alt_in.seek(0, 0) #start at beginning of file -- no header
                    if 'reference_seqs' in alt_head_els:
                        allele_seq_ind = alt_head_els.index('reference_seqs')
                    if 'reference_names' in alt_head_els:
                        allele_name_ind = alt_head_els.index('reference_names')
                    for line in alt_in:
                        line_els = line.rstrip().split("\t")
                        alternate_alleles[line_els[region_name_ind]] = (line_els[allele_seq_ind], line_els[allele_name_ind])

            info('Demultiplex reads and run CRISPResso on each amplicon...')
            n_reads_aligned_amplicons=[]
            crispresso_cmds = []
            for idx, row in df_template.iterrows():
                this_n_reads = get_n_reads_fastq(row['Demultiplexed_fastq.gz_filename'])
                n_reads_aligned_amplicons.append(this_n_reads)
                info('\n Processing:%s with %d reads'%(idx,this_n_reads))
                this_amp_seq = row['amplicon_seq']
                this_amp_name_string = ""
                if idx in alternate_alleles:
                    new_refs, new_names = alternate_alleles[idx]
                    this_amp_seq = new_refs
                    this_amp_name_string = "-an " + new_names

                crispresso_cmd= args.crispresso_command + ' -r1 %s -a %s %s -o %s --name %s' % (row['Demultiplexed_fastq.gz_filename'], this_amp_seq, this_amp_name_string, OUTPUT_DIRECTORY, idx)

                if n_reads_aligned_amplicons[-1] > args.min_reads_to_use_region:
                    this_run_args_from_amplicons_file = {}
                    for column_name in default_input_amplicon_headers:
                        if column_name in df_template.columns and row[column_name] and not pd.isnull(row[column_name]):
                            this_run_args_from_amplicons_file[column_name] = row[column_name]

                    #first, set the general CRISPResso options for this sub-run (e.g. plotting options, etc)
                    #note that the crispresso_options_for_pooled doesn't include e.g. amplicon_seq so when someone calls CRISPRessoPooled with -a that won't get passed on here
                    crispresso_cmd = CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, crispresso_options_for_pooled, args)
                    #next set the per-amplicon options we read from the Amplicons file (and are stored in this_run_args_from_amplicons_file)
                    crispresso_cmd = CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, default_input_amplicon_headers, this_run_args_from_amplicons_file)

                    crispresso_cmds.append(crispresso_cmd)

                else:
                    warn('Skipping amplicon [%s] because no reads align to it\n'% idx)
            print(123454321)
            print(crispresso_cmds)
            CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, n_processes_for_pooled, 'amplicon', args.skip_failed, start_end_percent=(16, 80))
            # Initialize array to track failed runs
            failed_batch_arr = []
            failed_batch_arr_desc = []
            for cmd in crispresso_cmds:

                # Extract the folder name from the CRISPResso command
                folder_name_regex = re.search(r'-o\s+\S+\s+--name\s+(\S+)', cmd)
                if folder_name_regex:
                    folder_name = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso_on_%s' % folder_name_regex.group(1))
                    failed_run_bool, failed_status_string = CRISPRessoShared.check_if_failed_run(folder_name, info)
                    if failed_run_bool:
                        failed_batch_arr.append(folder_name_regex.group(1))
                        failed_batch_arr_desc.append(failed_status_string)

            # Store the failed runs in crispresso2_info for later use
            crispresso2_info['results']['failed_batch_arr'] = failed_batch_arr
            crispresso2_info['results']['failed_batch_arr_desc'] = failed_batch_arr_desc

            df_template['n_reads']=n_reads_aligned_amplicons
            df_template['n_reads_aligned_%']=df_template['n_reads']/float(N_READS_ALIGNED)*100
            df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_AMPLICONS.txt'), sep='\t')



        if RUNNING_MODE=='AMPLICONS_AND_GENOME':
            info('Mapping amplicons to the reference genome...')

            filename_amplicon_aligned_locations = _jp('CRISPResso_amplicon_aligned_locations.txt')
            filename_aligned_amplicons_sam = _jp('CRISPResso_amplicons_aligned.sam')
            filename_aligned_amplicons_sam_log = _jp('CRISPResso_amplicons_aligned.sam.log')
            filename_amplicon_seqs_fasta = _jp('CRISPResso_amplicons_to_align.fa')

            if can_finish_incomplete_run and 'mapping_amplicons_to_reference_genome' in crispresso2_info['running_info']['finished_steps']:
                info('Reading previously-computed alignment of amplicons to genome')
                additional_columns_df = pd.read_csv(filename_amplicon_aligned_locations, sep="\t")
                additional_columns_df.set_index('amplicon_name', inplace=True)
            else:
                #write amplicons as fastq for alignment
                with open(filename_amplicon_seqs_fasta, 'w') as fastas:
                    for idx, row in df_template.iterrows():
                        fastas.write('>%s\n%s\n'%(idx, row.amplicon_seq))

                aligner_command= 'bowtie2 -x %s -p %s %s -f -U %s --no-hd --no-sq 2> %s > %s ' %(args.bowtie2_index, n_processes_for_pooled, bowtie2_options_string, \
                    filename_amplicon_seqs_fasta, filename_aligned_amplicons_sam_log, filename_aligned_amplicons_sam)
                bowtie_status=sb.call(aligner_command, shell=True)
                if bowtie_status:
                        raise Bowtie2Exception('Bowtie2 failed to align amplicons to the genome, please check the output file.')

                additional_columns = []
                with open (filename_aligned_amplicons_sam) as aln:
                    for line in aln.readlines():
                        line_els = line.split("\t")
                        if line_els[2] == "*":
                            info('The amplicon [%s] is not mappable to the reference genome provided!' % line_els[0])
                            additional_columns.append([line_els[0], 'NOT_ALIGNED', 0, -1, '+', ''])
                        else:
                            aln_len = CRISPRessoShared.get_ref_length_from_cigar(line_els[5])
                            seq_start = int(line_els[3])
                            seq_stop = seq_start + aln_len
                            strand = "-" if (int(line_els[1]) & 0x10) else "+"
                            additional_columns.append([line_els[0], line_els[2], seq_start, seq_stop, strand, line_els[9]])
                            info('The amplicon [%s] was mapped to: %s:%d-%d ' % (line_els[0], line_els[2], seq_start, seq_stop))
                additional_columns_df = pd.DataFrame(additional_columns, columns=['amplicon_name', 'chr_id', 'bpstart', 'bpend', 'strand', 'reference_seq']).set_index('amplicon_name')
                additional_columns_df.to_csv(filename_amplicon_aligned_locations, sep="\t", index_label='amplicon_name')

                crispresso2_info['running_info']['finished_steps']['mapping_amplicons_to_reference_genome'] = True
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info
                )

            files_to_remove.append(filename_amplicon_seqs_fasta)
            files_to_remove.append(filename_aligned_amplicons_sam)

            df_template=df_template.join(additional_columns_df)

            df_template.bpstart=df_template.bpstart.astype(int)
            df_template.bpend=df_template.bpend.astype(int)

            #Check reference is the same otherwise throw a warning
            for idx, row in df_template.iterrows():
                if row.amplicon_seq != row.reference_seq and row.amplicon_seq != CRISPRessoShared.reverse_complement(row.reference_seq):
                    warn('The amplicon sequence %s provided:\n%s\n\nis different from the reference sequence(both strands):\n\n%s\n\n%s\n' %(row.name, row.amplicon_seq, row.amplicon_seq, CRISPRessoShared.reverse_complement(row.amplicon_seq)))


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

                cmd_to_uncompress='bowtie2-inspect %s > %s 2>>%s' % (args.bowtie2_index, uncompressed_reference, log_filename)
                sb.call(cmd_to_uncompress, shell=True)

                info('Indexing fasta file with samtools...')
                #!samtools faidx {uncompressed_reference}
                sb.call('samtools faidx %s 2>>%s ' % (uncompressed_reference, log_filename), shell=True)

        # align reads to the genome in an unbiased way
        if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
            bam_filename_genome = _jp('{0}_GENOME_ALIGNED.bam'.format(normalize_name(
                args.name, args.fastq_r1, args.fastq_r2, args.aligned_pooled_bam,
            )))
            # if input bam is provided, don't align reads to the genome and use that bam
            if args.aligned_pooled_bam is not None:
                bam_filename_genome = args.aligned_pooled_bam

            if can_finish_incomplete_run and 'n_reads_aligned_genome' in crispresso2_info['running_info']['finished_steps']:
                info('Using previously-computed alignment of reads to genome')
                N_READS_ALIGNED = crispresso2_info['running_info']['finished_steps']['n_reads_aligned_genome']
            # if aligned bam is provided, count reads aligned to genome
            elif args.aligned_pooled_bam is not None:
                def rreplace(s, old, new):
                    li = s.rsplit(old)
                    return new.join(li)
                if os.path.exists(rreplace(args.aligned_pooled_bam, ".bam", ".bai")) or os.path.exists(args.aligned_pooled_bam+'.bai'):
                    info('Index file for input .bam file exists, skipping generation.')
                else:
                    info('Index file for input .bam file does not exist. Generating bam index file.')
                    sb.call('samtools index %s' % bam_filename_genome, shell=True)

                N_READS_ALIGNED = get_n_aligned_bam(bam_filename_genome)
                # save progress up to this point
                crispresso2_info['running_info']['finished_steps']['n_reads_aligned_genome'] = N_READS_ALIGNED
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info,
                )
            # otherwise, align reads to the genome and count reads
            else:
                info('Aligning reads to the provided genome index...')
                aligner_command = 'bowtie2 -x %s -p %s %s -U %s 2>>%s| samtools view -bS - | samtools sort -@ %d - -o %s' %(args.bowtie2_index, n_processes_for_pooled,
                    bowtie2_options_string, processed_output_filename, log_filename, n_processes_for_pooled, bam_filename_genome)
                if args.debug:
                    info('Aligning with command: ' + aligner_command)
                sb.call(aligner_command, shell=True)

                sb.call('samtools index %s' % bam_filename_genome, shell=True)

                N_READS_ALIGNED = get_n_aligned_bam(bam_filename_genome)

                # save progress up to this point
                crispresso2_info['running_info']['finished_steps']['n_reads_aligned_genome'] = N_READS_ALIGNED
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info,
                )

            MAPPED_REGIONS = _jp('MAPPED_REGIONS/')
            REPORT_ALL_DEPTH = _jp('REPORT_READS_ALIGNED_TO_GENOME_ALL_DEPTHS.txt')

            if can_finish_incomplete_run and 'genome_demultiplexing' in crispresso2_info['running_info']['finished_steps'] and os.path.isfile(REPORT_ALL_DEPTH):
                info('Using previously-computed demultiplexing of genomic reads')
                df_all_demux = pd.read_csv(REPORT_ALL_DEPTH, sep='\t')
                df_all_demux['loc'] = df_all_demux['chr_id']+' ' + df_all_demux['start'].apply(str) + ' '+df_all_demux['end'].apply(str)
                df_all_demux.set_index(['loc'], inplace=True)
            else:
                #REDISCOVER LOCATIONS and DEMULTIPLEX READS

                # first get rid of all files in the output directory
                if os.path.exists(MAPPED_REGIONS):
                    info('Deleting partially-completed demultiplexing in %s...'%MAPPED_REGIONS)
                    cmd = "rm -rf %s" % MAPPED_REGIONS
                    p = sb.call(cmd, shell=True)

                # make the output directory
                os.mkdir(MAPPED_REGIONS)

                # if we should only demultiplex where amplicons aligned... (as opposed to the whole genome)
                if RUNNING_MODE=='AMPLICONS_AND_GENOME' and args.demultiplex_only_at_amplicons:
                    s1 = r'''samtools view -F 0x0004 %s __REGIONCHR__:__REGIONSTART__-__REGIONEND__ 2>>%s |''' % (bam_filename_genome, log_filename)+\
                    r'''awk 'BEGIN{OFS="\t";num_records=0;fastq_filename="__OUTPUTPATH__REGION___REGIONCHR_____REGIONSTART_____REGIONEND__.fastq";} \
                        { \
                            print "@"$1"\n"$10"\n+\n"$11 > fastq_filename; \
                            num_records++; \
                        } \
                    END{ \
                      close(fastq_filename); \
                      system("gzip -f "fastq_filename);  \
                      record_log_str = "__REGIONCHR__\t__REGIONSTART__\t__REGIONEND__\t"num_records"\t"fastq_filename".gz\n"; \
                      print record_log_str > "__DEMUX_CHR_LOGFILENAME__"; \
                    } ' '''
                    cmd = (s1).replace('__OUTPUTPATH__', MAPPED_REGIONS)
                    cmd = cmd.replace("__MIN_READS__", str(args.min_reads_to_use_region))
                    with open(REPORT_ALL_DEPTH, 'w') as f:
                        f.write('chr_id\tstart\tend\tnumber of reads\toutput filename\n')

                    info('Preparing to demultiplex reads aligned to positions overlapping amplicons in the genome...')
                    # make command for each amplicon

                    chr_commands = []
                    chr_output_filenames = []
                    for idx, row in df_template.iterrows():
                        chr_output_filename = _jp('MAPPED_REGIONS/chr%s_%s_%s.info' % (row.chr_id, row.bpstart, row.bpend))
                        sub_chr_command = cmd.replace('__REGIONCHR__', str(row.chr_id)).replace('__REGIONSTART__',str(row.bpstart)).replace('__REGIONEND__',str(row.bpend)).replace("__DEMUX_CHR_LOGFILENAME__", chr_output_filename)
                        chr_commands.append(sub_chr_command)
                        chr_output_filenames.append(chr_output_filename)

                # if we should demultiplex everwhere (not just where amplicons aligned)
                else:
                    # next, create the general demux command
                    # variables like __CHR__ will be subbed out below for each iteration
                    s1 = r'''samtools view -F 0x0004 %s __CHR____REGION__ 2>>%s |''' % (bam_filename_genome, log_filename) + \
                    r'''awk 'BEGIN {OFS="\t"} {bpstart=$4;  bpend=bpstart; split ($6,a,"[MIDNSHP]"); n=0;\
                    for (i=1; i in a; i++){\
                        n+=1+length(a[i]);\
                        if (substr($6,n,1)=="S"){\
                            if (bpend==$4)\
                                bpstart-=a[i];\
                            else \
                                bpend+=a[i]; \
                            }\
                        else if( (substr($6,n,1)!="I")  && (substr($6,n,1)!="H") )\
                                bpend+=a[i];\
                        }\
                        if (($2 % 32)>=16)\
                            print $3,bpstart,bpend,"-",$1,$10,$11;\
                        else\
                            print $3,bpstart,bpend,"+",$1,$10,$11;}' | '''

                    s2 = r'''  sort -k1,1 -k2,2n  | awk \
                     'BEGIN{chr_id="NA";bpstart=-1;bpend=-1; fastq_filename="NA";num_records=0;fastq_records="";fastq_record_sep="";record_log_str = ""}\
                    { if ( (chr_id!=$1) || (bpstart!=$2) || (bpend!=$3) )\
                        {\
                        if (fastq_filename!="NA") {if (num_records < __MIN_READS__){\
                            record_log_str = record_log_str chr_id"\t"bpstart"\t"bpend"\t"num_records"\tNA\n"} \
                    else{print(fastq_records)>fastq_filename;close(fastq_filename); system("gzip -f "fastq_filename); record_log_str = record_log_str chr_id"\t"bpstart"\t"bpend"\t"num_records"\t"fastq_filename".gz\n"} \
                        }\
                        chr_id=$1; bpstart=$2; bpend=$3;\
                        fastq_filename=sprintf("__OUTPUTPATH__REGION_%s_%s_%s.fastq",$1,$2,$3);\
                        num_records = 0;\
                        fastq_records="";\
                        fastq_record_sep="";\
                        }\
                    fastq_records=fastq_records fastq_record_sep "@"$5"\n"$6"\n+\n"$7; \
                    fastq_record_sep="\n"; \
                    num_records++; \
                    } \
                    END{ \
                        if (fastq_filename!="NA") {if (num_records < __MIN_READS__){\
                            record_log_str = record_log_str chr_id"\t"bpstart"\t"bpend"\t"num_records"\tNA\n"} \
                    else{printf("%s",fastq_records)>fastq_filename;close(fastq_filename); system("gzip -f "fastq_filename); record_log_str = record_log_str chr_id"\t"bpstart"\t"bpend"\t"num_records"\t"fastq_filename".gz\n"} \
                        }\
                        print record_log_str > "__DEMUX_CHR_LOGFILENAME__" \
                    }' '''
                    cmd = (s1+s2).replace('__OUTPUTPATH__', MAPPED_REGIONS)
                    cmd = cmd.replace("__MIN_READS__", str(args.min_reads_to_use_region))

                    info('Preparing to demultiplex reads aligned to the genome...')
                    # next, get all of the chromosome names (for parallelization)
                    enumerate_chr_cmd = "samtools view -H %s" % bam_filename_genome
                    p = sb.Popen(enumerate_chr_cmd, shell=True, stdout=sb.PIPE)
                    chr_lines = p.communicate()[0].decode('utf-8').split("\n")
                    chrs = []
                    chr_lens = {}
                    for chr_line in chr_lines:
                        m = re.match(r'@SQ\s+SN:(\S+)\s+LN:(\d+)', chr_line)
                        if m:
                            chrs.append(m.group(1))
                            chr_lens[m.group(1)] = int(m.group(2))

                    chr_commands = []
                    chr_output_filenames = []
                    for chr_str in chrs:
                        chr_cmd = cmd.replace('__CHR__', chr_str)
                        # if we have a lot of reads, split up the chrs too
                        # with a step size of 10M, there are about 220 regions in hg19
                        # with a step size of 5M, there are about 368 regions in hg19
                        chr_step_size = 5000000 #step size for splitting up chrs
                        chr_len = chr_lens[chr_str]
                        if N_READS_ALIGNED > 10000000 and chr_len > chr_step_size*2:
                            curr_pos = 0
                            curr_end = curr_pos + chr_step_size
                            while curr_end < chr_len:
                                # make sure there aren't any reads at this breakpoint
                                n_reads_at_end = get_n_aligned_bam_region(bam_filename_genome, chr_str, curr_end-5, curr_end+5)
                                while n_reads_at_end > 0:
                                    curr_end += 500  # look for another place with no reads
                                    if curr_end >= chr_len:
                                        curr_end = chr_len
                                        break
                                    n_reads_at_end = get_n_aligned_bam_region(bam_filename_genome, chr_str, curr_end-5, curr_end+5)

                                chr_output_filename = _jp('MAPPED_REGIONS/%s_%s_%s.info' % (chr_str, curr_pos, curr_end))
                                sub_chr_command = chr_cmd.replace("__REGION__", ":%d-%d "%(curr_pos, curr_end)).replace("__DEMUX_CHR_LOGFILENAME__",chr_output_filename)
                                chr_commands.append(sub_chr_command)
                                chr_output_filenames.append(chr_output_filename)
                                curr_pos = curr_end
                                curr_end = curr_pos + chr_step_size
                            if curr_end < chr_len:
                                chr_output_filename = _jp('MAPPED_REGIONS/%s_%s_%s.info' % (chr_str, curr_pos, chr_len))
                                sub_chr_command = chr_cmd.replace("__REGION__", ":%d-%d "%(curr_pos, chr_len)).replace("__DEMUX_CHR_LOGFILENAME__",chr_output_filename)
                                chr_commands.append(sub_chr_command)
                                chr_output_filenames.append(chr_output_filename)

                        else:
                            # otherwise do the whole chromosome
                            chr_output_filename = _jp('MAPPED_REGIONS/%s.info' % (chr_str))
                            sub_chr_command = chr_cmd.replace("__REGION__", "").replace("__DEMUX_CHR_LOGFILENAME__",chr_output_filename)
                            chr_commands.append(sub_chr_command)
                            chr_output_filenames.append(chr_output_filename)

                if args.debug:
                    demux_file = _jp('DEMUX_COMMANDS.txt')
                    with open(demux_file, 'w') as fout:
                        fout.write("\n\n\n".join(chr_commands))
                    info('Wrote demultiplexing commands to ' + demux_file)

                info('Demultiplexing reads by location (%d genomic regions)...'%len(chr_commands), {'percent_complete': 85})
                CRISPRessoMultiProcessing.run_parallel_commands(chr_commands, n_processes=n_processes_for_pooled, descriptor='Demultiplexing reads by location', continue_on_fail=args.skip_failed)

                with open(REPORT_ALL_DEPTH, 'w') as f:
                    f.write('chr_id\tstart\tend\tnumber of reads\toutput filename\n')
                    for chr_output_filename in chr_output_filenames:
                        with open(chr_output_filename, 'r') as f_in:
                            for line in f_in:
                                f.write(line)

                df_all_demux = pd.read_csv(REPORT_ALL_DEPTH, sep='\t')
                df_all_demux.sort_values(by=['chr_id', 'start'], inplace=True)
                sum_aligned_reads = df_all_demux["number of reads"].sum()
                # write the sorted file
                df_all_demux.to_csv(REPORT_ALL_DEPTH, sep="\t", index=False, na_rep="NA")
                df_all_demux['loc'] = df_all_demux['chr_id']+' ' + df_all_demux['start'].apply(str) + ' '+df_all_demux['end'].apply(str)
                df_all_demux.set_index(['loc'], inplace=True)

                if sum_aligned_reads == 0:
                    raise NoReadsAlignedException("No reads aligned to the specified genome")

                crispresso2_info['running_info']['finished_steps']['genome_demultiplexing'] = True
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info,
                )

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

        if RUNNING_MODE == 'AMPLICONS_AND_GENOME':
            files_to_match = list(df_all_demux['output filename'].dropna())
            n_reads_aligned_genome = []
            fastq_region_filenames = []

            if can_finish_incomplete_run and 'crispresso_amplicons_and_genome' in crispresso2_info['running_info']['finished_steps']:
                info('Using previously-computed crispresso runs')
                (n_reads_aligned_genome, fastq_region_filenames, files_to_match) = crispresso2_info['running_info']['finished_steps']['crispresso_amplicons_and_genome'];
            else:
                crispresso_cmds = []
                for idx, row in df_template.iterrows():

                    info('Processing amplicon: %s' % idx )

                    # check if we have reads
                    demux_key = row['chr_id'] + ' ' + str(row['bpstart']) + ' ' + str(row['bpend'])
                    if demux_key in df_all_demux.index:
                        demux_row = df_all_demux.loc[demux_key]
                        N_READS = demux_row['number of reads']
                        n_reads_aligned_genome.append(N_READS)
                        fastq_filename_region = str(demux_row['output filename'])
                        if fastq_filename_region == "nan":
                            fastq_filename_region = ""
                        else:
                            if fastq_filename_region in files_to_match:
                                files_to_match.remove(fastq_filename_region)
                        fastq_region_filenames.append(fastq_filename_region)
                        #else:
                             #info('Warning: Fastq filename ' + fastq_filename_region + ' is not in ' + str(files_to_match))
                             #debug here??
                        if N_READS >= args.min_reads_to_use_region and fastq_filename_region != "":
                            info('\nThe amplicon [%s] has enough reads (%d) mapped to it! Running CRISPResso!\n' % (idx, N_READS))

                            crispresso_cmd = args.crispresso_command + ' -r1 %s -o %s --name %s' % (fastq_filename_region, OUTPUT_DIRECTORY, idx)

                            this_run_args_from_amplicons_file = {}
                            for column_name in default_input_amplicon_headers:
                                if column_name in df_template.columns and row[column_name] and not pd.isnull(row[column_name]):
                                    this_run_args_from_amplicons_file[column_name] = row[column_name]

                            #first, set the general CRISPResso options for this sub-run (e.g. plotting options, etc)
                            #note that the crispresso_options_for_pooled doesn't include e.g. amplicon_seq so when someone calls CRISPRessoPooled with -a that won't get passed on here
                            crispresso_cmd = CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, crispresso_options_for_pooled, args)
                            #next set the per-amplicon options we read from the Amplicons file (and are stored in this_run_args_from_amplicons_file)
                            crispresso_cmd = CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, default_input_amplicon_headers, this_run_args_from_amplicons_file)
                            info('Running CRISPResso:%s' % crispresso_cmd)
                            crispresso_cmds.append(crispresso_cmd)

                        else:
                            warn('The amplicon [%s] has too few reads (%d) mapped to it! Skipping the execution of CRISPResso!' % (idx, N_READS))
                    else:
                        fastq_region_filenames.append('')
                        n_reads_aligned_genome.append(0)
                        warn("The amplicon %s doesn't have any reads mapped to it!\n Please check your amplicon sequence." %  idx)

                CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, n_processes_for_pooled, 'amplicon', args.skip_failed, start_end_percent=(15, 85))

                crispresso2_info['running_info']['finished_steps']['crispresso_amplicons_and_genome'] = (n_reads_aligned_genome, fastq_region_filenames, files_to_match)
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info,
                )

            df_template['Amplicon_Specific_fastq.gz_filename'] = fastq_region_filenames
            df_template['n_reads'] = n_reads_aligned_genome
            df_template['n_reads_aligned_%'] = df_template['n_reads']/float(N_READS_ALIGNED)*100

            if args.gene_annotations:
                df_template = df_template.apply(lambda row: find_overlapping_genes(row, df_genes), axis=1)

            df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_GENOME_AND_AMPLICONS.txt'), sep='\t')

            #write another file with the not amplicon regions

            if args.skip_reporting_problematic_regions:
                df_regions = pd.DataFrame(columns=['chr_id', 'bpstart', 'bpend', 'fastq_file', 'n_reads', 'Reference_sequence'])
            else:
                filename_problematic_regions = _jp('REPORTS_READS_ALIGNED_TO_GENOME_NOT_MATCHING_AMPLICONS.txt')
                if can_finish_incomplete_run and 'reporting_problematic_regions' in crispresso2_info['running_info']['finished_steps']:
                    info('Skipping previously-computed reporting of problematic regions')
                    df_regions = pd.read_csv(filename_problematic_regions, sep='\t')
                else:
                    info('Reporting problematic regions...')
                    summarize_region_fastq_input = [f+" "+uncompressed_reference for f in files_to_match] #pass both params to parallel function
                    coordinates = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(summarize_region_fastq_input, summarize_region_fastq_chunk, n_processes=n_processes_for_pooled)
                    df_regions = pd.DataFrame(coordinates, columns=['chr_id', 'bpstart', 'bpend', 'fastq_file', 'n_reads', 'Reference_sequence'])
                    df_regions.dropna(inplace=True) #remove regions in chrUn

                    df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
                    df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])
                    df_regions['n_reads'] = pd.to_numeric(df_regions['n_reads'])

                    df_regions.bpstart = df_regions.bpstart.astype(int)
                    df_regions.bpend = df_regions.bpend.astype(int)

                    df_regions['n_reads_aligned_%'] = df_regions['n_reads']/float(N_READS_ALIGNED)*100

                    if args.gene_annotations:
                        info('Checking overlapping genes...')
                        df_regions = df_regions.apply(lambda row: find_overlapping_genes(row, df_genes), axis=1)

                    df_regions.sort_values(by='n_reads', ascending=False, inplace=True)

                    df_regions.fillna('NA').to_csv(filename_problematic_regions, sep='\t', index=None)

                    crispresso2_info['running_info']['finished_steps']['reporting_problematic_regions'] = True
                    CRISPRessoShared.write_crispresso_info(
                        crispresso2_info_file, crispresso2_info,
                    )

        if RUNNING_MODE=='ONLY_GENOME' :
            # Load regions and build REFERENCE TABLES
            filename_reads_aligned_to_genome_only = _jp('REPORT_READS_ALIGNED_TO_GENOME_ONLY.txt')
            if can_finish_incomplete_run and 'demultiplexing_genome_only_regions' in crispresso2_info['running_info']['finished_steps'] and os.path.exists(filename_reads_aligned_to_genome_only):
                info('Using previously-computed extraction of aligned regions')
                df_regions = pd.read_csv(filename_reads_aligned_to_genome_only, sep="\t")
            else:
                info('Parsing the demultiplexed files and extracting locations and reference sequences...')
                files_to_match = list(df_all_demux['output filename'].dropna())
                summarize_region_fastq_input = [f+" "+uncompressed_reference for f in files_to_match] #pass both params to parallel function
                coordinates = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(summarize_region_fastq_input, summarize_region_fastq_chunk, n_processes=n_processes_for_pooled)
                df_regions = pd.DataFrame(coordinates, columns=['chr_id', 'bpstart', 'bpend', 'fastq_file', 'n_reads', 'sequence'])

                df_regions.dropna(inplace=True)  #remove regions in chrUn

                df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
                df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])
                df_regions['n_reads'] = pd.to_numeric(df_regions['n_reads'])

                df_regions.bpstart = df_regions.bpstart.astype(int)
                df_regions.bpend = df_regions.bpend.astype(int)

                df_regions['n_reads_aligned_%'] = df_regions['n_reads']/float(N_READS_ALIGNED)*100

                if args.gene_annotations:
                    info('Checking overlapping genes...')
                    df_regions = df_regions.apply(lambda row: find_overlapping_genes(row, df_genes), axis=1)

                df_regions.sort_values(by='n_reads', ascending=False, inplace=True)

                df_regions.fillna('NA').to_csv(filename_reads_aligned_to_genome_only, sep='\t', index=None)

                crispresso2_info['running_info']['finished_steps']['demultiplexing_genome_only_regions'] = True
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info,
                )

            # run CRISPResso (last step of genome-only mode)
            if can_finish_incomplete_run and 'crispresso_genome_only' in crispresso2_info['running_info']['finished_steps']:
                info('Using previously-computed crispresso runs')
            else:
                info('Running CRISPResso on the regions discovered...')
                crispresso_cmds = []
                for idx, row in df_regions.iterrows():

                    if row.n_reads > args.min_reads_to_use_region:
                        info('\nRunning CRISPResso on: %s-%d-%d...'%(row.chr_id, row.bpstart, row.bpend))
                        if pd.isna(row.sequence):
                            raise Exception('Cannot extract sequence from input reference ' + uncompressed_reference)
                        crispresso_cmd = args.crispresso_command + ' -r1 %s -a %s -o %s' %(row.fastq_file, row.sequence, OUTPUT_DIRECTORY)
                        crispresso_cmd = CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, crispresso_options_for_pooled, args)
                        crispresso_cmds.append(crispresso_cmd)
                    else:
                        info('Skipping region: %s-%d-%d , not enough reads (%d)' %(row.chr_id, row.bpstart, row.bpend, row.n_reads))
                CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, n_processes_for_pooled, 'region', args.skip_failed, start_end_percent=(15, 85))

                crispresso2_info['running_info']['finished_steps']['crispresso_genome_only'] = True
                CRISPRessoShared.write_crispresso_info(
                    crispresso2_info_file, crispresso2_info,
                )

        # write alignment statistics
        with open(_jp('MAPPING_STATISTICS.txt'), 'w+') as outfile:
            outfile.write('READS IN INPUTS:%d\nREADS AFTER PREPROCESSING:%d\nREADS ALIGNED:%d' % (N_READS_INPUT, N_READS_AFTER_PREPROCESSING, N_READS_ALIGNED))

        quantification_summary = []

        if RUNNING_MODE == 'ONLY_AMPLICONS' or RUNNING_MODE == 'AMPLICONS_AND_GENOME':
            df_final_data = df_template
        else:
            df_final_data = df_regions

        all_region_names = []
        all_region_read_counts = {}
        good_region_names = []
        good_region_folders = {}
        header = 'Name\tUnmodified%\tModified%\tReads_total\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions'
        header_els = header.split("\t")
        header_el_count = len(header_els)
        empty_line_els = [np.nan]*(header_el_count-1)
        n_reads_index = header_els.index('Reads_total') - 1
        for idx, row in df_final_data.iterrows():
                run_name = idx
                if RUNNING_MODE=='ONLY_AMPLICONS' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
                    run_name=idx
                else:
                    run_name='REGION_%s_%d_%d' %(row.chr_id, row.bpstart, row.bpend )
                folder_name = 'CRISPResso_on_%s'%run_name

                all_region_names.append(run_name)
                all_region_read_counts[run_name] = row.n_reads

                run_data = None
                try:
                    run_data = CRISPRessoShared.load_crispresso_info(_jp(folder_name))
                except:
                    warn('Skipping the folder %s: not enough reads, incomplete, or empty folder.'% folder_name)
                    this_els = empty_line_els[:]
                    this_els[n_reads_index] = row.n_reads
                    to_add = [run_name]
                    to_add.extend(this_els)
                    quantification_summary.append(to_add)
                else:
                    n_tot = row.n_reads
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
                    vals.extend([round(unmod_pct, 8), round(mod_pct, 8), n_tot, n_aligned, n_unmod, n_mod, n_discarded, n_insertion, n_deletion, n_substitution, n_only_insertion, n_only_deletion, n_only_substitution, n_insertion_and_deletion, n_insertion_and_substitution, n_deletion_and_substitution, n_insertion_and_deletion_and_substitution])
                    quantification_summary.append(vals)

                    good_region_names.append(run_name)
                    good_region_folders[run_name] = folder_name


        samples_quantification_summary_filename = _jp('SAMPLES_QUANTIFICATION_SUMMARY.txt')

        df_summary_quantification = pd.DataFrame(quantification_summary, columns=header_els)
        if args.crispresso1_mode:
            crispresso1_columns = ['amplicon_name', 'Unmodified%', 'Modified%', 'Reads_aligned', 'Reads_total']
            crispresso1_print_columns = ['Amplicon_Name', 'Unmodified%', 'Modified%', 'Reads_aligned', 'Reads_total']
            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename, sep='\t', index=None, columns=crispresso1_columns, header=crispresso1_print_columns)
        else:
            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename, sep='\t', index=None)

        crispresso2_info['results']['alignment_stats']['samples_quantification_summary_filename'] = os.path.basename(samples_quantification_summary_filename)
        crispresso2_info['results']['final_data'] = df_final_data
        crispresso2_info['results']['all_region_names'] = all_region_names
        crispresso2_info['results']['all_region_read_counts'] = all_region_read_counts
        crispresso2_info['results']['good_region_names'] = good_region_names
        crispresso2_info['results']['good_region_folders'] = good_region_folders
        crispresso2_info['running_info']['running_mode'] = RUNNING_MODE

        crispresso2_info['results']['general_plots']['summary_plot_names'] = []
        crispresso2_info['results']['general_plots']['summary_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_datas'] = {}


        df_summary_quantification.set_index('Name')

        save_png = True
        if args.suppress_report:
            save_png = False

        if not args.suppress_plots:
            plot_root = _jp("CRISPRessoPooled_reads_summary")

            debug('Plotting reads summary', {'percent_complete': 90})
            CRISPRessoPlot.plot_reads_total(plot_root, df_summary_quantification, save_png, args.min_reads_to_use_region)
            plot_name = os.path.basename(plot_root)
            crispresso2_info['results']['general_plots']['summary_plot_root'] = plot_name
            crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
            crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'CRISPRessoPooled Read Allocation Summary'
            crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads allocated to each amplicon. The vertical line shows the cutoff for analysis, set using the --min_reads_to_use_region parameter.'
            crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [('CRISPRessoPooled summary', os.path.basename(samples_quantification_summary_filename))]

            plot_root = _jp("CRISPRessoPooled_modification_summary")
            debug('Plotting modification summary', {'percent_complete': 95})
            CRISPRessoPlot.plot_unmod_mod_pcts(plot_root, df_summary_quantification, save_png, args.min_reads_to_use_region)
            plot_name = os.path.basename(plot_root)
            crispresso2_info['results']['general_plots']['summary_plot_root'] = plot_name
            crispresso2_info['results']['general_plots']['summary_plot_names'].append(plot_name)
            crispresso2_info['results']['general_plots']['summary_plot_titles'][plot_name] = 'CRISPRessoPooled Modification Summary'
            crispresso2_info['results']['general_plots']['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads aligned to each amplicon, divided into the reads that are modified and unmodified. The vertical line shows the cutoff for analysis, set using the --min_reads_to_use_region parameter.'
            crispresso2_info['results']['general_plots']['summary_plot_datas'][plot_name] = [('CRISPRessoPooled summary', os.path.basename(samples_quantification_summary_filename))]




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
                warn('Less than half (%d/%d) of reads aligned to amplicons. Finding most frequent unaligned reads.'%(tot_reads, N_READS_INPUT))
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
                p = sb.Popen(cmd, shell=True, stdout=sb.PIPE, preexec_fn=default_sigpipe)
                top_unaligned = p.communicate()[0].decode('utf-8')
                top_unaligned_filename=_jp('CRISPRessoPooled_TOP_UNALIGNED.txt')

                with open(top_unaligned_filename, 'w') as outfile:
                    outfile.write(top_unaligned)
                warn('Perhaps one or more of the given amplicon sequences were incomplete or incorrect. Below is a list of the most frequent unaligned reads (in the first 10000 unaligned reads). Check this list to see if an amplicon is among these reads.\n%s'%top_unaligned)


        #cleaning up
        if not args.keep_intermediate:
             info('Removing Intermediate files...')

             if not args.aligned_pooled_bam:
                if args.fastq_r2!='':
                    files_to_remove+=[processed_output_filename, flash_hist_filename, flash_histogram_filename,\
                                flash_not_combined_1_filename, flash_not_combined_2_filename]
                    if args.force_merge_pairs:
                        files_to_remove.append(new_merged_filename)
                        files_to_remove.append(old_flashed_filename)
                else:
                    files_to_remove+=[processed_output_filename]

             if args.trim_sequences and args.fastq_r2!='':
                 files_to_remove+=[output_forward_paired_filename, output_reverse_paired_filename,\
                                                   output_forward_unpaired_filename, output_reverse_unpaired_filename]

             if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
                 if args.aligned_pooled_bam is None:
                     files_to_remove+=[bam_filename_genome]
                     files_to_remove+=[bam_filename_genome+".bai"]

             if RUNNING_MODE=='ONLY_AMPLICONS':
                files_to_remove+=[bam_filename_amplicons, amplicon_fa_filename]
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
            CRISPRessoReport.make_pooled_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        if args.compile_postrun_references:
            postrun_references = []
            names_arr = crispresso2_info['results']['good_region_names']
            for name in names_arr:
                folder_name = 'CRISPResso_on_%s' % name
                sub_folder = os.path.join(OUTPUT_DIRECTORY, folder_name)
                run_data = None
                try:
                    run_data = CRISPRessoShared.load_crispresso_info(sub_folder)
                except Exception as e:
                    raise Exception('CRISPResso run %s is not complete. Cannot read CRISPResso2_info.json file.'% sub_folder)
                ref_sequences = [run_data['results']['refs'][ref_name]['sequence'] for ref_name in run_data['results']['ref_names']]
                allele_frequency_table_zip_filename = os.path.join(sub_folder, run_data['running_info']['allele_frequency_table_zip_filename'])
                if not os.path.exists(allele_frequency_table_zip_filename):
                    raise Exception('CRISPResso run %s is not complete. Cannot read allele frequency table.'% sub_folder)
                this_alleles = []
                this_freqs = []
                this_names = []
                with zipfile.ZipFile(allele_frequency_table_zip_filename, 'r') as archive:
                    with archive.open(run_data['running_info']['allele_frequency_table_filename'], 'r') as f:
                        head = f.readline().decode('UTF-8')
                        head_els = head.rstrip().split("\t")
                        allele_ind = head_els.index('Aligned_Sequence')
                        freq_ind = head_els.index('%Reads')

                        new_allele_idx = 1
                        for line in f:
                            line_els = line.decode('UTF-8').split('\t')
                            allele_seq = line_els[allele_ind].replace('-', '')
                            allele_freq = float(line_els[freq_ind])
                            #add first allele -- then add other alleles if they are more frequent than the cutoff
                            if len(this_alleles) > 0 and allele_freq < args.compile_postrun_reference_allele_cutoff:
                                break
                            if allele_seq not in this_alleles:
                                this_alleles.append(allele_seq)
                                this_freqs.append(allele_freq)
                                this_ref_name = ""
                                for ref_name in run_data['results']['ref_names']:
                                    if allele_seq == run_data['results']['refs'][ref_name]['sequence']:
                                        this_ref_name = ref_name
                                if this_ref_name == "":
                                    this_ref_name = "Alt_" + str(new_allele_idx)
                                    new_allele_idx += 1
                                this_names.append(this_ref_name)
                    postrun_references.append([name, ",".join(this_alleles), ",".join(this_names), ",".join([str(x) for x in this_freqs])])

            postrun_reference_file = _jp("CRISPRessoPooled_postrun_references.txt")
            pd.DataFrame(postrun_references, columns=['region_name', 'reference_seqs', 'reference_names', 'reference_frequencies']).to_csv(postrun_reference_file, sep="\t", index=False)
            crispresso2_info['postrun_reference_file'] = os.path.basename(postrun_reference_file)
            info('Produced postrun reference file: ' + postrun_reference_file)

        end_time =  datetime.now()
        end_time_string =  end_time.strftime('%Y-%m-%d %H:%M:%S')
        running_time = end_time - start_time
        running_time_string =  str(running_time)

        crispresso2_info['running_info']['end_time'] = end_time
        crispresso2_info['running_info']['end_time_string'] = end_time_string
        crispresso2_info['running_info']['running_time'] = running_time
        crispresso2_info['running_info']['running_time_string'] = running_time_string

        CRISPRessoShared.write_crispresso_info(
            crispresso2_info_file, crispresso2_info,
        )

        if args.zip_output:
            CRISPRessoShared.zip_results(OUTPUT_DIRECTORY)

        info('All Done!', {'percent_complete': 100})
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
