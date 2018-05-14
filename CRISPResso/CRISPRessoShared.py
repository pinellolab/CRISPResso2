'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2017 The General Hospital Corporation. All Rights Reserved.
'''

import argparse
import os
import string
import pandas as pd
from collections import defaultdict
class CRISPRessoException(Exception):
    pass

__version__ = "2.0.04b"

##dict to lookup abbreviated params
crispresso_options_lookup = {
    'r1':'fastq_r1',
    'r2':'fastq_r2',
    'a':'amplicon_seq',
    'an':'amplicon_name',
    'amas':'amplicon_min_alignment_score',
    'g':'guide_seq',
    'c':'coding_seq',
    'q':'min_average_read_quality',
    's':'min_single_bp_quality',
    'n':'name',
    'o':'output_folder',
    'w':'window_around_sgrna',
    }


def getCRISPRessoArgParser(_ROOT, parserTitle = "CRISPResso Parameters",requiredParams={}):
    parser = argparse.ArgumentParser(description=parserTitle,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r1','--fastq_r1', type=str,  help='First fastq file',default='Fastq filename',required='fastq_r1' in requiredParams)
    parser.add_argument('-r2','--fastq_r2', type=str,  help='Second fastq file for paired end reads',default='')

    parser.add_argument('-a','--amplicon_seq', type=str,  help='Amplicon Sequence (can be comma-separated list of multiple sequences)', required='amplicon_seq' in requiredParams)

    parser.add_argument('-an','--amplicon_name', type=str,  help='Amplicon Name (can be comma-separated list of multiple names, corresponding to amplicon sequences given in --amplicon_seq', default='Amplicon')
    parser.add_argument('-amas','--amplicon_min_alignment_score', type=str,  help='Amplicon Minimum Alignment Score; score between 0 and 100; sequences must have at least this homology score with the amplicon to be aligned (can be comma-separated list of multiple scores, corresponding to amplicon sequences given in --amplicon_seq)', default="60")
    parser.add_argument('--default_min_aln_score','--min_identity_score',  type=int, help='Default minimum homology score for a read to align to a reference amplicon', default=60)
    parser.add_argument('--expand_ambiguous_alignments', help='If multiple reference amplicons are given, reads that align to multiple reference amplicons will count equally toward each amplicon. Default behavior is to exclude ambiguous alignments.', action='store_true')
    parser.add_argument('-g','--guide_seq',  help="sgRNA sequence, if more than one, please separate by comma/s. Note that the sgRNA needs to be input as the guide RNA sequence (usually 20 nt) immediately adjacent to but not including the PAM sequence (5' of NGG for SpCas9). If the PAM is found on the opposite strand with respect to the Amplicon Sequence, ensure the sgRNA sequence is also found on the opposite strand. The CRISPResso convention is to depict the expected cleavage position using the value of the parameter cleavage_offset nt  3' from the end of the guide. In addition, the use of alternate nucleases to SpCas9 is supported. For example, if using the Cpf1 system, enter the sequence (usually 20 nt) immediately 3' of the PAM sequence and explicitly set the cleavage_offset parameter to 1, since the default setting of -3 is suitable only for SpCas9.", default='')
    parser.add_argument('-c','--coding_seq',  help='Subsequence/s of the amplicon sequence covering one or more coding sequences for the frameshift analysis. If more than one (for example, split by intron/s), please separate by comma.', default='')
    parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
    parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
    parser.add_argument('--min_bp_quality_or_N', type=int, help='Bases with a quality score (phred33) less than this value will be set to "N"', default=0)
    parser.add_argument('-n','--name',  help='Output name', default='')
    parser.add_argument('-o','--output_folder',  help='', default='')

    ## read preprocessing params
    parser.add_argument('--split_paired_end',help='Splits a single fastq file containing paired end reads in two files before running CRISPResso',action='store_true')
    parser.add_argument('--trim_sequences',help='Enable the trimming of Illumina adapters with Trimmomatic',action='store_true')
    parser.add_argument('--trimmomatic_options_string', type=str, help='Override options for Trimmomatic',default=' ILLUMINACLIP:%s:0:90:10:0:true MINLEN:40' % os.path.join(_ROOT, 'data', 'NexteraPE-PE.fa'))
    parser.add_argument('--min_paired_end_reads_overlap',  type=int, help='Minimum required overlap length between two reads to provide a confident overlap. ', default=10)
    parser.add_argument('--left_adapter_umi_trim_seq', type=str, help='If deduplicating by umis, this is the end of the sequence on the left part of the read that should be trimmed after deduplicating')
    parser.add_argument('--right_adapter_umi_trim_seq', type=str, help='If deduplicating by umis, this is the start of the sequence on the right part of the read that should be trimmed after deduplicating')

    parser.add_argument('-w','--window_around_sgrna', type=int, help='Window(s) in bp around the cleavage position (half on on each side) as determined by the provide guide RNA sequence to quantify the indels. Any indels outside this window are excluded. A value of 0 disables this filter.', default=1)
    parser.add_argument('--cleavage_offset', type=int, help="Cleavage offset to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. The default is -3 and is suitable for the SpCas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1.", default=-3)
    parser.add_argument('--exclude_bp_from_left', type=int, help='Exclude bp from the left side of the amplicon sequence for the quantification of the indels', default=15)
    parser.add_argument('--exclude_bp_from_right', type=int, help='Exclude bp from the right side of the amplicon sequence for the quantification of the indels', default=15)

    parser.add_argument('--ignore_substitutions',help='Ignore substitutions events for the quantification and visualization',action='store_true')
    parser.add_argument('--ignore_insertions',help='Ignore insertions events for the quantification and visualization',action='store_true')
    parser.add_argument('--ignore_deletions',help='Ignore deletions events for the quantification and visualization',action='store_true')

    parser.add_argument('--needleman_wunsch_gap_open',type=int,help='Gap open option for Needleman-Wunsch alignment',default=-20)
    parser.add_argument('--needleman_wunsch_gap_extend',type=int,help='Gap extend option for Needleman-Wunsch alignment',default=-2)
    parser.add_argument('--needleman_wunsch_gap_incentive',type=int,help='Gap incentive value for inserting indels at cut sites',default=1)

    parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
    parser.add_argument('--dump',help='Dump numpy arrays and pandas dataframes to file for debugging purposes',action='store_true')
    parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')
    parser.add_argument('--offset_around_cut_to_plot',  type=int, help='Offset to use to summarize alleles around the cut site in the alleles table plot.', default=20)
    parser.add_argument('--min_frequency_alleles_around_cut_to_plot', type=float, help='Minimum %% reads required to report an allele in the alleles table plot.', default=0.2)
    parser.add_argument('--max_rows_alleles_around_cut_to_plot',  type=int, help='Maximum number of rows to report in the alleles table plot. ', default=50)

    parser.add_argument('--conversion_nuc_from',  help='For base editor plots, this is the nucleotide targeted by the base editor',default='C')
    parser.add_argument('--conversion_nuc_to',  help='For base editor plots, this is the nucleotide produced by the base editor',default='T')

    parser.add_argument('--base_editor_mode', help='Sets defaults for base editing experiments: analysis window is same size as sgRNA',action='store_true')
    parser.add_argument('-aw','--analysis_window_coordinates', type=str, help='Bp positions in the amplicon sequence specifying the analysis window. Any indels outside this window are excluded. Ranges are separted by the dash sign like "start-stop", and multiple ranges can be separated by the underscore (_). ' +
        'A value of 0 disables this filter. (can be comma-separated list of values, corresponding to amplicon sequences given in --amplicon_seq e.g. 5-10,5-10_20-30 would specify the 5th-10th bp in the first reference and the 5th-10th and 20th-30th bp in the second reference)', default=None)

    parser.add_argument('--crispresso1_mode', help='Parameter usage as in CRISPResso 1',action='store_true')
    parser.add_argument('--auto', help='Infer amplicon sequence from most common reads',action='store_true')
    parser.add_argument('--debug', help='Show debug messages', action='store_true')



    return parser

#######
# Nucleotide functions
#######
nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def reverse(seq):
    return "".join(c for c in seq.upper()[-1::-1])

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference(set(['A','T','C','G','N'])))

def capitalize_sequence(x):
    return str(x).upper() if not pd.isnull(x) else x


######
# File functions
######
def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise CRISPRessoException('The file: "'+filename+'" cannot be opened')


#get a clean name that we can use for a filename
validFilenameChars = "+-_.() %s%s" % (string.ascii_letters, string.digits)

def clean_filename(filename):
    cleanedFilename = unicodedata.normalize('NFKD', unicode(filename)).encode('ASCII', 'ignore')
    return ''.join(c for c in cleanedFilename if c in validFilenameChars)

def force_symlink(src, dst):

    if os.path.exists(dst) and os.path.samefile(src,dst):
        return

    try:
        os.symlink(src, dst)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)

def parse_count_file(fileName):
    if os.path.exists(fileName):
        with open(fileName) as infile:
            lines = infile.readlines()
            ampSeq = lines[0].rstrip().split("\t")
            ampSeq.pop(0) #get rid of 'Amplicon' at the beginning of line
            ampSeq = "".join(ampSeq)
            lab_freqs={}
            for i in range(1,len(lines)):
                line = lines[i].rstrip()
                lab_freq_arr = line.split()
                lab = lab_freq_arr.pop(0)
                lab_freqs[lab] = lab_freq_arr
        return ampSeq,lab_freqs
    else:
        print("Cannot find output file '%s'"%fileName)
        return None,None

def parse_alignment_file(fileName):
    if os.path.exists(fileName):
        with open(fileName) as infile:
            lines = infile.readlines()
            ampSeq = lines[0].rstrip().split("\t")
            ampSeq.pop(0) #get rid of 'Amplicon' at the beginning of line
            ampSeq = "".join(ampSeq)
            lab_freqs={}
            for i in range(1,len(lines)):
                line = lines[i].rstrip()
                lab_freq_arr = line.split()
                lab = lab_freq_arr.pop(0)
                lab_freqs[lab] = lab_freq_arr
        return ampSeq,lab_freqs
    else:
        print("Cannot find output file '%s'"%fileName)
        return None,None

def check_output_folder(output_folder):
    """
    Checks to see that the CRISPResso run has completed, and gathers the amplicon info for that run
    returns:
    - quantification file = CRISPResso_quantification_of_editing_frequency.txt for this run
    - amplicons = a list of amplicons analyzed in this run
    - amplicon_info = a dict of attributes found in quantification_file for each amplicon
    """
    quantification_file=os.path.join(output_folder,'CRISPResso_quantification_of_editing_frequency.txt')
    amplicon_info = {}
    amplicons = []
    all_files = os.listdir(output_folder)
    if os.path.exists(quantification_file):
        with open(quantification_file) as quant_file:
            head_line = quant_file.readline()
            head_line_els = head_line.split("\t")
            for line in quant_file:
                line_els = line.split("\t")
                amplicon_name = line_els[0]
                amplicon_quant_file = os.path.join(output_folder,amplicon_name + '.effect_vector_combined.txt')
                amplicon_info[amplicon_name] = {}
                if os.path.exists(amplicon_quant_file):
                    amplicons.append(amplicon_name)
                else:
                    raise OutputFolderIncompleteException('The folder %s  is not a valid CRISPResso2 output folder. Cannot find quantification file %s for amplicon %s.' % (output_folder,amplicon_quant_file,amplicon_name))
                amplicon_info[amplicon_name]['quantification_file'] = amplicon_quant_file

                amplicon_mod_count_file = os.path.join(output_folder,amplicon_name + '.target_modification_count_vectors.txt')
                if not os.path.exists(amplicon_mod_count_file):
                    raise OutputFolderIncompleteException('The folder %s  is not a valid CRISPResso2 output folder. Cannot find modification count vector file %s for amplicon %s.' % (output_folder,amplicon_mod_count_file,amplicon_name))
                amplicon_info[amplicon_name]['modification_count_file'] = amplicon_mod_count_file

                allele_files = []
                for file in all_files:
                    if re.match(amplicon_name + '.Alleles_frequency_table_around_cut_site_for_',file):
                        allele_files.append(os.path.join(output_folder,file))

                amplicon_info[amplicon_name]['allele_files'] = allele_files

                for idx,el in enumerate(head_line_els):
                    amplicon_info[amplicon_name][el] = line_els[idx]
        return quantification_file,amplicons,amplicon_info
    else:
        raise OutputFolderIncompleteException("The folder %s  is not a valid CRISPResso2 output folder. Cannot find quantification file '%s'." %(output_folder,quantification_file))

######
# allele modification functions
######

def get_row_around_cut(row,cut_point,offset):
    cut_idx=row['ref_positions'].index(cut_point)
    #cut_start = max(0,cut_idx-offset+1) #don't overflow the amplicon sequence
    #cut_end = min(len(row['ref_positions']),cut_idx+offset+1)
    #return row['Aligned_Sequence'][cut_start:cut_end],row['Reference_Sequence'][cut_start:cut_end],row['Read_Status']=='UNMODIFIED',row['n_deleted'],row['n_inserted'],row['n_mutated'],row['#Reads'], row['%Reads']
    #don't check overflow -- it was checked when program started
    return row['Aligned_Sequence'][cut_idx-offset+1:cut_idx+offset+1],row['Reference_Sequence'][cut_idx-offset+1:cut_idx+offset+1],row['Read_Status']=='UNMODIFIED',row['n_deleted'],row['n_inserted'],row['n_mutated'],row['#Reads'], row['%Reads']


def get_dataframe_around_cut(df_alleles, cut_point,offset):
    df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut(row,cut_point,offset),axis=1).values),
                        columns=['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated','#Reads','%Reads'])
    df_alleles_around_cut=df_alleles_around_cut.groupby(['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by='%Reads',inplace=True,ascending=False)
    df_alleles_around_cut['Unedited']=df_alleles_around_cut['Unedited']>0
    return df_alleles_around_cut


######
# misc functions
######
def get_crispresso_logo():
    return (r'''
     _
    '  )
    .-'
   (____
C)|     \
  \     /
   \___/
''')

def get_crispresso_header(description,header_str):
    """
    Creates the CRISPResso header string with the header_str between two crispresso mugs
    """
    term_width = 80

    logo = get_crispresso_logo()
    logo_lines = logo.splitlines()
    max_logo_width = max([len(x) for x in logo_lines])

    output_line = ""
    if header_str is not None:
        header_str = header_str.strip()

        header_lines = header_str.splitlines()
        while(len(header_lines) < len(logo_lines)):
            header_lines = [""] + header_lines

        max_header_width = max([len(x) for x in header_lines])


        pad_space = (term_width - (max_logo_width*2) - max_header_width)/4 - 1
        pad_string = " " * pad_space

        for i in range(len(logo_lines))[::-1]:
            output_line = (logo_lines[i].ljust(max_logo_width) + pad_string + header_lines[i].ljust(max_header_width) + pad_string + logo_lines[i].ljust(max_logo_width)).center(term_width) + "\n" + output_line

    else:
        pad_space = (term_width - max_logo_width)/2 - 1
        pad_string = " " * pad_space
        for i in range(len(logo_lines))[::-1]:
            output_line = (pad_string + logo_lines[i].ljust(max_logo_width) + pad_string).center(term_width) + "\n" + output_line

    output_line += '\n'+('[CRISPresso version ' + __version__ + ']').center(term_width) + '\n' + ('[Kendell Clement and Luca Pinello 2018]').center(term_width) + "\n" + ('[For support, contact kclement@mgh.harvard.edu]').center(term_width) + "\n"

    description_str = ""
    for str in description:
        str = str.strip()
        description_str += str.center(term_width) + "\n"

    return "\n" + description_str + output_line

def get_crispresso_footer():
    logo = get_crispresso_logo()
    logo_lines = logo.splitlines()

    max_logo_width = max([len(x) for x in logo_lines])
    pad_space = (80 - (max_logo_width))/2 - 1
    pad_string = " " * pad_space

    output_line = ""
    for i in range(len(logo_lines))[::-1]:
        output_line = pad_string + logo_lines[i].ljust(max_logo_width) + pad_string + "\n" + output_line

    return output_line
