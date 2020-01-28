# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''


import os
import errno
import sys
from copy import deepcopy
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

#if a reference index is provided align the reads to it
#extract region
def get_region_from_fa(chr_id,bpstart,bpend,uncompressed_reference):
    region='%s:%d-%d' % (chr_id,bpstart,bpend-1)
    p = sb.Popen("samtools faidx %s %s |   grep -v ^\> | tr -d '\n'" %(uncompressed_reference,region), shell=True,stdout=sb.PIPE)
    return p.communicate()[0]

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
            genes_overlapping.append( '%s' % (row_g['#name']))
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
        parser.add_argument('-p','--n_processes',type=int, help='Specify the number of processes to use for Bowtie2.\
        Please use with caution since increasing this parameter will increase significantly the memory required to run CRISPResso.',default=1)
        parser.add_argument('-x','--bowtie2_index', type=str, help='Basename of Bowtie2 index for the reference genome', default='')
        parser.add_argument('--bowtie2_options_string', type=str, help='Override options for the Bowtie2 alignment command',default=' -k 1 --end-to-end -N 0 --np 0 ')
        parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=1000)
        parser.add_argument('--skip_failed',  help='Continue with pooled analysis even if one sample fails',action='store_true')
        parser.add_argument('--crispresso_command', help='CRISPResso command to call',default='CRISPResso')

        args = parser.parse_args()

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = set(['fastq_r1','fastq_r2','amplicon_seq','amplicon_name','output_folder','name'])
        crispresso_options_for_pooled = list(crispresso_options-options_to_ignore)


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

        crispresso2WGS_info_file = os.path.join(OUTPUT_DIRECTORY,'CRISPResso2Pooled_info.pickle')
        crispresso2_info = {} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['version'] = CRISPRessoShared.__version__
        crispresso2_info['args'] = deepcopy(args)

        crispresso2_info['log_filename'] = os.path.basename(log_filename)

        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

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





        #count reads
        N_READS_INPUT=get_n_reads_fastq(args.fastq_r1)
        N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(processed_output_filename)


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
            aligner_command= 'bowtie2 -x %s -p %s %s -U %s 2>>%s | samtools view -bS - > %s' %(custom_index_filename,args.n_processes,args.bowtie2_options_string,processed_output_filename,log_filename,bam_filename_amplicons)


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

            CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,args.n_processes,'amplicon',args.skip_failed)

            df_template['n_reads']=n_reads_aligned_amplicons
            df_template['n_reads_aligned_%']=df_template['n_reads']/float(N_READS_ALIGNED)*100
            df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_AMPLICONS.txt'),sep='\t')



        if RUNNING_MODE=='AMPLICONS_AND_GENOME':
            print 'Mapping amplicons to the reference genome...'
            #find the locations of the amplicons on the genome and their strand and check if there are mutations in the reference genome
            additional_columns=[]
            for idx,row in df_template.iterrows():
                fields_to_append=list(np.take(get_align_sequence(row.Amplicon_Sequence, args.bowtie2_index).split('\t'),[0,1,2,3,5]))
                if fields_to_append[0]=='*':
                    info('The amplicon [%s] is not mappable to the reference genome provided!' % idx )
                    additional_columns.append([idx,'NOT_ALIGNED',0,-1,'+',''])
                else:
                    additional_columns.append([idx]+fields_to_append)
                    info('The amplicon [%s] was mapped to: %s ' % (idx,' '.join(fields_to_append[:3]) ))


            df_template=df_template.join(pd.DataFrame(additional_columns,columns=['Name','chr_id','bpstart','bpend','strand','Reference_Sequence']).set_index('Name'))

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


        #####CORRECT ONE####
        #align in unbiased way the reads to the genome
        if RUNNING_MODE=='ONLY_GENOME' or RUNNING_MODE=='AMPLICONS_AND_GENOME':
            info('Aligning reads to the provided genome index...')
            bam_filename_genome = _jp('%s_GENOME_ALIGNED.bam' % database_id)
            aligner_command= 'bowtie2 -x %s -p %s %s -U %s 2>>%s| samtools view -bS - > %s' %(args.bowtie2_index,args.n_processes,args.bowtie2_options_string,processed_output_filename,log_filename,bam_filename_genome)
            info('aligning with command: ' + aligner_command)
            sb.call(aligner_command,shell=True)

            N_READS_ALIGNED=get_n_aligned_bam(bam_filename_genome)

            #REDISCOVER LOCATIONS and DEMULTIPLEX READS
            MAPPED_REGIONS=_jp('MAPPED_REGIONS/')
            if not os.path.exists(MAPPED_REGIONS):
                os.mkdir(MAPPED_REGIONS)

            s1=r'''samtools view -F 0x0004 %s 2>>%s |''' % (bam_filename_genome,log_filename)+\
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
            cmd=s1+s2.replace('__OUTPUTPATH__',MAPPED_REGIONS)

            info('Demultiplexing reads by location...')
            sb.call(cmd,shell=True)

            #gzip the missing ones
            sb.call('gzip -f %s/*.fastq' % MAPPED_REGIONS,shell=True)

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

            crispresso_cmds = []
            for idx,row in df_template.iterrows():

                info('Processing amplicon: %s' % idx )

                #check if we have reads
                fastq_filename_region=os.path.join(MAPPED_REGIONS,'REGION_%s_%s_%s.fastq.gz' % (row['chr_id'],row['bpstart'],row['bpend']))

                if os.path.exists(fastq_filename_region):

                    N_READS=get_n_reads_fastq(fastq_filename_region)
                    n_reads_aligned_genome.append(N_READS)
                    fastq_region_filenames.append(fastq_filename_region)
                    files_to_match.remove(fastq_filename_region)
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

            CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,args.n_processes,'amplicon',args.skip_failed)

            df_template['Amplicon_Specific_fastq.gz_filename']=fastq_region_filenames
            df_template['n_reads']=n_reads_aligned_genome
            df_template['n_reads_aligned_%']=df_template['n_reads']/float(N_READS_ALIGNED)*100

            if args.gene_annotations:
                df_template=df_template.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)

            df_template.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_GENOME_AND_AMPLICONS.txt'),sep='\t')

            #write another file with the not amplicon regions

            info('Reporting problematic regions...')
            coordinates=[]
            for region in files_to_match:
                #region format: REGION_chr8_1077_1198.fastq.gz
                #But if the chr has underscores, it could look like this:
                #    REGION_chr8_KI270812v1_alt_1077_1198.fastq.gz
                region_info = os.path.basename(region).replace('.fastq.gz','').replace('.fastq','').split('_')
                chr_string = "_".join(region_info[1:len(region_info)-2]) #in case there are underscores
                coordinates.append([chr_string] + region_info[-2:]+[region,get_n_reads_fastq(region)])

            df_regions=pd.DataFrame(coordinates,columns=['chr_id','bpstart','bpend','fastq_file','n_reads'])

            df_regions.dropna(inplace=True) #remove regions in chrUn

            df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
            df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])
            df_regions['n_reads'] = pd.to_numeric(df_regions['n_reads'])

            df_regions.bpstart=df_regions.bpstart.astype(int)
            df_regions.bpend=df_regions.bpend.astype(int)

            df_regions['n_reads_aligned_%']=df_regions['n_reads']/float(N_READS_ALIGNED)*100

            df_regions['Reference_sequence']=df_regions.apply(lambda row: get_region_from_fa(row.chr_id,row.bpstart,row.bpend,uncompressed_reference),axis=1)


            if args.gene_annotations:
                info('Checking overlapping genes...')
                df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)

            if np.sum(np.array(map(int,pd.__version__.split('.')))*(100,10,1))< 170:
                df_regions.sort('n_reads',ascending=False,inplace=True)
            else:
                df_regions.sort_values(by='n_reads',ascending=False,inplace=True)


            df_regions.fillna('NA').to_csv(_jp('REPORTS_READS_ALIGNED_TO_GENOME_NOT_MATCHING_AMPLICONS.txt'),sep='\t',index=None)


        if RUNNING_MODE=='ONLY_GENOME' :
            #Load regions and build REFERENCE TABLES
            info('Parsing the demultiplexed files and extracting locations and reference sequences...')
            coordinates=[]
            for region in glob.glob(os.path.join(MAPPED_REGIONS,'REGION*.fastq.gz')):
                coord_from_filename = os.path.basename(region).replace('.fastq.gz','').split('_')[1:4]
#                print('ccord from filename: ' + str(coord_from_filename))
                if not (coord_from_filename[1].isdigit() and coord_from_filename[2].isdigit()):
                    warn('Skipping region [%s] because the region name cannot be parsed\n'% region)
                    continue
                coordinates.append(coord_from_filename+[region,get_n_reads_fastq(region)])

            df_regions=pd.DataFrame(coordinates,columns=['chr_id','bpstart','bpend','fastq_file','n_reads'])

            df_regions.dropna(inplace=True) #remove regions in chrUn

            df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
            df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])
            df_regions['n_reads'] = pd.to_numeric(df_regions['n_reads'])

            df_regions.bpstart=df_regions.bpstart.astype(int)
            df_regions.bpend=df_regions.bpend.astype(int)
            df_regions['sequence']=df_regions.apply(lambda row: get_region_from_fa(row.chr_id,row.bpstart,row.bpend,uncompressed_reference),axis=1)

            df_regions['n_reads_aligned_%']=df_regions['n_reads']/float(N_READS_ALIGNED)*100

            if args.gene_annotations:
                info('Checking overlapping genes...')
                df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)

            if np.sum(np.array(map(int,pd.__version__.split('.')))*(100,10,1))< 170:
                df_regions.sort('n_reads',ascending=False,inplace=True)
            else:
                df_regions.sort_values(by='n_reads',ascending=False,inplace=True)


            df_regions.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_GENOME_ONLY.txt'),sep='\t',index=None)


            #run CRISPResso
            #demultiplex reads in the amplicons and call crispresso!
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
            CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,args.n_processes,'region',args.skip_failed)

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
    			warn('Less than half (%d/%d) of reads aligned. Finding most frequent unaligned reads.'%(tot_reads,N_READS_INPUT))
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
                 files_to_remove=[processed_output_filename,flash_hist_filename,flash_histogram_filename,\
                              flash_not_combined_1_filename,flash_not_combined_2_filename]
                 if args.force_merge_pairs:
                    files_to_remove.append(new_merged_filename)
                    files_to_remove.append(old_flashed_filename)
             else:
                 files_to_remove=[processed_output_filename]

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

        if not args.suppress_report:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2Pooled_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_pooled_report_from_folder(report_name,crispresso2_info,OUTPUT_DIRECTORY,_ROOT)
            crispresso2_info['report_location'] = report_name
            crispresso2_info['report_filename'] = os.path.basename(report_name)

        cp.dump(crispresso2_info, open(crispresso2WGS_info_file, 'wb' ) )

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
