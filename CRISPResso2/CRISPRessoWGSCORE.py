# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''


import argparse
import gzip
import os
import re
import string
import subprocess as sb
import sys
import traceback
import unicodedata
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing


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
                error('You need to install %s module to use CRISPRessoWGS!' % library_name)
                sys.exit(1)

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference(set(['A','T','C','G','N'])))

def capitalize_sequence(x):
    return str(x).upper() if not pd.isnull(x) else x

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)

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
        sys.stdout.write('\nCRISPRessoWGS requires samtools')
        sys.stdout.write('\n\nPlease install it and add to your path following the instruction at: http://www.htslib.org/download/')
        return False


def check_bowtie2():

    cmd_path1=which('bowtie2')
    cmd_path2=which('bowtie2-inspect')

    if cmd_path1 and cmd_path2:
        return True
    else:
        sys.stdout.write('\nCRISPRessoWGS requires Bowtie2!')
        sys.stdout.write('\n\nPlease install it and add to your path following the instruction at: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2')
        return False

#if a reference index is provided aligne the reads to it
#extract region
def get_region_from_fa(chr_id,bpstart,bpend,uncompressed_reference):
    region='%s:%d-%d' % (chr_id,bpstart,bpend-1)
    p = sb.Popen("samtools faidx %s %s |   grep -v ^\> | tr -d '\n'" %(uncompressed_reference,region), shell=True,stdout=sb.PIPE)
    return p.communicate()[0]


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
        genes_overlapping.append( '%s (%s)' % (row_g.name2,row_g['name']))

    row['gene_overlapping']=','.join(genes_overlapping)

    return row


def find_last(mylist,myvalue):
    return len(mylist) - mylist[::-1].index(myvalue) -1

def get_reference_positions( pos, cigar,full_length=True):
    positions = []

    ops = re.findall(r'(\d+)(\w)', cigar)

    for c in ops:
        l,op=c
        l=int(l)

        if op == 'S' or op == 'I':
            if full_length:
                for i in range(0,l):
                    positions.append(None)
        elif op == 'M':
            for i in range(pos,pos+l):
                positions.append(i)
            pos += l
        elif op == 'D' or op == 'N':
            pos += l

    return positions



def write_trimmed_fastq(in_bam_filename,bpstart,bpend,out_fastq_filename):
    p = sb.Popen(
                'samtools view %s | cut -f1,4,6,10,11' % in_bam_filename,
                stdout = sb.PIPE,
                stderr = sb.STDOUT,
                shell=True
                )

    output=p.communicate()[0]
    n_reads=0

    with gzip.open(out_fastq_filename,'w+') as outfile:

        for line in output.split('\n'):
            if line:
                (name,pos,cigar,seq,qual)=line.split()
                #print name,pos,cigar,seq
                pos=int(pos)
                positions=get_reference_positions(pos,cigar)

                if bpstart in positions and bpend in positions:# and positions[0]<=bpstart and  positions[-1]>=bpend:

                    st=positions.index(bpstart)
                    en=find_last(positions,bpend)
                    #print st,en,seq,seq[st:en]
                    n_reads+=1
                    #print '>%s\n%s\n+\n%s\n' %(name,seq[st:en],qual[st:en])
                    outfile.write('@%s_%d\n%s\n+\n%s\n' %(name,n_reads,seq[st:en],qual[st:en]))
    return n_reads





pd=check_library('pandas')
np=check_library('numpy')

###EXCEPTIONS############################

class AmpliconsNamesNotUniqueException(Exception):
    pass

class SgRNASequenceException(Exception):
    pass

class NTException(Exception):
    pass

class ExonSequenceException(Exception):
    pass

def main():
    def print_stacktrace_if_debug():
        debug_flag = False
        if 'args' in vars() and 'debug' in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)
            error(traceback.format_exc())
    try:
        description = ['~~~CRISPRessoWGS~~~','-Analysis of CRISPR/Cas9 outcomes from WGS data-']
        wgs_string = r'''
 ____________
|     __  __ |
||  |/ _ (_  |
||/\|\__)__) |
|____________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description,wgs_string))

        parser = CRISPRessoShared.getCRISPRessoArgParser(parserTitle = 'CRISPRessoWGS Parameters',requiredParams={})

        #tool specific optional
        parser.add_argument('-b','--bam_file', type=str,  help='WGS aligned bam file', required=True,default='bam filename' )
        parser.add_argument('-f','--region_file', type=str,  help='Regions description file. A BED format  file containing the regions to analyze, one per line. The REQUIRED\
        columns are: chr_id(chromosome name), bpstart(start position), bpend(end position), the optional columns are:name (an unique indentifier for the region), guide_seq, expected_hdr_amplicon_seq,coding_seq, see CRISPResso help for more details on these last 3 parameters)', required=True)
        parser.add_argument('-r','--reference_file', type=str, help='A FASTA format reference file (for example hg19.fa for the human genome)', default='',required=True)
        parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=10)
        parser.add_argument('--skip_failed',  help='Continue with pooled analysis even if one sample fails',action='store_true')
        parser.add_argument('--gene_annotations', type=str, help='Gene Annotation Table from UCSC Genome Browser Tables (http://genome.ucsc.edu/cgi-bin/hgTables?command=start), \
        please select as table "knowGene", as output format "all fields from selected table" and as file returned "gzip compressed"', default='')
        parser.add_argument('-p','--n_processes',type=int, help='Specify the number of processes to use for the quantification.\
        Please use with caution since increasing this parameter will increase the memory required to run CRISPResso.',default=1)
        parser.add_argument('--crispresso_command', help='CRISPResso command to call',default='CRISPResso')

        args = parser.parse_args()

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = set(['fastq_r1','fastq_r2','amplicon_seq','amplicon_name','output_folder','name'])
        crispresso_options_for_wgs = list(crispresso_options-options_to_ignore)

        info('Checking dependencies...')

        if check_samtools() and check_bowtie2():
            info('\n All the required dependencies are present!')
        else:
            sys.exit(1)

        #check files
        check_file(args.bam_file)

        check_file(args.reference_file)

        check_file(args.region_file)

        if args.gene_annotations:
            check_file(args.gene_annotations)


        #INIT
        get_name_from_bam=lambda  x: os.path.basename(x).replace('.bam','')

        if not args.name:
            database_id='%s' % get_name_from_bam(args.bam_file)
        else:
            database_id=args.name


        OUTPUT_DIRECTORY='CRISPRessoWGS_on_%s' % database_id

        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory

        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoWGS_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\n%s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        #check if bam has the index already
        if os.path.exists(args.bam_file+'.bai'):
            info('Index file for input .bam file exists, skipping generation.')
        else:
            info('Creating index file for input .bam file...')
            sb.call('samtools index %s ' % (args.bam_file),shell=True)


        #load gene annotation
        if args.gene_annotations:
            info('Loading gene coordinates from annotation file: %s...' % args.gene_annotations)
            try:
                df_genes=pd.read_table(args.gene_annotations,compression='gzip')
                df_genes.txEnd=df_genes.txEnd.astype(int)
                df_genes.txStart=df_genes.txStart.astype(int)
                df_genes.head()
            except:
                info('Failed to load the gene annotations file.')


        #Load and validate the REGION FILE
        df_regions=pd.read_csv(args.region_file,names=[
                'chr_id','bpstart','bpend','Name','sgRNA',
                'Expected_HDR','Coding_sequence'],comment='#',sep='\t',dtype={'Name':str})


        #remove empty amplicons/lines
        df_regions.dropna(subset=['chr_id','bpstart','bpend'],inplace=True)

        df_regions.Expected_HDR=df_regions.Expected_HDR.apply(capitalize_sequence)
        df_regions.sgRNA=df_regions.sgRNA.apply(capitalize_sequence)
        df_regions.Coding_sequence=df_regions.Coding_sequence.apply(capitalize_sequence)


        #check or create names
        for idx,row in df_regions.iterrows():
            if pd.isnull(row.Name):
                df_regions.ix[idx,'Name']='_'.join(map(str,[row['chr_id'],row['bpstart'],row['bpend']]))


        if not len(df_regions.Name.unique())==df_regions.shape[0]:
            raise Exception('The amplicon names should be all distinct!')

        df_regions=df_regions.set_index('Name')
        #df_regions.index=df_regions.index.str.replace(' ','_')
        df_regions.index=df_regions.index.to_series().str.replace(' ','_')

        #extract sequence for each region
        uncompressed_reference=args.reference_file

        if os.path.exists(uncompressed_reference+'.fai'):
            info('The index for the reference fasta file is already present! Skipping generation.')
        else:
            info('Indexing reference file... Please be patient!')
            sb.call('samtools faidx %s >>%s 2>&1' % (uncompressed_reference,log_filename),shell=True)

        df_regions['sequence']=df_regions.apply(lambda row: get_region_from_fa(row.chr_id,row.bpstart,row.bpend,uncompressed_reference),axis=1)

        for idx,row in df_regions.iterrows():

            if not pd.isnull(row.sgRNA):

                cut_points=[]

                for current_guide_seq in row.sgRNA.strip().upper().split(','):

                    wrong_nt=find_wrong_nt(current_guide_seq)
                    if wrong_nt:
                        raise NTException('The sgRNA sequence %s contains wrong characters:%s'  % (current_guide_seq, ' '.join(wrong_nt)))

                    offset_fw=args.quantification_window_center+len(current_guide_seq)-1
                    offset_rc=(-args.quantification_window_center)-1
                    cut_points+=[m.start() + offset_fw for \
                                m in re.finditer(current_guide_seq,  row.sequence)]+[m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse_complement(current_guide_seq),  row.sequence)]

                if not cut_points:
                    df_regions.ix[idx,'sgRNA']=''


        df_regions=df_regions.convert_objects(convert_numeric=True)

        df_regions.bpstart=df_regions.bpstart.astype(int)
        df_regions.bpend=df_regions.bpend.astype(int)

        if args.gene_annotations:
            df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes),axis=1)


        #extract reads with samtools in that region and create a bam
        #create a fasta file with all the trimmed reads
        info('\nProcessing each regions...')

        ANALYZED_REGIONS=_jp('ANALYZED_REGIONS/')
        if not os.path.exists(ANALYZED_REGIONS):
            os.mkdir(ANALYZED_REGIONS)

        df_regions['n_reads']=0
        df_regions['bam_file_with_reads_in_region']=''
        df_regions['fastq.gz_file_trimmed_reads_in_region']=''

        for idx,row in df_regions.iterrows():

            if row['sequence']:

                fastq_gz_filename=os.path.join(ANALYZED_REGIONS,'%s.fastq.gz' % clean_filename('REGION_'+str(idx)))
                bam_region_filename=os.path.join(ANALYZED_REGIONS,'%s.bam' % clean_filename('REGION_'+str(idx)))

                #create place-holder fastq files
                open(fastq_gz_filename, 'w+').close()

                region='%s:%d-%d' % (row.chr_id,row.bpstart,row.bpend-1)
                info('\nExtracting reads in:%s and create the .bam file: %s' % (region,bam_region_filename))

                #extract reads in region
                cmd=r'''samtools view -b -F 4 %s %s > %s ''' % (args.bam_file, region, bam_region_filename)
                #print cmd
                sb.call(cmd,shell=True)


                #index bam file
                cmd=r'''samtools index %s ''' % (bam_region_filename)
                #print cmd
                sb.call(cmd,shell=True)

                info('Trim reads and create a fastq.gz file in: %s' % fastq_gz_filename)
                #trim reads in bam and convert in fastq
                n_reads=write_trimmed_fastq(bam_region_filename,row['bpstart'],row['bpend'],fastq_gz_filename)
                df_regions.ix[idx,'n_reads']=n_reads
                df_regions.ix[idx,'bam_file_with_reads_in_region']=bam_region_filename
                df_regions.ix[idx,'fastq.gz_file_trimmed_reads_in_region']=fastq_gz_filename


        df_regions.fillna('NA').to_csv(_jp('REPORT_READS_ALIGNED_TO_SELECTED_REGIONS_WGS.txt'),sep='\t')

        #Run Crispresso
        info('\nRunning CRISPResso on each regions...')
        crispresso_cmds = []
        for idx,row in df_regions.iterrows():

               if row['n_reads']>=args.min_reads_to_use_region:
                    info('\nThe region [%s] has enough reads (%d) mapped to it!' % (idx,row['n_reads']))

                    crispresso_cmd= args.crispresso_command + ' -r1 %s -a %s -o %s --name %s' %\
                    (row['fastq.gz_file_trimmed_reads_in_region'],row['sequence'],OUTPUT_DIRECTORY,idx)

                    if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                        crispresso_cmd+=' -g %s' % row['sgRNA']

                    if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                        crispresso_cmd+=' -e %s' % row['Expected_HDR']

                    if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                        crispresso_cmd+=' -c %s' % row['Coding_sequence']

                    crispresso_cmd=CRISPRessoShared.propagate_crispresso_options(crispresso_cmd,crispresso_options_for_wgs,args)
                    crispresso_cmds.append(crispresso_cmd)
#                    info('Running CRISPResso:%s' % crispresso_cmd)
#                    sb.call(crispresso_cmd,shell=True)

               else:
                    info('\nThe region [%s] has too few reads mapped to it (%d)! Skipping the running of CRISPResso!' % (idx,row['n_reads']))

        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,args.n_processes,'region',args.skip_failed)

        quantification_summary=[]

        for idx,row in df_regions.iterrows():

            folder_name='CRISPResso_on_%s' % idx
            try:
                quantification_file,amplicon_names,amplicon_info=CRISPRessoShared.check_output_folder(_jp(folder_name))
                for amplicon_name in amplicon_names:
                    N_TOTAL = float(amplicon_info[amplicon_name]['Total'])
                    N_UNMODIFIED = float(amplicon_info[amplicon_name]['Unmodified'])
                    N_MODIFIED = float(amplicon_info[amplicon_name]['Modified'])
                    quantification_summary.append([idx,amplicon_name,N_UNMODIFIED/N_TOTAL*100,N_MODIFIED/N_TOTAL*100,N_TOTAL,row.n_reads])
            except CRISPRessoShared.OutputFolderIncompleteException as e:
                quantification_summary.append([idx,"",np.nan,np.nan,np.nan,row.n_reads])
                warn('Skipping the folder %s: not enough reads, incomplete, or empty folder.'% folder_name)

        df_summary_quantification=pd.DataFrame(quantification_summary,columns=['Name','Amplicon','Unmodified%','Modified%','Reads_aligned','Reads_total'])
        df_summary_quantification.fillna('NA').to_csv(_jp('SAMPLES_QUANTIFICATION_SUMMARY.txt'),sep='\t',index=None)


        info('All Done!')
        print(CRISPRessoShared.get_crispresso_footer())
        sys.exit(0)

    except Exception as e:
        print_stacktrace_if_debug()
        error('\n\nERROR: %s' % e)
        sys.exit(-1)

if __name__ == '__main__':
    main()
