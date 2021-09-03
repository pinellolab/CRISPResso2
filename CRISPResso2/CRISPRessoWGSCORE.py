# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''


import argparse
from datetime import datetime
import gzip
import os
import re
from copy import deepcopy
import string
import subprocess as sb
import sys
import traceback
import unicodedata
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoReport
from CRISPResso2 import CRISPRessoPlot


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
    return list(set(sequence.upper()).difference({'A', 'T', 'C', 'G', 'N'}))

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
        sys.stdout.write('\n\nPlease install samtools and add it to your path following the instructions at: http://www.htslib.org/download/')
        return False


def check_bowtie2():

    cmd_path1=which('bowtie2')
    cmd_path2=which('bowtie2-inspect')

    if cmd_path1 and cmd_path2:
        return True
    else:
        sys.stdout.write('\nCRISPRessoWGS requires Bowtie2!')
        sys.stdout.write('\n\nPlease install Bowtie2 and add it to your path following the instructions at: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2')
        return False

#if a reference index is provided aligne the reads to it
#extract region
def get_region_from_fa(chr_id, bpstart, bpend, uncompressed_reference):
    region='%s:%d-%d' % (chr_id, bpstart, bpend-1)
    p = sb.Popen("samtools faidx %s %s |   grep -v ^\> | tr -d '\n'" %(uncompressed_reference, region), shell=True, stdout=sb.PIPE)
    return p.communicate()[0].decode('utf-8').upper()


#get a clean name that we can use for a filename
validFilenameChars = "+-_.() %s%s" % (string.ascii_letters, string.digits)

def clean_filename(filename):
    cleanedFilename = unicodedata.normalize('NFKD', filename)
    return ''.join(c for c in cleanedFilename if c in validFilenameChars)


def find_overlapping_genes(row, df_genes):
    df_genes_overlapping=df_genes.ix[(df_genes.chrom==row.chr_id) &
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


def find_last(mylist, myvalue):
    return len(mylist) - mylist[::-1].index(myvalue) -1

def get_reference_positions( pos, cigar,full_length=True):
    positions = []

    ops = re.findall(r'(\d+)(\w)', cigar)

    for c in ops:
        l, op=c
        l=int(l)

        if op == 'S' or op == 'I':
            if full_length:
                for i in range(0, l):
                    positions.append(None)
        elif op == 'M':
            for i in range(pos, pos+l):
                positions.append(i)
            pos += l
        elif op == 'D' or op == 'N':
            pos += l

    return positions

def write_trimmed_fastq(in_bam_filename, bpstart, bpend, out_fastq_filename):
    p = sb.Popen(
                'samtools view %s | cut -f1,4,6,10,11' % in_bam_filename,
                stdout = sb.PIPE,
                stderr = sb.STDOUT,
                shell=True
                )

    output=p.communicate()[0].decode('utf-8')
    n_reads=0

    with gzip.open(out_fastq_filename, 'wt') as outfile:
        for line in output.split('\n'):
            if line:
                (name, pos, cigar, seq, qual)=line.split()
                #print name,pos,cigar,seq
                pos=int(pos)
                positions=get_reference_positions(pos, cigar)

                if bpstart in positions and bpend in positions:# and positions[0]<=bpstart and  positions[-1]>=bpend:

                    st=positions.index(bpstart)
                    en=find_last(positions, bpend)
                    #print st,en,seq,seq[st:en]
                    n_reads+=1
                    #print '>%s\n%s\n+\n%s\n' %(name,seq[st:en],qual[st:en])
                    outfile.write('@%s_%d\n%s\n+\n%s\n' %(name, n_reads, seq[st:en], qual[st:en]))
    return n_reads

pd=check_library('pandas')
np=check_library('numpy')

def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename, shell=True, stdout=sb.PIPE)
     n_reads = int(float(p.communicate()[0])/4.0)
     return n_reads

def extract_reads(row):
    if row.sequence:
        #create place-holder fastq files
        open(row.fastq_file_trimmed_reads_in_region, 'w+').close()

        region='%s:%d-%d' % (row.chr_id, row.bpstart, row.bpend-1)

        info('Extracting reads in:%s and creating .bam file: %s' % (region, row.bam_file_with_reads_in_region))

        cmd=r'''samtools view -b -F 4 %s %s > %s ''' % (row.original_bam, region, row.bam_file_with_reads_in_region)
        sb.call(cmd, shell=True)

        cmd=r'''samtools index %s ''' % (row.bam_file_with_reads_in_region)
        sb.call(cmd, shell=True)

        #trim reads in bam and convert in fastq
        row.n_reads=write_trimmed_fastq(row.bam_file_with_reads_in_region, row.bpstart, row.bpend, row.fastq_file_trimmed_reads_in_region)
    else:
        row.n_reads = 0
        row.bam_file_with_reads_in_region = ''
        row.fastq_file_trimmed_reads_in_region = ''

    return row

def extract_reads_chunk(df):
    new_df = pd.DataFrame(columns=df.columns)
    for i in range(len(df)):
        new_df = new_df.append(extract_reads(df.iloc[i].copy()))
    return(new_df)


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
        start_time =  datetime.now()
        start_time_string =  start_time.strftime('%Y-%m-%d %H:%M:%S')

        description = ['~~~CRISPRessoWGS~~~', '-Analysis of CRISPR/Cas9 outcomes from WGS data-']
        wgs_string = r'''
 ____________
|     __  __ |
||  |/ _ (_  |
||/\|\__)__) |
|____________|
        '''
        print(CRISPRessoShared.get_crispresso_header(description, wgs_string))

        parser = CRISPRessoShared.getCRISPRessoArgParser(parserTitle = 'CRISPRessoWGS Parameters', requiredParams={})

        #tool specific optional
        parser.add_argument('-b', '--bam_file', type=str,  help='WGS aligned bam file', required=True, default='bam filename' )
        parser.add_argument('-f', '--region_file', type=str,  help='Regions description file. A BED format  file containing the regions to analyze, one per line. The REQUIRED\
        columns are: chr_id(chromosome name), bpstart(start position), bpend(end position), the optional columns are:name (an unique indentifier for the region), guide_seq, expected_hdr_amplicon_seq,coding_seq, see CRISPResso help for more details on these last 3 parameters)', required=True)
        parser.add_argument('-r', '--reference_file', type=str, help='A FASTA format reference file (for example hg19.fa for the human genome)', default='', required=True)
        parser.add_argument('--min_reads_to_use_region',  type=float, help='Minimum number of reads that align to a region to perform the CRISPResso analysis', default=10)
        parser.add_argument('--skip_failed',  help='Continue with pooled analysis even if one sample fails', action='store_true')
        parser.add_argument('--gene_annotations', type=str, help='Gene Annotation Table from UCSC Genome Browser Tables (http://genome.ucsc.edu/cgi-bin/hgTables?command=start), \
        please select as table "knownGene", as output format "all fields from selected table" and as file returned "gzip compressed"', default='')
        parser.add_argument('--crispresso_command', help='CRISPResso command to call', default='CRISPResso')

        args = parser.parse_args()

        crispresso_options = CRISPRessoShared.get_crispresso_options()
        options_to_ignore = {'fastq_r1', 'fastq_r2', 'amplicon_seq', 'amplicon_name', 'output_folder', 'name'}
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

        # for computation performed in CRISPressoPooled (e.g. bowtie alignment, etc) use n_processes_for_wgs
        n_processes_for_wgs = 1
        if args.n_processes == "max":
            n_processes_for_wgs = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes_for_wgs = int(args.n_processes)

        # here, we set args.n_processes as another value because this value is propagated to sub-CRISPResso runs (not for usage in CRISPRessoWGS)
        args.n_processes = CRISPRessoShared.get_sub_n_processes(suppress_plots=args.suppress_plots, suppress_report=args.suppress_report, n_processes=args.n_processes)

        #INIT
        get_name_from_bam=lambda  x: os.path.basename(x).replace('.bam', '')

        if not args.name:
            database_id='%s' % get_name_from_bam(args.bam_file)
        else:
            clean_name = CRISPRessoShared.slugify(args.name)
            if args.name != clean_name:
                warn(
                     'The specified name {0} contained invalid characters and was changed to: {1}'.format(
                         args.name, clean_name,
                    ),
                )
            database_id = clean_name


        OUTPUT_DIRECTORY='CRISPRessoWGS_on_%s' % database_id

        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory

        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoWGS_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        crispresso2_info_file = os.path.join(OUTPUT_DIRECTORY, 'CRISPResso2WGS_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)
        crispresso2_info['running_info']['finished_steps'] = {}


        crispresso_cmd_to_write = ' '.join(sys.argv)
        if args.write_cleaned_report:
            cmd_copy = sys.argv[:]
            cmd_copy[0] = 'CRISPRessoWGS'
            for i in range(len(cmd_copy)):
                if os.sep in cmd_copy[i]:
                    cmd_copy[i] = os.path.basename(cmd_copy[i])

            crispresso_cmd_to_write = ' '.join(cmd_copy) #clean command doesn't show the absolute path to the executable or other files
        crispresso2_info['running_info']['command_used'] = crispresso_cmd_to_write

        with open(log_filename, 'w+') as outfile:
            outfile.write('CRISPResso version %s\n[Command used]:\n%s\n\n[Execution log]:\n' %(CRISPRessoShared.__version__, crispresso_cmd_to_write))

        #keep track of args to see if it is possible to skip computation steps on rerun
        can_finish_incomplete_run = False
        if args.no_rerun:
            if os.path.exists(crispresso2_info_file):
                previous_run_data = CRISPRessoShared.load_crispresso_info(OUTPUT_DIRECTORY)
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
                            if 'finished_steps' in previous_run_data['running_info']:
                                for key in previous_run_data['running_info']['finished_steps'].keys():
                                    crispresso2_info['running_info']['finished_steps'][key] = previous_run_data['running_info']['finished_steps'][key]
                                    if args.debug:
                                        info('finished: ' + key)
                else:
                    info('The no_rerun flag is set, but this analysis will be rerun because the existing run was performed using an old version of CRISPResso (' + str(previous_run_data['running_info']['version']) + ').')

        #write this file early on so we can check the params if we have to rerun
        CRISPRessoShared.write_crispresso_info(
            crispresso2_info_file, crispresso2_info,
        )

        def rreplace(s, old, new):
            li = s.rsplit(old)
            return new.join(li)

        bam_index = ''
        #check if bam has the index already
        if os.path.exists(rreplace(args.bam_file, ".bam", ".bai")):
            info('Index file for input .bam file exists, skipping generation.')
            bam_index = args.bam_file.replace(".bam", ".bai")
        elif os.path.exists(args.bam_file+'.bai'):
            info('Index file for input .bam file exists, skipping generation.')
            bam_index = args.bam_file+'.bai'
        else:
            info('Creating index file for input .bam file...')
            sb.call('samtools index %s ' % (args.bam_file), shell=True)
            bam_index = args.bam_file+'.bai'


        #load gene annotation
        if args.gene_annotations:
            info('Loading gene coordinates from annotation file: %s...' % args.gene_annotations)
            try:
                df_genes=pd.read_csv(args.gene_annotations, compression='gzip', sep="\t")
                df_genes.txEnd=df_genes.txEnd.astype(int)
                df_genes.txStart=df_genes.txStart.astype(int)
                df_genes.head()
            except:
                raise Exception('Failed to load the gene annotations file.')


        #Load and validate the REGION FILE
        df_regions=pd.read_csv(args.region_file, names=[
                'chr_id', 'bpstart', 'bpend', 'Name', 'sgRNA',
                'Expected_HDR', 'Coding_sequence'], comment='#', sep='\t', dtype={'Name':str, 'chr_id':str})


        #remove empty amplicons/lines
        df_regions.dropna(subset=['chr_id', 'bpstart', 'bpend'], inplace=True)

        df_regions.Expected_HDR=df_regions.Expected_HDR.apply(capitalize_sequence)
        df_regions.sgRNA=df_regions.sgRNA.apply(capitalize_sequence)
        df_regions.Coding_sequence=df_regions.Coding_sequence.apply(capitalize_sequence)


        #check or create names
        for idx, row in df_regions.iterrows():
            if pd.isnull(row.Name):
                df_regions.ix[idx, 'Name']='_'.join(map(str, [row['chr_id'], row['bpstart'], row['bpend']]))


        if not len(df_regions.Name.unique())==df_regions.shape[0]:
            raise Exception('The amplicon names should be all distinct!')

        df_regions.set_index('Name', inplace=True)
        #df_regions.index=df_regions.index.str.replace(' ','_')
        df_regions.index=df_regions.index.to_series().str.replace(' ', '_')

        #extract sequence for each region
        uncompressed_reference=args.reference_file

        if os.path.exists(uncompressed_reference+'.fai'):
            info('The index for the reference fasta file is already present! Skipping generation.')
        else:
            info('Indexing reference file... Please be patient!')
            sb.call('samtools faidx %s >>%s 2>&1' % (uncompressed_reference, log_filename), shell=True)

        info('Retrieving reference sequences for amplicons and checking for sgRNAs')
        df_regions['sequence']=df_regions.apply(lambda row: get_region_from_fa(row.chr_id, row.bpstart, row.bpend, uncompressed_reference), axis=1)

        for idx, row in df_regions.iterrows():

            if not pd.isnull(row.sgRNA):

                cut_points=[]
                guides = row.sgRNA.strip().upper().split(',')
                guide_qw_centers = CRISPRessoShared.set_guide_array(args.quantification_window_center, guides, 'guide quantification center')
                for idx, current_guide_seq in enumerate(guides):

                    wrong_nt=find_wrong_nt(current_guide_seq)
                    if wrong_nt:
                        raise NTException('The sgRNA sequence %s contains wrong characters:%s'  % (current_guide_seq, ' '.join(wrong_nt)))

                    offset_fw=guide_qw_centers[idx]+len(current_guide_seq)-1
                    offset_rc=(-guide_qw_centers[idx])-1
                    cut_points+=[m.start() + offset_fw for \
                                m in re.finditer(current_guide_seq,  row.sequence)]+[m.start() + offset_rc for m in re.finditer(CRISPRessoShared.reverse_complement(current_guide_seq),  row.sequence)]

                if not cut_points:
                    df_regions.ix[idx, 'sgRNA']=''
                    info('Cannot find guide ' + str(row.sgRNA) + ' in amplicon ' + str(idx) + ' (' + str(row) + ')')

        df_regions['bpstart'] = pd.to_numeric(df_regions['bpstart'])
        df_regions['bpend'] = pd.to_numeric(df_regions['bpend'])

        df_regions.bpstart=df_regions.bpstart.astype(int)
        df_regions.bpend=df_regions.bpend.astype(int)

        if args.gene_annotations:
            df_regions=df_regions.apply(lambda row: find_overlapping_genes(row, df_genes), axis=1)


        #extract reads with samtools in that region and create a bam
        #create a fasta file with all the trimmed reads
        info('\nProcessing each region...')

        ANALYZED_REGIONS=_jp('ANALYZED_REGIONS/')
        if not os.path.exists(ANALYZED_REGIONS):
            os.mkdir(ANALYZED_REGIONS)


        df_regions['region_number'] = np.arange(len(df_regions))

        def set_filenames(row):
            row_fastq_exists = False
            fastq_gz_filename=os.path.join(ANALYZED_REGIONS, '%s.fastq.gz' % clean_filename('REGION_'+str(row.region_number)))
            bam_region_filename=os.path.join(ANALYZED_REGIONS, '%s.bam' % clean_filename('REGION_'+str(row.region_number)))
            #if bam file already exists, don't regenerate it
            if os.path.isfile(fastq_gz_filename):
                row_fastq_exists = True
            return bam_region_filename, fastq_gz_filename, row_fastq_exists

        df_regions['bam_file_with_reads_in_region'], df_regions['fastq_file_trimmed_reads_in_region'], df_regions['row_fastq_exists'] = zip(*df_regions.apply(set_filenames, axis=1))
        df_regions['n_reads'] = 0
        df_regions['original_bam'] = args.bam_file #stick this in the df so we can parallelize the analysis and not pass params


        report_reads_aligned_filename = _jp('REPORT_READS_ALIGNED_TO_SELECTED_REGIONS_WGS.txt')
        num_rows_without_fastq = len(df_regions[df_regions.row_fastq_exists == False])

        if can_finish_incomplete_run and num_rows_without_fastq == 0 and os.path.isfile(report_reads_aligned_filename) and 'generation_of_fastq_files_for_each_amplicon' in crispresso2_info['running_info']['finished_steps']:
            info('Skipping generation of fastq files for each amplicon.')
            df_regions = pd.read_csv(report_reads_aligned_filename, comment='#', sep='\t', dtype={'Name':str, 'chr_id':str})
            df_regions.set_index('Name', inplace=True)

        else:
            #run region extraction here
            df_regions = CRISPRessoMultiProcessing.run_pandas_apply_parallel(df_regions, extract_reads_chunk, n_processes_for_wgs)
            df_regions.sort_values('region_number', inplace=True)
            cols_to_print = ["chr_id", "bpstart", "bpend", "sgRNA", "Expected_HDR", "Coding_sequence", "sequence", "n_reads", "bam_file_with_reads_in_region", "fastq_file_trimmed_reads_in_region"]
            if args.gene_annotations:
                cols_to_print.append('gene_overlapping')
            df_regions.fillna('NA').to_csv(report_reads_aligned_filename, sep='\t', columns = cols_to_print, index_label="Name")

            #save progress
            crispresso2_info['running_info']['finished_steps']['generation_of_fastq_files_for_each_amplicon'] = True
            CRISPRessoShared.write_crispresso_info(
                crispresso2_info_file, crispresso2_info,
            )

        #Run Crispresso
        info('Running CRISPResso on each region...')
        crispresso_cmds = []
        for idx, row in df_regions.iterrows():
            if row['n_reads']>=args.min_reads_to_use_region:
                info('\nThe region [%s] has enough reads (%d) mapped to it!' % (idx, row['n_reads']))

                crispresso_cmd= args.crispresso_command + ' -r1 %s -a %s -o %s --name %s' %\
                (row['fastq_file_trimmed_reads_in_region'], row['sequence'], OUTPUT_DIRECTORY, idx)

                if row['sgRNA'] and not pd.isnull(row['sgRNA']):
                    crispresso_cmd+=' -g %s' % row['sgRNA']

                if row['Expected_HDR'] and not pd.isnull(row['Expected_HDR']):
                    crispresso_cmd+=' -e %s' % row['Expected_HDR']

                if row['Coding_sequence'] and not pd.isnull(row['Coding_sequence']):
                    crispresso_cmd+=' -c %s' % row['Coding_sequence']

                crispresso_cmd=CRISPRessoShared.propagate_crispresso_options(crispresso_cmd, crispresso_options_for_wgs, args)

                #logging like this causes the multiprocessing step to not block for some reason #mysteriesOfThPythonUniverse
                #log_name = _jp("CRISPResso_on_"+idx) +".log"
                #crispresso_cmd += " &> %s"%log_name

                crispresso_cmds.append(crispresso_cmd)
#                    info('Running CRISPResso:%s' % crispresso_cmd)
#                    sb.call(crispresso_cmd,shell=True)

            else:
                info('\nThe region [%s] has too few reads mapped to it (%d)! Not running CRISPResso!' % (idx, row['n_reads']))

        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds, args.n_processes, 'region', args.skip_failed)

        quantification_summary=[]
        all_region_names = []
        all_region_read_counts = {}
        good_region_names = []
        good_region_folders = {}
        header = 'Name\tUnmodified%\tModified%\tReads_total\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions'
        header_els = header.split("\t")
        header_el_count = len(header_els)
        empty_line_els = [np.nan]*(header_el_count-1)
        n_reads_index = header_els.index('Reads_total') - 1
        for idx, row in df_regions.iterrows():
            folder_name='CRISPResso_on_%s' % idx
            run_name = idx

            all_region_names.append(run_name)
            all_region_read_counts[run_name] = row.n_reads

            run_file = os.path.join(_jp(folder_name), 'CRISPResso2_info.json')
            if not os.path.exists(run_file):
                warn('Skipping the folder %s: not enough reads, incomplete, or empty folder.'% folder_name)
                this_els = empty_line_els[:]
                this_els[n_reads_index] = row.n_reads
                to_add = [run_name]
                to_add.extend(this_els)
                quantification_summary.append(to_add)
            else:
                run_data = CRISPRessoShared.load_crispresso_info(
                    _jp(folder_name),
                )
                ref_name = run_data['results']['ref_names'][0] #only expect one amplicon sequence
                n_tot = row.n_reads
                n_aligned = run_data['results']['alignment_stats']['counts_total'][ref_name]
                n_unmod = run_data['results']['alignment_stats']['counts_unmodified'][ref_name]
                n_mod = run_data['results']['alignment_stats']['counts_modified'][ref_name]
                n_discarded = run_data['results']['alignment_stats']['counts_discarded'][ref_name]

                n_insertion = run_data['results']['alignment_stats']['counts_insertion'][ref_name]
                n_deletion = run_data['results']['alignment_stats']['counts_deletion'][ref_name]
                n_substitution = run_data['results']['alignment_stats']['counts_substitution'][ref_name]
                n_only_insertion = run_data['results']['alignment_stats']['counts_only_insertion'][ref_name]
                n_only_deletion = run_data['results']['alignment_stats']['counts_only_deletion'][ref_name]
                n_only_substitution = run_data['results']['alignment_stats']['counts_only_substitution'][ref_name]
                n_insertion_and_deletion = run_data['results']['alignment_stats']['counts_insertion_and_deletion'][ref_name]
                n_insertion_and_substitution = run_data['results']['alignment_stats']['counts_insertion_and_substitution'][ref_name]
                n_deletion_and_substitution = run_data['results']['alignment_stats']['counts_deletion_and_substitution'][ref_name]
                n_insertion_and_deletion_and_substitution = run_data['results']['alignment_stats']['counts_insertion_and_deletion_and_substitution'][ref_name]

                unmod_pct = "NA"
                mod_pct = "NA"
                if n_aligned > 0:
                    unmod_pct = 100*n_unmod/float(n_aligned)
                    mod_pct = 100*n_mod/float(n_aligned)

                vals = [run_name]
                vals.extend([round(unmod_pct, 8), round(mod_pct, 8), n_aligned, n_tot, n_unmod, n_mod, n_discarded, n_insertion, n_deletion, n_substitution, n_only_insertion, n_only_deletion, n_only_substitution, n_insertion_and_deletion, n_insertion_and_substitution, n_deletion_and_substitution, n_insertion_and_deletion_and_substitution])
                quantification_summary.append(vals)

                good_region_names.append(idx)
                good_region_folders[idx] = folder_name
        samples_quantification_summary_filename = _jp('SAMPLES_QUANTIFICATION_SUMMARY.txt')

        df_summary_quantification=pd.DataFrame(quantification_summary, columns=header_els)
        if args.crispresso1_mode:
            crispresso1_columns=['Name', 'Unmodified%', 'Modified%', 'Reads_aligned', 'Reads_total']
            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename, sep='\t', index=None, columns=crispresso1_columns)
        else:
            df_summary_quantification.fillna('NA').to_csv(samples_quantification_summary_filename, sep='\t', index=None)

        crispresso2_info['results']['alignment_stats']['samples_quantification_summary_filename'] = os.path.basename(samples_quantification_summary_filename)
        crispresso2_info['results']['regions'] = df_regions
        crispresso2_info['results']['all_region_names'] = all_region_names
        crispresso2_info['results']['all_region_read_counts'] = all_region_read_counts
        crispresso2_info['results']['good_region_names'] = good_region_names
        crispresso2_info['results']['good_region_folders'] = good_region_folders

        crispresso2_info['summary_plot_names'] = []
        crispresso2_info['summary_plot_titles'] = {}
        crispresso2_info['summary_plot_labels'] = {}
        crispresso2_info['summary_plot_datas'] = {}

        df_summary_quantification.set_index('Name')

        save_png = True
        if args.suppress_report:
            save_png = False

        if not args.suppress_plots:
            plot_root = _jp("CRISPRessoWGS_reads_summary")
            CRISPRessoPlot.plot_reads_total(plot_root, df_summary_quantification, save_png, args.min_reads_to_use_region)
            plot_name = os.path.basename(plot_root)
            crispresso2_info['reads_summary_plot'] = plot_name
            crispresso2_info['summary_plot_names'].append(plot_name)
            crispresso2_info['summary_plot_titles'][plot_name] = 'CRISPRessoWGS Read Allocation Summary'
            crispresso2_info['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads allocated to each amplicon. The vertical line shows the cutoff for analysis, set using the --min_reads_to_use_region parameter.'
            crispresso2_info['summary_plot_datas'][plot_name] = [('CRISPRessoWGS summary', os.path.basename(samples_quantification_summary_filename))]

            plot_root = _jp("CRISPRessoWGS_modification_summary")
            CRISPRessoPlot.plot_unmod_mod_pcts(plot_root, df_summary_quantification, save_png, args.min_reads_to_use_region)
            plot_name = os.path.basename(plot_root)
            crispresso2_info['modification_summary_plot'] = plot_name
            crispresso2_info['summary_plot_names'].append(plot_name)
            crispresso2_info['summary_plot_titles'][plot_name] = 'CRISPRessoWGS Modification Summary'
            crispresso2_info['summary_plot_labels'][plot_name] = 'Each bar shows the total number of reads aligned to each amplicon, divided into the reads that are modified and unmodified. The vertical line shows the cutoff for analysis, set using the --min_reads_to_use_region parameter.'
            crispresso2_info['summary_plot_datas'][plot_name] = [('CRISPRessoWGS summary', os.path.basename(samples_quantification_summary_filename))]

        if not args.suppress_report and not args.suppress_plots:
            if (args.place_report_in_output_folder):
                report_name = _jp("CRISPResso2WGS_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_wgs_report_from_folder(report_name, crispresso2_info, OUTPUT_DIRECTORY, _ROOT)
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

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

        info('Analysis Complete!')
        print(CRISPRessoShared.get_crispresso_footer())
        sys.exit(0)

    except Exception as e:
        print_stacktrace_if_debug()
        error('\n\nERROR: %s' % e)
        sys.exit(-1)

if __name__ == '__main__':
    main()
