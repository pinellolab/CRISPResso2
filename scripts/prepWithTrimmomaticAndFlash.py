import argparse
from concurrent.futures import process
import os
import numpy as np
import pandas as pd
import logging
import sys
import subprocess as sb
from CRISPResso2 import CRISPRessoPlot, CRISPRessoShared, CRISPRessoCORE

def main():
    logging.basicConfig(
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    parser = CRISPRessoShared.getCRISPRessoArgParser(
        parserTitle = "Prepare reads for CRISPResso run using legacy trimming and merging software (i.e., Trimmomatic and Flash)",
        requiredParams=['fastq_r1','name'])

    try: #some of these args may exist in the parser already, so catch those exceptions
        parser.add_argument('--split_interleaved_input', '--split_paired_end',
                            help='Splits a single fastq file containing paired end reads into two files before running CRISPResso',
                            action='store_true')
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--trim_sequences', help='Enable the trimming of Illumina adapters with Trimmomatic',
                            action='store_true')
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--trimmomatic_command', type=str, help='Command to run trimmomatic', default='trimmomatic')
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--trimmomatic_options_string', type=str,
                            help='Override options for Trimmomatic, e.g. "ILLUMINACLIP:/data/NexteraPE-PE.fa:0:90:10:0:true"',
                            default='')
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--flash_command', type=str, help='Command to run flash', default='flash')
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--min_paired_end_reads_overlap', type=int,
                            help='Parameter for the FLASH read merging step. Minimum required overlap length between two reads to provide a confident overlap. ',
                            default=10)
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--max_paired_end_reads_overlap', type=int,
                            help='Parameter for the FLASH merging step.  Maximum overlap length expected in approximately 90%% of read pairs. Please see the FLASH manual for more information.',
                            default=100)
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--stringent_flash_merging',
                            help='Use stringent parameters for flash merging. In the case where flash could merge R1 and R2 reads ambiguously, the expected overlap is calculated as 2*average_read_length - amplicon_length. The flash parameters for --min-overlap and --max-overlap will be set to prefer merged reads with length within 10bp of the expected overlap. These values override the --min_paired_end_reads_overlap or --max_paired_end_reads_overlap CRISPResso parameters.',
                            action='store_true')
    except argparse.ArgumentError as e:
        pass

    try:
        parser.add_argument('--force_merge_pairs', help=argparse.SUPPRESS,
                            action='store_true')  # help=Force-merges R1 and R2 if they cannot be merged using flash (use with caution -- may create non-biological apparent indels at the joining site)
    except argparse.ArgumentError as e:
        pass

    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)


    if args.bam_input:
        raise CRISPRessoShared.BadParameterException('This script is not compatible with bam input')

    if args.stringent_flash_merging and args.amplicon_seq == None:
        raise CRISPRessoShared.BadParameterException('Amplicon sequences must be provided for stringent flash merging')

    clean_name=CRISPRessoShared.slugify(args.name)
    if args.name!= clean_name:
        logger.warning('The specified name %s contained invalid characters and was changed to: %s' % (args.name, clean_name))

    OUTPUT_DIRECTORY=clean_name

    if args.output_folder:
        OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)
    OUTPUT_DIRECTORY=clean_name

    if args.output_folder:
        OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder), OUTPUT_DIRECTORY)

    clean_file_prefix = ""
    if args.file_prefix != "":
        clean_file_prefix = CRISPRessoShared.slugify(args.file_prefix)
        if not clean_file_prefix.endswith("."):
            clean_file_prefix += "."

    _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY, clean_file_prefix + filename) 

    files_to_remove = [] #these files will be deleted at the end of the run

    try:
        os.makedirs(OUTPUT_DIRECTORY)
        logger.info('Creating Folder %s' % OUTPUT_DIRECTORY)
#            info('Done!') #crispresso2 doesn't announce that the folder is created... save some electricity here
    except:
        logger.warning('Folder %s already exists.' % OUTPUT_DIRECTORY)

    log_filename=_jp('log.txt')

    N_READS_INPUT = CRISPRessoCORE.get_n_reads_fastq(args.fastq_r1)

    if N_READS_INPUT == 0:
        raise CRISPRessoShared.BadParameterException('The input contains 0 reads.')

    if args.split_interleaved_input:
        if args.fastq_r2!='':
            raise CRISPRessoShared.BadParameterException('The option --split_interleaved_input is available only when a single fastq file is specified!')
        else:
            logger.info('Splitting paired end single fastq file into two files...')
            args.fastq_r1, args.fastq_r2=CRISPRessoCORE.split_interleaved_fastq(args.fastq_r1,
                output_filename_r1=_jp(os.path.basename(args.fastq_r1.replace('.fastq', '')).replace('.gz', '')+'_splitted_r1.fastq.gz'),
                output_filename_r2=_jp(os.path.basename(args.fastq_r1.replace('.fastq', '')).replace('.gz', '')+'_splitted_r2.fastq.gz'),)
            files_to_remove+=[args.fastq_r1, args.fastq_r2]
            N_READS_INPUT = N_READS_INPUT/2

            logger.info('Done!')

    #Trim and merge reads
    if args.bam_input != '' and args.trim_sequences:
        raise CRISPRessoShared.BadParameterException('Read trimming options are not available with bam input')
    elif args.fastq_r1 != '' and args.fastq_r2 == '': #single end reads
        if not args.trim_sequences: #no trimming or merging required
            output_forward_filename=args.fastq_r1
        else:
            logger.info('Trimming single-end sequences with Trimmomatic...')
            CRISPRessoCORE.check_program(args.trimmomatic_command)
            output_forward_filename=_jp('reads.trimmed.fq.gz')
            #Trimming with trimmomatic
            cmd='%s SE -phred33 %s  %s %s >>%s 2>&1'\
            % (args.trimmomatic_command, args.fastq_r1,
                output_forward_filename,
                args.trimmomatic_options_string.replace('NexteraPE-PE.fa', 'TruSeq3-SE.fa'),
                log_filename)
            #print cmd
            trimmomatic_sb =sb.run(cmd, shell=True)

            if trimmomatic_sb.returncode:
                logger.error(trimmomatic_sb.stdout)
                logger.error(trimmomatic_sb.stderr)
                raise CRISPRessoShared.TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

            files_to_remove += [output_forward_filename]

        processed_output_filename=output_forward_filename

    elif args.fastq_r1 != '' and args.fastq_r2 != '':#paired end reads
        if not args.trim_sequences:
            output_forward_paired_filename=args.fastq_r1
            output_reverse_paired_filename=args.fastq_r2
        else:
            logger.info('Trimming paired-end sequences with Trimmomatic...')
            CRISPRessoCORE.check_program(args.trimmomatic_command)
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
            trimmomatic_sb =sb.run(cmd, shell=True)
            if trimmomatic_sb.returncode:
                raise CRISPRessoShared.TrimmomaticException('TRIMMOMATIC failed to run, please check the log file.')

            files_to_remove += [output_forward_paired_filename]
            files_to_remove += [output_reverse_paired_filename]

            logger.info('Done!')

        #for paired-end reads, merge them
        logger.info('Estimating average read length...')
        if args.debug:
            logger.info('Checking average read length from ' + output_forward_paired_filename)
        if CRISPRessoCORE.get_n_reads_fastq(output_forward_paired_filename):
            avg_read_length=CRISPRessoCORE.get_avg_read_length_fastq(output_forward_paired_filename)
            if args.debug:
                logger.info('Average read length is ' + str(avg_read_length) + ' from ' + output_forward_paired_filename)
        else:
            raise CRISPRessoShared.NoReadsAfterQualityFilteringException('No reads survived the average or single bp quality filtering.')

        #Merging with Flash
        logger.info('Merging paired sequences with Flash...')
        CRISPRessoCORE.check_program(args.flash_command)

        max_amplicon_len = 0 #for flash
        min_amplicon_len = 99**99 #for flash

        amplicon_seq_arr = []
        if args.amplicon_seq is not None:
            amplicon_seq_arr = args.amplicon_seq.split(",")
        if args.expected_hdr_amplicon_seq != "":
            amplicon_seq_arr.append(args.expected_hdr_amplicon_seq)

        for idx, seq in enumerate(amplicon_seq_arr):
            this_seq = seq.strip().upper()
            this_seq_length = len(this_seq)
            if this_seq_length > max_amplicon_len:
                max_amplicon_len = this_seq_length
            if this_seq_length < min_amplicon_len:
                min_amplicon_len = this_seq_length

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
                logger.info('Warning: Reads are longer than amplicon.')
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

        logger.info('Running FLASH command: ' + cmd)
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
                    logger.info('Wrote force-merged reads to ' + new_merged_filename)

        logger.info('Done!')

    #count reads
    N_READS_AFTER_PREPROCESSING=CRISPRessoCORE.get_n_reads_fastq(processed_output_filename)
    
    logger.debug(f'Initial reads: {N_READS_INPUT}')
    logger.debug(f'Reads after preprocessing: {N_READS_AFTER_PREPROCESSING}')
    logger.info('Finished read preparation. Please pass the following file as input to CRISPResso: '+processed_output_filename)


class FlashException(Exception):
    pass
class TrimmomaticException(Exception):
    pass

if __name__ == "__main__":
    main()

