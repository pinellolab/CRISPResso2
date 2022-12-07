# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import argparse
from copy import deepcopy
import os
import sys
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoReport
import traceback


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


def check_library(library_name):
    try:
        return __import__(library_name)
    except:
        error('You need to install {0} module to use CRISPRessoPooledWGSCompare!'.format(library_name))
        sys.exit(1)


def check_PooledWGS_output_folder(output_folder):
    quantification_summary_file = os.path.join(
        output_folder, 'SAMPLES_QUANTIFICATION_SUMMARY.txt',
    )

    if os.path.exists(quantification_summary_file):
        return quantification_summary_file
    else:
        raise PooledWGSOutputFolderIncompleteException(
            'The folder {0} is not a valid CRISPRessoPooled or CRISPRessoWGS output folder.'.format(
                output_folder,
            ),
        )


pd = check_library('pandas')


###EXCEPTIONS############################
class PooledWGSOutputFolderIncompleteException(Exception):
    pass


_ROOT = os.path.abspath(os.path.dirname(__file__))
CRISPResso_compare_to_call = 'CRISPRessoCompare'


def main():
    try:
        description = [
            '~~~CRISPRessoPooledWGSCompare~~~',
            '-Comparison of two CRISPRessoPooled or CRISPRessoWGS analyses-',
        ]

        compare_header = r'''
 ____________________________________
| __  __  __     __ __        __  __ |
||__)/  \/  \|  |_ |  \ /|  |/ _ (_  |
||   \__/\__/|__|__|__// |/\|\__)__) |
|   __ __      __      __  __        |
|  /  /  \|\/||__) /\ |__)|_         |
|  \__\__/|  ||   /--\| \ |__        |
|____________________________________|
        '''
        compare_header = CRISPRessoShared.get_crispresso_header(
            description, compare_header,
        )
        print(compare_header)

        parser = argparse.ArgumentParser(
            description='CRISPRessoPooledWGSCompare Parameters',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument(
            'crispresso_pooled_wgs_output_folder_1',
            type=str,
            help='First output folder with CRISPRessoPooled or CRISPRessoWGS analysis',
        )
        parser.add_argument(
            'crispresso_pooled_wgs_output_folder_2',
            type=str,
            help='Second output folder with CRISPRessoPooled or CRISPRessoWGS analysis',
        )

        #OPTIONALS
        parser.add_argument('-n', '--name',  help='Output name', default='')
        parser.add_argument(
            '-n1',
            '--sample_1_name',
            help='Sample 1 name',
            default='Sample_1',
        )
        parser.add_argument(
            '-n2',
            '--sample_2_name',
            help='Sample 2 name',
            default='Sample_2',
        )
        parser.add_argument('-o', '--output_folder',  help='', default='')
        parser.add_argument(
            '-p',
            '--n_processes',
            type=str,
            help="""
Specify the number of processes to use for analysis.
Please use with caution since increasing this parameter will significantly
increase the memory required to run CRISPResso. Can be set to 'max'.
            """,
            default='1',
        )
        parser.add_argument(
            '--reported_qvalue_cutoff',
            type=float,
            help='Q-value cutoff for signifance in tests for differential editing. Each base position is tested (for insertions, deletions, substitutions, and all modifications) using Fisher\'s exact test, followed by Bonferonni correction. The number of bases with a significance below this threshold in the quantification window are counted and reported in the output summary.',
            default=0.05
        )
        parser.add_argument(
            '--min_frequency_alleles_around_cut_to_plot',
            type=float,
            help='Minimum %% reads required to report an allele in the alleles table plot.',
            default=0.2,
        )
        parser.add_argument(
            '--max_rows_alleles_around_cut_to_plot',
            type=int,
            help='Maximum number of rows to report in the alleles table plot. ',
            default=50,
        )
        parser.add_argument(
            '--place_report_in_output_folder',
            help='If true, report will be written inside the CRISPResso output folder. By default, the report will be written one directory up from the report output.',
            action='store_true',
        )
        parser.add_argument(
            '--suppress_report',
            help='Suppress output report',
            action='store_true',
        )
        parser.add_argument(
            '--debug',
            help='Show debug messages',
            action='store_true',
        )
        parser.add_argument('--zip_output', help="If set, the output will be placed in a zip folder.", action='store_true')

        args = parser.parse_args()
        debug_flag = args.debug

        crispresso_compare_options = [
            'reported_qvalue_cutoff',
            'min_frequency_alleles_around_cut_to_plot',
            'max_rows_alleles_around_cut_to_plot',
            'place_report_in_output_folder',
            'suppress_report',
            'debug',
        ]

        if args.zip_output and not args.place_report_in_output_folder:
            logger.warn('Invalid arguement combination: If zip_output is True then place_report_in_output_folder must also be True. Setting place_report_in_output_folder to True.')
            args.place_report_in_output_folder = True

        sample_1_name = CRISPRessoShared.slugify(args.sample_1_name)
        sample_2_name = CRISPRessoShared.slugify(args.sample_2_name)

        n_processes = 1
        if args.n_processes == 'max':
            n_processes = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes = int(args.n_processes)

        # check that the CRISPRessoPooled output is present
        quantification_summary_file_1 = check_PooledWGS_output_folder(
            args.crispresso_pooled_wgs_output_folder_1,
        )
        quantification_summary_file_2 = check_PooledWGS_output_folder(
            args.crispresso_pooled_wgs_output_folder_2,
        )

        # create outputfolder and initialize the log
        get_name_from_folder = lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPRessoPooled_on_', '').replace('CRISPRessoWGS_on_', '')

        if not args.name:
            database_id = '{0}_VS_{1}'.format(
                get_name_from_folder(
                    args.crispresso_pooled_wgs_output_folder_1,
                ),
                get_name_from_folder(
                    args.crispresso_pooled_wgs_output_folder_2,
                ),
            )
        else:
            database_id = CRISPRessoShared.slugify(args.name)

        OUTPUT_DIRECTORY = 'CRISPRessoPooledWGSCompare_on_{0}'.format(database_id)

        if args.output_folder:
            OUTPUT_DIRECTORY = os.path.join(
                os.path.abspath(args.output_folder), OUTPUT_DIRECTORY,
            )

        _jp = lambda filename: os.path.join(OUTPUT_DIRECTORY, filename) #handy function to put a file in the output directory
        log_filename = _jp('CRISPRessoPooledWGSCompare_RUNNING_LOG.txt')

        try:
            info('Creating Folder %s' % OUTPUT_DIRECTORY)
            os.makedirs(OUTPUT_DIRECTORY)
            info('Done!')
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename = _jp('CRISPRessoPooledWGSCompare_RUNNING_LOG.txt')
        logger.addHandler(logging.FileHandler(log_filename))
        logger.addHandler(CRISPRessoShared.StatusHandler(_jp('CRISPResso_status.txt')))

        with open(log_filename, 'w+') as outfile:
            outfile.write(
                '[Command used]:\nCRISPRessoPooledWGSCompare {0}\n\n[Execution log]:\n'.format(
                    ' '.join(sys.argv),
                ),
            )

        crispresso2Compare_info_file = os.path.join(OUTPUT_DIRECTORY,'CRISPResso2PooledWGSCompare_info.json')
        crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} #keep track of all information for this run to be saved at the end of the run
        crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        crispresso2_info['running_info']['args'] = deepcopy(args)

        crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        crispresso2_info['results']['general_plots']['summary_plot_names'] = []
        crispresso2_info['results']['general_plots']['summary_plot_titles'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_labels'] = {}
        crispresso2_info['results']['general_plots']['summary_plot_datas'] = {}

        save_png = True
        if args.suppress_report:
            save_png = False

        # load data and calculate the difference
        df_quant_1 = pd.read_csv(quantification_summary_file_1, sep='\t')
        df_quant_2 = pd.read_csv(quantification_summary_file_2, sep='\t')
#        df_comp=df_quant_1.set_index(['Name','Amplicon']).join(df_quant_2.set_index(['Name','Amplicon']),lsuffix='_%s' % args.sample_1_name,rsuffix='_%s' % args.sample_2_name)
        df_comp = df_quant_1.set_index('Name').join(
            df_quant_2.set_index('Name'),
            lsuffix='_{0}'.format(sample_1_name),
            rsuffix='_{0}'.format(sample_2_name),
        )

        print('looking for ' + '({0}-{1})_Unmodified%'.format(sample_1_name, sample_2_name))
        df_comp[
            '({0}-{1})_Unmodified%'.format(sample_1_name, sample_2_name)
        ] = df_comp['Unmodified%_{0}'.format(sample_1_name)] - df_comp[
            'Unmodified%_{0}'.format(sample_2_name)
        ]

        df_comp.fillna('NA').to_csv(_jp('COMPARISON_SAMPLES_QUANTIFICATION_SUMMARIES.txt'), sep='\t')

        # now run CRISPRessoCompare for the pairs for wich we have data in both folders
        crispresso_cmds = []
        processed_regions = []
        processed_region_folder_names = {}
        processed_region_html_files = {}
        for idx, row in df_comp.iterrows():
            if idx in processed_regions:
                continue
            if row.isnull().any():
                warn('Skipping sample {0} since it was not processed in one or both conditions'.format(idx))
            else:
                processed_regions.append(idx)
                crispresso_output_folder_1 = os.path.join(
                    args.crispresso_pooled_wgs_output_folder_1,
                    'CRISPResso_on_{0}'.format(idx),
                )
                crispresso_output_folder_2 = os.path.join(
                    args.crispresso_pooled_wgs_output_folder_2,
                    'CRISPResso_on_{0}'.format(idx),
                )
                compare_output_name = '{0}_{1}_VS_{2}'.format(
                    idx, sample_1_name, sample_2_name,
                )
                crispresso_compare_cmd = CRISPResso_compare_to_call + \
                    ' "{0}" "{1}" -o "{2}" -n {3} -n1 "{4}" -n2 "{5}" '.format(
                      crispresso_output_folder_1,
                      crispresso_output_folder_2,
                      OUTPUT_DIRECTORY,
                      compare_output_name,
                      args.sample_1_name + '_' + idx,
                      args.sample_2_name + '_' + idx,
                    )

                crispresso_compare_cmd = CRISPRessoShared.propagate_crispresso_options(
                    crispresso_compare_cmd, crispresso_compare_options, args,
                )
                info('Running CRISPRessoCompare:%s' % crispresso_compare_cmd)
                crispresso_cmds.append(crispresso_compare_cmd)

                sub_folder = os.path.join(
                    OUTPUT_DIRECTORY,
                    'CRISPRessoCompare_on_' + compare_output_name,
                )
                this_sub_html_file = os.path.basename(sub_folder)+".html"
                if args.place_report_in_output_folder:
                    this_sub_html_file = os.path.join(
                        os.path.basename(sub_folder),
                        "CRISPResso2Compare_report.html",
                    )
                processed_region_html_files[idx] = this_sub_html_file
                processed_region_folder_names[idx] = compare_output_name

        CRISPRessoMultiProcessing.run_crispresso_cmds(
            crispresso_cmds, n_processes, 'Comparison',
        )
        crispresso2_info['results']['processed_regions'] = processed_regions
        crispresso2_info['results']['processed_region_folder_names'] = processed_region_folder_names

        header_string = ''
        sig_count_summary_lines = []
        for region in processed_regions:
            processed_region_folder_name = processed_region_folder_names[region]
            processed_region_folder = os.path.join(OUTPUT_DIRECTORY, 'CRISPRessoCompare_on_'+processed_region_folder_name)
            run_info = CRISPRessoShared.load_crispresso_info(processed_region_folder, 'CRISPResso2Compare_info.json')
            sig_counts_filename = run_info['running_info']['sig_counts_report_location']
            this_sig_filepath = os.path.join(processed_region_folder, sig_counts_filename)
            if os.path.exists(this_sig_filepath):
                with open(this_sig_filepath, 'r') as fin:
                    header_string = fin.readline()
                    for line in fin:
                        sig_count_summary_lines.append(region + "\t" + line)
        sig_count_summary_file = _jp("CRISPRessoPooledWGSCompare_significant_base_count_summary.txt")
        with open(sig_count_summary_file, 'w') as fout:
            fout.write('Sample\t'+header_string)
            fout.writelines(sig_count_summary_lines)

        if not args.suppress_report:
            if args.place_report_in_output_folder:
                report_name = _jp("CRISPResso2PooledWGSCompare_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_multi_report(
                processed_regions,
                processed_region_html_files,
                report_name,
                OUTPUT_DIRECTORY,
                _ROOT,
                'CRISPREssoPooledWGSCompare Report<br>{0} vs {1}'.format(
                    sample_1_name, sample_2_name,
                ),
            )
            crispresso2_info['running_info']['report_location'] = report_name
            crispresso2_info['running_info']['report_filename'] = os.path.basename(report_name)

        CRISPRessoShared.write_crispresso_info(
            crispresso2Compare_info_file, crispresso2_info,
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
