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
import cPickle as cp


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


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoPooledWGSCompare!' % library_name)
                sys.exit(1)


def check_PooledWGS_output_folder(output_folder):
    quantification_summary_file=os.path.join(output_folder,'SAMPLES_QUANTIFICATION_SUMMARY.txt')

    if os.path.exists(quantification_summary_file):
        return quantification_summary_file
    else:
        raise PooledWGSOutputFolderIncompleteException('The folder %s is not a valid CRISPRessoPooled or CRISPRessoWGS output folder.' % output_folder)



pd=check_library('pandas')

###EXCEPTIONS############################
class PooledWGSOutputFolderIncompleteException(Exception):
    pass


_ROOT = os.path.abspath(os.path.dirname(__file__))
#CRISPResso_compare_to_call = os.path.join(os.path.dirname(_ROOT),'CRISPRessoCompare.py')
CRISPResso_compare_to_call ='CRISPRessoCompare'



def main():
    try:
        description = ['~~~CRISPRessoPooledWGSCompare~~~','-Comparison of two CRISPRessoPooled or CRISPRessoWGS analyses-']

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
        compare_header = CRISPRessoShared.get_crispresso_header(description,compare_header)
        print(compare_header)


        parser = argparse.ArgumentParser(description='CRISPRessoPooledWGSCompare Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('crispresso_pooled_wgs_output_folder_1', type=str,  help='First output folder with CRISPRessoPooled or CRISPRessoWGS analysis')
        parser.add_argument('crispresso_pooled_wgs_output_folder_2', type=str,  help='Second output folder with CRISPRessoPooled or CRISPRessoWGS analysis')

        #OPTIONALS
        parser.add_argument('-n','--name',  help='Output name', default='')
        parser.add_argument('-n1','--sample_1_name',  help='Sample 1 name', default='Sample_1')
        parser.add_argument('-n2','--sample_2_name',  help='Sample 2 name', default='Sample_2')
        parser.add_argument('-o','--output_folder',  help='', default='')
        parser.add_argument('-p','--n_processes',type=str, help='Specify the number of processes to use for analysis.\
        Please use with caution since increasing this parameter will significantly increase the memory required to run CRISPResso. Can be set to \'max\'.',default='1')
        parser.add_argument('--min_frequency_alleles_around_cut_to_plot', type=float, help='Minimum %% reads required to report an allele in the alleles table plot.', default=0.2)
        parser.add_argument('--max_rows_alleles_around_cut_to_plot',  type=int, help='Maximum number of rows to report in the alleles table plot. ', default=50)
        parser.add_argument('--place_report_in_output_folder',  help='If true, report will be written inside the CRISPResso output folder. By default, the report will be written one directory up from the report output.', action='store_true')
        parser.add_argument('--suppress_report',  help='Suppress output report', action='store_true')
        parser.add_argument('--debug', help='Show debug messages', action='store_true')

        args = parser.parse_args()
        debug_flag = args.debug

        crispresso_compare_options=['min_frequency_alleles_around_cut_to_plot','max_rows_alleles_around_cut_to_plot','place_report_in_output_folder','suppress_report','debug']

        sample_1_name = CRISPRessoShared.slugify(args.sample_1_name)
        sample_2_name = CRISPRessoShared.slugify(args.sample_2_name)

        n_processes = 1
        if args.n_processes == "max":
            n_processes = CRISPRessoMultiProcessing.get_max_processes()
        else:
            n_processes = int(args.n_processes)

        #check that the CRISPRessoPooled output is present
        quantification_summary_file_1=check_PooledWGS_output_folder(args.crispresso_pooled_wgs_output_folder_1)
        quantification_summary_file_2=check_PooledWGS_output_folder(args.crispresso_pooled_wgs_output_folder_2)

        #create outputfolder and initialize the log
        get_name_from_folder=lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPRessoPooled_on_','').replace('CRISPRessoWGS_on_','')

        if not args.name:
                 database_id='%s_VS_%s' % (get_name_from_folder(args.crispresso_pooled_wgs_output_folder_1),get_name_from_folder(args.crispresso_pooled_wgs_output_folder_2))
        else:
                 database_id=CRISPRessoShared.slugify(args.name)


        OUTPUT_DIRECTORY='CRISPRessoPooledWGSCompare_on_%s' % database_id

        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoPooledWGSCompare_RUNNING_LOG.txt')


        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoPooledWGSCompare_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\nCRISPRessoPooledWGSCompare %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        crispresso2Compare_info_file = os.path.join(OUTPUT_DIRECTORY,'CRISPResso2PooledWGSCompare_info.pickle')
        crispresso2_info = {} #keep track of all information for this run to be pickled and saved at the end of the run
        crispresso2_info['version'] = CRISPRessoShared.__version__
        crispresso2_info['args'] = deepcopy(args)

        crispresso2_info['log_filename'] = os.path.basename(log_filename)

        crispresso2_info['summary_plot_names'] = []
        crispresso2_info['summary_plot_titles'] = {}
        crispresso2_info['summary_plot_labels'] = {}
        crispresso2_info['summary_plot_datas'] = {}

        save_png = True
        if args.suppress_report:
            save_png = False


        #load data and calculate the difference
        df_quant_1=pd.read_csv(quantification_summary_file_1,sep="\t")
        df_quant_2=pd.read_csv(quantification_summary_file_2,sep="\t")
#        df_comp=df_quant_1.set_index(['Name','Amplicon']).join(df_quant_2.set_index(['Name','Amplicon']),lsuffix='_%s' % args.sample_1_name,rsuffix='_%s' % args.sample_2_name)
        df_comp=df_quant_1.set_index('Name').join(df_quant_2.set_index('Name'),lsuffix='_%s' % sample_1_name,rsuffix='_%s' % sample_2_name)

        df_comp['(%s-%s)_Unmodified%%' % (sample_1_name,args.sample_2_name)]=df_comp['Unmodified%%_%s' % sample_1_name]-df_comp['Unmodified%%_%s' % args.sample_2_name]

        df_comp.fillna('NA').to_csv(_jp('COMPARISON_SAMPLES_QUANTIFICATION_SUMMARIES.txt'),sep='\t')


        #now run CRISPRessoCompare for the pairs for wich we have data in both folders
        crispresso_cmds = []
        processed_regions = []
        processed_region_folder_names = {}
        processed_region_html_files = {}
        for idx,row in df_comp.iterrows():
            if idx in processed_regions:
                continue
            if row.isnull().any():
                warn('Skipping sample %s since it was not processed in one or both conditions' % idx)
            else:
                processed_regions.append(idx)
                #crispresso_output_folder_1=os.path.join(args.crispresso_pooled_wgs_output_folder_1,'CRISPResso_on_%s' % idx)
                #crispresso_output_folder_2=os.path.join(args.crispresso_pooled_wgs_output_folder_2,'CRISPResso_on_%s' % idx)
                crispresso_output_folder_1=os.path.join(args.crispresso_pooled_wgs_output_folder_1,'CRISPResso_on_%s' % idx)
                crispresso_output_folder_2=os.path.join(args.crispresso_pooled_wgs_output_folder_2,'CRISPResso_on_%s' % idx)
                compare_output_name = "%s_%s_VS_%s"%(idx,sample_1_name,sample_2_name)
                crispresso_compare_cmd=CRISPResso_compare_to_call +' "%s" "%s" -o "%s" -n %s -n1 "%s" -n2 "%s" ' % (crispresso_output_folder_1,
                                                                   crispresso_output_folder_2,
                                                                   OUTPUT_DIRECTORY,
                                                                   compare_output_name,
                                                                   sample_1_name+'_'+idx,
                                                                   sample_2_name+'_'+idx
                                                                  )
                crispresso_compare_cmd=CRISPRessoShared.propagate_crispresso_options(crispresso_compare_cmd,crispresso_compare_options,args)
                info('Running CRISPRessoCompare:%s' % crispresso_compare_cmd)
                crispresso_cmds.append(crispresso_compare_cmd)
                sub_folder = os.path.join(OUTPUT_DIRECTORY,"CRISPRessoCompare_on_"+compare_output_name)
                this_sub_html_file = os.path.basename(sub_folder)+".html"
                if args.place_report_in_output_folder:
                    this_sub_html_file = os.path.join(os.path.basename(sub_folder),"CRISPResso2Compare_report.html")
                processed_region_html_files[idx] = this_sub_html_file
                processed_region_folder_names[idx] = compare_output_name


        CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_cmds,n_processes,'Comparison')
        crispresso2_info['processed_regions'] = processed_regions
        crispresso2_info['processed_region_folder_names'] = processed_region_folder_names

        if not args.suppress_report:
            if args.place_report_in_output_folder:
                report_name = _jp("CRISPResso2PooledWGSCompare_report.html")
            else:
                report_name = OUTPUT_DIRECTORY+'.html'
            CRISPRessoReport.make_multi_report(processed_regions,processed_region_html_files,report_name,OUTPUT_DIRECTORY,_ROOT,'CRISPREssoPooledWGSCompare Report<br>'+sample_1_name+" vs" + sample_2_name)
            crispresso2_info['report_location'] = report_name
            crispresso2_info['report_filename'] = os.path.basename(report_name)

        with open(crispresso2Compare_info_file,"wb") as info_file:
            cp.dump(crispresso2_info, info_file)

        info('All Done!')
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
