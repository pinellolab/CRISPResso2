'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import os
import sys
from jinja2 import Environment, FileSystemLoader
import shutil

running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

if running_python3:
    import pickle as cp #python 3
else:
    import cPickle as cp #python 2.7

def make_report_from_folder(crispresso_report_file,crispresso_folder,_ROOT):
    info_file = os.path.join(crispresso_folder,'CRISPResso2_info.pickle')
    if not os.path.exists(info_file):
        raise Exception('CRISPResso run is not complete. Cannot create report for run at ' + crispresso_folder)

    run_data = cp.load(open(info_file,'rb'))
    make_report(run_data,crispresso_report_file,crispresso_folder,_ROOT)

def make_report(run_data,crispresso_report_file,crispresso_folder,_ROOT):

    #dicts for each amplicon fig_names[amp_name] = [list of fig names]
    #                        fig_locs[amp_name][fig_name] = figure location
    fig_names = {} #all except for the figure 1 (which is common to all amplicons)
    fig_2b_names = {} #multiple -- one per sgRNA
    fig_9_names = {} # multiple -- one per sgRNA
    fig_locs = {}
    fig_titles = {}
    fig_captions = {}
    data_files = {}
#    print('crispresso_report file: ' + crispresso_report_file + ' crispresso_folder : ' + crispresso_folder + ' root: ' + _ROOT)

    def add_fig_if_exists(fig_name,fig_root,fig_title,fig_caption,
        amplicon_fig_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions):
            """
            Helper function to add figure if the file exists
            if fig at filename exists,
            amplicon_figs[figname] is set to that file
            """
            #fullpath=os.path.join(crispresso_folder,fig_root+'.png')
            fullpath=os.path.join(crispresso_folder,fig_root+'.png')
#            print('adding file ' + fig_root + ' at ' + fullpath)
            if os.path.exists(fullpath):
                amplicon_fig_names.append(fig_name)
                #amplicon_fig_locs[fig_name]=os.path.basename(fig_root+'.png')
                amplicon_fig_locs[fig_name]=os.path.basename(fig_root)
                amplicon_fig_titles[fig_name] = fig_title
                amplicon_fig_captions[fig_name] = fig_caption

    global_fig_names= []
    for fig in ['1a','1b','1c','5a','6a','8a']:
        fig_name = 'plot_'+ fig
        if fig_name + '_root' in run_data:
            add_fig_if_exists(fig_name,run_data[fig_name + '_root'],'Figure ' + fig,run_data[fig_name + '_caption'],
                global_fig_names,fig_locs,fig_titles,fig_captions)


    amplicons = []
    for amplicon_name in run_data['ref_names']:
        amplicons.append(amplicon_name)
        amplicon_fig_names = []
        amplicon_fig_locs = {}
        amplicon_fig_titles = {}
        amplicon_fig_captions = {}

        #fig 2b's
        amplicon_fig_2b_names = []
        if 'plot_2b_roots' in run_data['refs'][amplicon_name]:
            for idx,plot_2b_root in enumerate(run_data['refs'][amplicon_name]['plot_2b_roots']):
                fig_name = "plot_2b_" + str(idx)
                add_fig_if_exists(fig_name,plot_2b_root,'Figure ' + fig_name + ' sgRNA ' + str(idx+1),run_data['refs'][amplicon_name]['plot_2b_captions'][idx],
                    amplicon_fig_2b_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions)


        for fig in ['2a','3a','3b','4a','4b','4c','4d','4e','4f','5','6','7','8','10a','10b','10c','10d','10e','10f','10g']:
            fig_name = 'plot_'+ fig
            if fig_name + '_root' in run_data['refs'][amplicon_name]:
                add_fig_if_exists(fig_name,run_data['refs'][amplicon_name][fig_name + '_root'],'Figure ' + fig_name,run_data['refs'][amplicon_name][fig_name + '_caption'],
                        amplicon_fig_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions)

        #fig 9's
        amplicon_fig_9_names = []
        if 'plot_9_roots' in run_data['refs'][amplicon_name]:
            for idx,plot_9_root in enumerate(run_data['refs'][amplicon_name]['plot_9_roots']):
                fig_name = "plot_9_" + str(idx)
                add_fig_if_exists(fig_name,plot_9_root,'Figure ' + fig_name + ' sgRNA ' + str(idx+1),run_data['refs'][amplicon_name]['plot_9_captions'][idx],
                    amplicon_fig_9_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions)



        fig_names[amplicon_name] = amplicon_fig_names
        fig_2b_names[amplicon_name] = amplicon_fig_2b_names
        fig_9_names[amplicon_name] = amplicon_fig_9_names

        fig_locs[amplicon_name] = amplicon_fig_locs
        fig_titles[amplicon_name] = amplicon_fig_titles
        fig_captions[amplicon_name] = amplicon_fig_captions

    report_display_name = ""
    if run_data['args'].name != "":
        report_display_name = run_data['args'].name

    report_data={'amplicons':amplicons,'fig_names':fig_names,'fig_2b_names':fig_2b_names,'fig_9_names':fig_9_names,
            'fig_locs':fig_locs,'fig_titles':fig_titles,'fig_captions':fig_captions,'run_data':run_data,
            'command_used':run_data['command_used'],'params':run_data['args_string'],'report_display_name':report_display_name}


    j2_env = Environment(loader=FileSystemLoader(os.path.join(_ROOT,'templates')))
    template = j2_env.get_template('report.html')

    dest_dir = os.path.dirname(crispresso_report_file)
    shutil.copy2(os.path.join(_ROOT,'templates','CRISPResso_justcup.png'),dest_dir)
    shutil.copy2(os.path.join(_ROOT,'templates','favicon.ico'),dest_dir)

    outfile = open(crispresso_report_file,'w')
    outfile.write(template.render(report_data=report_data))
    outfile.close()

def make_batch_report_from_folder(crispressoBatch_report_file,batch_folder,_ROOT):
    all_files = os.listdir(batch_folder)
    window_nuc_pct_quilts = [os.path.basename(x) for x in all_files if x.endswith('.Quantification_Window_Nucleotide_Percentage_Quilt.pdf')]
    nuc_pct_quilts = [os.path.basename(x) for x in all_files if x.endswith('.Nucleotide_Percentage_Quilt.pdf')]

    window_nuc_conv_plots = [os.path.basename(x) for x in all_files if x.endswith('.Quantification_Window_Nucleotide_Conversion.pdf')]
    nuc_conv_plots = [os.path.basename(x) for x in all_files if x.endswith('.Nucleotide_Conversion.pdf')]


    sub_folders = [x for x in all_files if x.startswith('CRISPResso_on_')]
    run_names = []
    sub_html_files = {}
    for sub_folder in sub_folders:
        info_file = os.path.join(batch_folder,sub_folder,'CRISPResso2_info.pickle')
        if not os.path.exists(info_file):
            raise Exception('CRISPResso run %s is not complete. Cannot add to batch report.'% sub_folder)
        run_data = cp.load(open(info_file,'rb'))
        if not 'report_filename' in run_data:
            raise Exception('CRISPResso run %s has no report. Cannot add to batch report.'% sub_folder)
        run_name = run_data['args'].name
        run_names.append(run_name)
        sub_html_files[run_name] = os.path.join(sub_folder,os.path.basename(run_data['report_filename']))

    make_multi_report(run_names,sub_html_files,crispressoBatch_report_file,_ROOT,'CRISPResso Batch',
        window_nuc_pct_quilts=window_nuc_pct_quilts,
        nuc_pct_quilts=nuc_pct_quilts,
        window_nuc_conv_plots=window_nuc_conv_plots,
        nuc_conv_plots=nuc_conv_plots)

def make_pooled_report_from_folder(crispressoPooled_report_file,pooled_folder,_ROOT):
    all_files = os.listdir(pooled_folder)

    sub_folders = [x for x in all_files if x.startswith('CRISPResso_on_')]
    run_names = []
    sub_html_files = {}

    sub_2a_labels = {}
    sub_2a_pdfs = {}

    for sub_folder in sub_folders:
        info_file = os.path.join(pooled_folder,sub_folder,'CRISPResso2_info.pickle')
        if not os.path.exists(info_file):
            raise Exception('CRISPResso run %s is not complete. Cannot add to pooled report.'% sub_folder)
        run_data = cp.load(open(info_file,'rb'))
        if not 'report_filename' in run_data:
            raise Exception('CRISPResso run %s has no report. Cannot add to pooled report.'% sub_folder)

        run_name = run_data['args'].name
        run_names.append(run_name)
        sub_html_files[run_name] = os.path.join(sub_folder,os.path.basename(run_data['report_filename']))

        this_sub_2a_labels = []
        this_sub_2a_pdfs = []
        for ref_name in run_data['ref_names']:
            this_sub_2a_labels.append("Nucleotide distribution across " + ref_name)
            this_sub_2a_pdfs.append(run_data['refs'][ref_name]['plot_2a_root']+".pdf")

        sub_2a_labels[run_name] = this_sub_2a_labels
        sub_2a_pdfs[run_name] = this_sub_2a_pdfs

    make_multi_report(run_names,sub_html_files,crispressoPooled_report_file,_ROOT,'CRISPResso Pooled')

def make_multi_report(run_names,sub_html_files,crispresso_multi_report_file,_ROOT,crispresso_mode,
    window_nuc_pct_quilts=[],
    nuc_pct_quilts=[],
    window_nuc_conv_plots=[],
    nuc_conv_plots=[]
):

        def dirname(path):
            return os.path.basename(os.path.dirname(path))
        j2_env = Environment(loader=FileSystemLoader(os.path.join(_ROOT,'templates')))
        j2_env.filters['dirname'] = dirname
        template = j2_env.get_template('multiReport.html')

        dest_dir = os.path.dirname(crispresso_multi_report_file)
        shutil.copy2(os.path.join(_ROOT,'templates','CRISPResso_justcup.png'),dest_dir)
        shutil.copy2(os.path.join(_ROOT,'templates','favicon.ico'),dest_dir)

        outfile = open(crispresso_multi_report_file,'w')
        outfile.write(template.render(window_nuc_pct_quilts=window_nuc_pct_quilts,nuc_pct_quilts=nuc_pct_quilts,
            window_nuc_conv_plots=window_nuc_conv_plots,nuc_conv_plots=nuc_conv_plots,run_names=run_names,sub_html_files=sub_html_files,crispresso_mode=crispresso_mode))
        outfile.close()
