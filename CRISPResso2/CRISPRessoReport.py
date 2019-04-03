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
    fig_datas = {}
#    print('crispresso_report file: ' + crispresso_report_file + ' crispresso_folder : ' + crispresso_folder + ' root: ' + _ROOT)

    def add_fig_if_exists(fig_name,fig_root,fig_title,fig_caption,fig_data,
        amplicon_fig_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions,amplicon_fig_datas):
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
                amplicon_fig_datas[fig_name] = []
                for (data_caption,data_file) in fig_data:
                    if os.path.exists(os.path.join(crispresso_folder,data_file)):
                        amplicon_fig_datas[fig_name].append((data_caption,data_file))

    global_fig_names= []
    for fig in ['1a','1b','1c','5a','6a','8a']:
        fig_name = 'plot_'+ fig
        if fig_name + '_root' in run_data:
            add_fig_if_exists(fig_name,run_data[fig_name + '_root'],'Figure ' + fig,run_data[fig_name + '_caption'],run_data[fig_name+'_data'],
                global_fig_names,fig_locs,fig_titles,fig_captions,fig_datas)


    amplicons = []
    for amplicon_name in run_data['ref_names']:
        amplicons.append(amplicon_name)
        amplicon_fig_names = []
        amplicon_fig_locs = {}
        amplicon_fig_titles = {}
        amplicon_fig_captions = {}
        amplicon_fig_datas = {}

        #fig 2b's
        amplicon_fig_2b_names = []
        if 'plot_2b_roots' in run_data['refs'][amplicon_name]:
            for idx,plot_2b_root in enumerate(run_data['refs'][amplicon_name]['plot_2b_roots']):
                fig_name = "plot_2b_" + str(idx)
                add_fig_if_exists(fig_name,plot_2b_root,'Figure ' + fig_name + ' sgRNA ' + str(idx+1),run_data['refs'][amplicon_name]['plot_2b_captions'][idx],run_data['refs'][amplicon_name]['plot_2b_datas'][idx],
                    amplicon_fig_2b_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions,amplicon_fig_datas)


        for fig in ['2a','3a','3b','4a','4b','4c','4d','4e','4f','5','6','7','8','10a','10b','10c','10d','10e','10f','10g']:
            fig_name = 'plot_'+ fig
            if fig_name + '_root' in run_data['refs'][amplicon_name]:
                add_fig_if_exists(fig_name,run_data['refs'][amplicon_name][fig_name + '_root'],'Figure ' + fig_name,run_data['refs'][amplicon_name][fig_name + '_caption'],run_data['refs'][amplicon_name][fig_name + '_data'],
                        amplicon_fig_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions,amplicon_fig_datas)

        #fig 9's
        amplicon_fig_9_names = []
        if 'plot_9_roots' in run_data['refs'][amplicon_name]:
            for idx,plot_9_root in enumerate(run_data['refs'][amplicon_name]['plot_9_roots']):
                fig_name = "plot_9_" + str(idx)
                add_fig_if_exists(fig_name,plot_9_root,'Figure ' + fig_name + ' sgRNA ' + str(idx+1),run_data['refs'][amplicon_name]['plot_9_captions'][idx],run_data['refs'][amplicon_name]['plot_9_datas'][idx],
                    amplicon_fig_9_names,amplicon_fig_locs,amplicon_fig_titles,amplicon_fig_captions,amplicon_fig_datas)



        fig_names[amplicon_name] = amplicon_fig_names
        fig_2b_names[amplicon_name] = amplicon_fig_2b_names
        fig_9_names[amplicon_name] = amplicon_fig_9_names

        fig_locs[amplicon_name] = amplicon_fig_locs
        fig_titles[amplicon_name] = amplicon_fig_titles
        fig_captions[amplicon_name] = amplicon_fig_captions
        fig_datas[amplicon_name] = amplicon_fig_datas

    report_display_name = ""
    if run_data['args'].name != "":
        report_display_name = run_data['args'].name

    report_data={'amplicons':amplicons,'fig_names':fig_names,'fig_2b_names':fig_2b_names,'fig_9_names':fig_9_names,
            'fig_locs':fig_locs,'fig_titles':fig_titles,'fig_captions':fig_captions,'fig_datas':fig_datas,'run_data':run_data,
            'command_used':run_data['command_used'],'params':run_data['args_string'],'report_display_name':report_display_name}

    j2_env = Environment(loader=FileSystemLoader(os.path.join(_ROOT,'templates')))
    template = j2_env.get_template('report.html')

#    dest_dir = os.path.dirname(crispresso_report_file)
#    shutil.copy2(os.path.join(_ROOT,'templates','CRISPResso_justcup.png'),dest_dir)
#    shutil.copy2(os.path.join(_ROOT,'templates','favicon.ico'),dest_dir)

    outfile = open(crispresso_report_file,'w')
    outfile.write(template.render(report_data=report_data))
    outfile.close()

def make_batch_report_from_folder(crispressoBatch_report_file,crispresso2_info,batch_folder,_ROOT):
    batch_names = crispresso2_info['completed_batch_arr']
    display_names = crispresso2_info['batch_input_names']

    window_nuc_pct_quilts = [x + ".pdf" for x in crispresso2_info['window_nuc_pct_quilt_plot_names']]
    nuc_pct_quilts = [x + ".pdf" for x in crispresso2_info['nuc_pct_quilt_plot_names']]

    window_nuc_conv_plots = [x + ".pdf" for x in crispresso2_info['window_nuc_conv_plot_names']]
    nuc_conv_plots = [x + ".pdf" for x in crispresso2_info['nuc_conv_plot_names']]


    sub_html_files = {}
    run_names = []
    for name in batch_names:
        display_name = display_names[name]
        sub_folder = 'CRISPResso_on_' + name
        info_file = os.path.join(batch_folder,sub_folder,'CRISPResso2_info.pickle')
        if not os.path.exists(info_file):
            raise Exception('CRISPResso run %s is not complete. Cannot add to batch report.'% sub_folder)
        run_data = cp.load(open(info_file,'rb'))
        if not 'report_filename' in run_data:
            raise Exception('CRISPResso run %s has no report. Cannot add to batch report.'% sub_folder)
        sub_html_files[display_name] = os.path.join(sub_folder,os.path.basename(run_data['report_filename']))
        run_names.append(display_name)

    make_multi_report(run_names,sub_html_files,crispressoBatch_report_file,_ROOT,'CRISPResso Batch Output',
        window_nuc_pct_quilts=window_nuc_pct_quilts,
        nuc_pct_quilts=nuc_pct_quilts,
        window_nuc_conv_plots=window_nuc_conv_plots,
        nuc_conv_plots=nuc_conv_plots)


def make_pooled_report_from_folder(crispresso_report_file,crispresso2_info,folder,_ROOT):
    names_arr = crispresso2_info['good_region_names']
    make_multi_report_from_folder(crispresso2_info,names_arr,'CRISPResso Pooled Output',crispresso_report_file,folder,_ROOT)

def make_meta_report_from_folder(crispresso_report_file,crispresso2_info,folder,_ROOT):
    names_arr = crispresso2_info['meta_names_arr']
    input_names = crispresso2_info['meta_input_names']
    make_multi_report_from_folder(crispresso2_info,names_arr,'CRISPResso Meta Output',crispresso_report_file,folder,_ROOT,display_names=input_names)

def make_wgs_report_from_folder(crispresso_report_file,crispresso2_info,folder,_ROOT):
    names_arr = crispresso2_info['good_region_names']
    make_multi_report_from_folder(crispresso2_info,names_arr,'CRISPResso WGS Output',crispresso_report_file,folder,_ROOT)

def make_multi_report_from_folder(crispresso2_info,names_arr,report_name,crispresso_report_file,folder,_ROOT,display_names=None):
    """
    Prepares information to make a report of multiple CRISPResso runs - like CRISPRessoWGS or CRISPRessoPooled

    Parameters:
    crispresso2_info (dict): information from the crispresso multi run
    names_arr (arr of strings): Names of the crispresso runs
    report_name (string): text to be shown at top of report
    crispresso_report_file (string): path to write report to
    folder (string): folder containing crispresso runs
    _ROOT (string): location of crispresso assets (images, templates, etc)
    display_names (dict): report_name->display_name; Titles to be shown for crispresso runs (if different from names_arr, e.g. if display_names have spaces or bad chars, they won't be the same as names_arr)

    Returns:
    Nothin
    """

    summary_plot_names = []
    if 'summary_plot_names' in crispresso2_info:
        summary_plot_names = crispresso2_info['summary_plot_names']
    summary_plot_titles = {}
    if 'summary_plot_titles' in crispresso2_info:
        summary_plot_titles = crispresso2_info['summary_plot_titles']
    summary_plot_labels = {}
    if 'summary_plot_labels' in crispresso2_info:
        summary_plot_labels = crispresso2_info['summary_plot_labels']
    summary_plot_datas = {}
    if 'summary_plot_datas' in crispresso2_info:
        summary_plot_datas = crispresso2_info['summary_plot_datas']

    run_names = []
    sub_html_files = {}
    sub_2a_labels = {}
    sub_2a_pdfs = {}

    for name in names_arr:
        display_name = name
        if display_names is not None:
            display_name = display_names[name]

        folder_name = 'CRISPResso_on_%s' % name
        sub_folder = os.path.join(folder,folder_name)
        info_file = os.path.join(sub_folder,'CRISPResso2_info.pickle')
        if not os.path.exists(info_file):
            raise Exception('CRISPResso run %s is not complete. Cannot add to report.'% sub_folder)
        run_data = cp.load(open(info_file,'rb'))
        if not 'report_filename' in run_data:
            raise Exception('CRISPResso run %s has no report. Cannot add to report.'% sub_folder)

        run_names.append(display_name)
        sub_html_files[display_name] = os.path.join(folder_name,os.path.basename(run_data['report_filename']))

        this_sub_2a_labels = []
        this_sub_2a_pdfs = []
        for ref_name in run_data['ref_names']:
            pdf_file = run_data['refs'][ref_name]['plot_2a_root']+".pdf"
            if os.path.exists(pdf_file):
                this_sub_2a_pdfs.append(run_data['refs'][ref_name]['plot_2a_root']+".pdf")
                this_sub_2a_labels.append("Nucleotide distribution across " + ref_name)

        sub_2a_labels[display_name] = this_sub_2a_labels
        sub_2a_pdfs[display_name] = this_sub_2a_pdfs


    make_multi_report(run_names,sub_html_files,crispresso_report_file,_ROOT,report_name,
            summary_plot_names=summary_plot_names,summary_plot_titles=summary_plot_titles,summary_plot_labels=summary_plot_labels,summary_plot_datas=summary_plot_datas)

def make_multi_report(run_names,sub_html_files,crispresso_multi_report_file,_ROOT,report_name,
    window_nuc_pct_quilts=[],
    nuc_pct_quilts=[],
    window_nuc_conv_plots=[],
    nuc_conv_plots=[],
    summary_plot_names=[],
    summary_plot_titles={},
    summary_plot_labels={},
    summary_plot_datas={}
):
        """
        Makes an HTML report for a run containing multiple crispresso runs

        Parameters:
        run_names (arr of strings): names of runs
        sub_html_files (dict): dict of run_name->file_loc
        crispresso_multi_report_file (string): path of file to write to
        report_name (string): description of report type to be shown at top of report

        summary_plot_names (list): list of plot names - keys for following dicts
        summary_plot_titles (dict): dict of plot_name->plot_title
        summary_plot_labels (dict): dict of plot_name->plot_label
        summary_plot_datas (dict): dict of plot_name->(datafile_description, data_filename)

        """

        def dirname(path):
            return os.path.basename(os.path.dirname(path))
        j2_env = Environment(loader=FileSystemLoader(os.path.join(_ROOT,'templates')))
        j2_env.filters['dirname'] = dirname
        template = j2_env.get_template('multiReport.html')

        outfile = open(crispresso_multi_report_file,'w')
        outfile.write(template.render(window_nuc_pct_quilts=window_nuc_pct_quilts,nuc_pct_quilts=nuc_pct_quilts,
            window_nuc_conv_plots=window_nuc_conv_plots,nuc_conv_plots=nuc_conv_plots,
            summary_plot_names=summary_plot_names,summary_plot_titles=summary_plot_titles,summary_plot_labels=summary_plot_labels,summary_plot_datas=summary_plot_datas,
            run_names=run_names,sub_html_files=sub_html_files,report_name=report_name))
        outfile.close()
