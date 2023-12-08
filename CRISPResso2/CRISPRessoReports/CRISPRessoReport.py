'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''

import os
from jinja2 import Environment, FileSystemLoader
from jinja_partials import generate_render_partial, render_partial
from CRISPResso2 import CRISPRessoShared


def render_template(template_name, jinja2_env, **data):
    """Render a template with partials.

    Parameters
    ----------
    template_name: str
        The name of the template to render. For example, if you have a template
        file called `templates/my_template.html` you would pass in
        `my_template.html`.
    jinja2_env: jinja2.Environment
        The Jinja2 environment being used.
    **data: keyword arguments of any type
        Additional keyword arguments that are passed to the template.

    Returns
    -------
    The rendered template.
    """
    def custom_partial_render(partial_template_name, **partial_data):
        template = jinja2_env.get_template(partial_template_name)
        partial_data.update(
            render_partial=generate_render_partial(
                custom_partial_render,
            ),
        )
        return template.render(**partial_data)
    return render_partial(
        template_name, custom_partial_render, **data,
    )


def make_report_from_folder(crispresso_report_file, crispresso_folder, _ROOT):
    """
    Makes an html report for a crispresso run

    Parameters:
    crispresso_report_file (string): name of the html file to create
    crispresso_folder (string): path to the crispresso output
    _ROOT (string): path to crispresso executables (for templates)

    Returns:
    Nothin
    """
    run_data = CRISPRessoShared.load_crispresso_info(crispresso_folder)
    make_report(run_data, crispresso_report_file, crispresso_folder, _ROOT)


def add_fig_if_exists(fig_name, fig_root, fig_title, fig_caption, fig_data,
                      amplicon_fig_names, amplicon_figures, crispresso_folder):
    """
        Helper function to add figure if the file exists
        if fig at filename exists,
        amplicon_figs[figname] is set to that file
        """
    # fullpath=os.path.join(crispresso_folder,fig_root+'.png')
    pngfullpath = os.path.join(crispresso_folder, fig_root + '.png')
    htmlfullpath = os.path.join(crispresso_folder, fig_root + '.html')
    #            print('adding file ' + fig_root + ' at ' + fullpath)
    if os.path.exists(pngfullpath) or os.path.exists(htmlfullpath):
        amplicon_fig_names.append(fig_name)
        # amplicon_fig_locs[fig_name]=os.path.basename(fig_root+'.png')
        amplicon_figures['locs'][fig_name] = os.path.basename(fig_root)
        amplicon_figures['titles'][fig_name] = fig_title
        amplicon_figures['captions'][fig_name] = fig_caption
        amplicon_figures['datas'][fig_name] = []
        for (data_caption, data_file) in fig_data:
            if os.path.exists(os.path.join(crispresso_folder, data_file)):
                amplicon_figures['datas'][fig_name].append((data_caption, data_file))
        if os.path.exists(htmlfullpath):
            with open(htmlfullpath, encoding="utf-8") as html:
                html_string = "<div align='center'>"
                html_string += html.read()
                html_string += "</div>"
            amplicon_figures['htmls'][fig_name] = html_string


def assemble_figs(run_data, crispresso_folder):
    """
        Helper function create the data structre for the figures
    """
    figures = {'names': {}, 'locs': {}, 'titles': {}, 'captions': {}, 'datas': {}, 'htmls': {}, 'sgRNA_based_names': {}}

    global_fig_names = []
    for fig in ['1a', '1b', '1c', '1d', '5a', '6a', '8a', '11c']:
        fig_name = 'plot_' + fig
        if fig_name + '_root' in run_data['results']['general_plots']:
            add_fig_if_exists(fig_name, run_data['results']['general_plots'][fig_name + '_root'], 'Figure ' + fig,
                              run_data['results']['general_plots'][fig_name + '_caption'],
                              run_data['results']['general_plots'][fig_name + '_data'],
                              global_fig_names, figures, crispresso_folder)

    amplicons = []
    for amplicon_name in run_data['results']['ref_names']:
        amplicons.append(amplicon_name)
        amplicon_figures = {'names': [], 'locs': {}, 'titles': {}, 'captions': {}, 'datas': {}, 'htmls': {}}

        for fig in ['2a', '3a', '3b', '4a', '4b', '4c', '4d', '4e', '4f', '4g', '5', '6', '7', '8', '10a', '10b', '10c',
                    '11a']:
            fig_name = 'plot_' + fig
            if fig_name + '_root' in run_data['results']['refs'][amplicon_name]:
                add_fig_if_exists(fig_name, run_data['results']['refs'][amplicon_name][fig_name + '_root'],
                                  'Figure ' + fig_name,
                                  run_data['results']['refs'][amplicon_name][fig_name + '_caption'],
                                  run_data['results']['refs'][amplicon_name][fig_name + '_data'],
                                  global_fig_names, amplicon_figures, crispresso_folder)

        this_sgRNA_based_fig_names = {}
        for fig in ['2b', '9', '10d', '10e', '10f', '10g', '11b']:
            # fig 2b's
            this_fig_names = []
            if 'plot_' + fig + '_roots' in run_data['results']['refs'][amplicon_name]:
                for idx, plot_root in enumerate(run_data['results']['refs'][amplicon_name]['plot_' + fig + '_roots']):
                    fig_name = "plot_" + fig + "_" + str(idx)
                    add_fig_if_exists(fig_name, plot_root, 'Figure ' + fig_name + ' sgRNA ' + str(idx + 1),
                                      run_data['results']['refs'][amplicon_name]['plot_' + fig + '_captions'][idx],
                                      run_data['results']['refs'][amplicon_name]['plot_' + fig + '_datas'][idx],
                                      this_fig_names, amplicon_figures, crispresso_folder)
            this_sgRNA_based_fig_names[fig] = this_fig_names

        figures['names'][amplicon_name] = amplicon_figures['names']
        figures['sgRNA_based_names'][amplicon_name] = this_sgRNA_based_fig_names

        figures['locs'][amplicon_name] = amplicon_figures['locs']
        figures['titles'][amplicon_name] = amplicon_figures['titles']
        figures['captions'][amplicon_name] = amplicon_figures['captions']
        figures['datas'][amplicon_name] = amplicon_figures['datas']
        figures['htmls'][amplicon_name] = amplicon_figures['htmls']
    data = {'amplicons': amplicons, 'figures': figures}
    return data


def make_report(run_data, crispresso_report_file, crispresso_folder, _ROOT):
    # dicts for each amplicon fig_names[amp_name] = [list of fig names]
    #                        fig_locs[amp_name][fig_name] = figure location
    #    print('crispresso_report file: ' + crispresso_report_file + ' crispresso_folder : ' + crispresso_folder + ' root: ' + _ROOT)
    data = assemble_figs(run_data, crispresso_folder)

    report_display_name = ""
    if run_data['running_info']['args'].name != "":
        report_display_name = run_data['running_info']['args'].name

    # find path between the report and the data (if the report is in another directory vs in the same directory as the data)
    crispresso_data_path = os.path.relpath(crispresso_folder, os.path.dirname(crispresso_report_file))
    if crispresso_data_path == ".":
        crispresso_data_path = ""
    else:
        crispresso_data_path += "/"

    report_data = {
        'amplicons': data['amplicons'],
        'figures': data['figures'],
        'run_data': run_data,
        'report_display_name': report_display_name,
        'crispresso_data_path': crispresso_data_path,
    }

    j2_env = Environment(loader=FileSystemLoader(os.path.join(_ROOT, 'CRISPRessoReports', 'templates')))

    #    dest_dir = os.path.dirname(crispresso_report_file)
    #    shutil.copy2(os.path.join(_ROOT,'templates','CRISPResso_justcup.png'),dest_dir)
    #    shutil.copy2(os.path.join(_ROOT,'templates','favicon.ico'),dest_dir)

    with open(crispresso_report_file, 'w', encoding="utf-8") as outfile:
        outfile.write(render_template(
            'report.html', j2_env, report_data=report_data,
        ))


def make_batch_report_from_folder(crispressoBatch_report_file, crispresso2_info, batch_folder, _ROOT):
    batch_names = crispresso2_info['results']['completed_batch_arr']
    failed_runs = crispresso2_info['results']['failed_batch_arr']
    display_names = crispresso2_info['results']['batch_input_names']

    window_nuc_pct_quilts = crispresso2_info['results']['general_plots']['window_nuc_pct_quilt_plot_names']
    nuc_pct_quilts = crispresso2_info['results']['general_plots']['nuc_pct_quilt_plot_names']

    window_nuc_conv_plots = crispresso2_info['results']['general_plots']['window_nuc_conv_plot_names']
    nuc_conv_plots = crispresso2_info['results']['general_plots']['nuc_conv_plot_names']

    summary_plot_names = []
    if 'summary_plot_names' in crispresso2_info['results']['general_plots']:
        summary_plot_names = crispresso2_info['results']['general_plots']['summary_plot_names']
    summary_plot_titles = {}
    if 'summary_plot_titles' in crispresso2_info['results']['general_plots']:
        summary_plot_titles = crispresso2_info['results']['general_plots']['summary_plot_titles']
    summary_plot_labels = {}
    if 'summary_plot_labels' in crispresso2_info['results']['general_plots']:
        summary_plot_labels = crispresso2_info['results']['general_plots']['summary_plot_labels']
    summary_plot_datas = {}
    if 'summary_plot_datas' in crispresso2_info['results']['general_plots']:
        summary_plot_datas = crispresso2_info['results']['general_plots']['summary_plot_datas']

    allele_modification_heatmap_plot = {}
    if 'allele_modification_heatmap_plot_names' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['names'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_names']
    else:
        allele_modification_heatmap_plot['names'] = []
    if 'allele_modification_heatmap_plot_paths' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['paths'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_paths']
    else:
        allele_modification_heatmap_plot['paths'] = {}
    if 'allele_modification_heatmap_plot_titles' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['titles'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_titles']
    else:
        allele_modification_heatmap_plot['titles'] = []
    if 'allele_modification_heatmap_plot_labels' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['labels'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_labels']
    else:
        allele_modification_heatmap_plot['labels'] = {}
    if 'allele_modification_heatmap_plot_datas' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['datas'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_datas']
    else:
        allele_modification_heatmap_plot['datas'] = {}

    allele_modification_line_plot = {}
    if 'allele_modification_line_plot_names' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['names'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_names']
    else:
        allele_modification_line_plot['names'] = []
    if 'allele_modification_line_plot_paths' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['paths'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_paths']
    else:
        allele_modification_line_plot['paths'] = {}
    if 'allele_modification_line_plot_titles' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['titles'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_titles']
    else:
        allele_modification_line_plot['titles'] = []
    if 'allele_modification_line_plot_labels' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['labels'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_labels']
    else:
        allele_modification_line_plot['labels'] = {}
    if 'allele_modification_line_plot_datas' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['datas'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_datas']
    else:
        allele_modification_line_plot['datas'] = {}

    allele_modification_heatmap_plot['htmls'] = {}
    for heatmap_plot_name, heatmap_plot_path in allele_modification_heatmap_plot['paths'].items():
        with open(heatmap_plot_path, encoding="utf-8") as fh:
            allele_modification_heatmap_plot['htmls'][heatmap_plot_name] = fh.read()

    allele_modification_line_plot['htmls'] = {}
    for line_plot_name, line_plot_path in allele_modification_line_plot['paths'].items():
        with open(line_plot_path, encoding="utf-8") as fh:
            allele_modification_line_plot['htmls'][line_plot_name] = fh.read()

    #find path between the report and the data (if the report is in another directory vs in the same directory as the data)
    crispresso_data_path = os.path.relpath(batch_folder, os.path.dirname(crispressoBatch_report_file))
    if crispresso_data_path == ".":
        crispresso_data_path = ""
    else:
        crispresso_data_path += "/"

    sub_html_files = {}
    run_names = []
    for name in batch_names:
        display_name = display_names[name]
        sub_folder = 'CRISPResso_on_' + name
        crispresso_folder = os.path.join(batch_folder, sub_folder)
        run_data = CRISPRessoShared.load_crispresso_info(crispresso_folder)
        if 'running_info' not in run_data:
            raise Exception(f'CRISPResso run {sub_folder} has no report. Cannot add to batch report.')

        this_sub_html_file = sub_folder + ".html"
        if run_data['running_info']['args'].place_report_in_output_folder:
            this_sub_html_file = os.path.join(sub_folder, run_data['running_info']['report_filename'])
        sub_html_files[display_name] = this_sub_html_file

        run_names.append(display_name)

    output_title = 'CRISPResso Batch Output'
    if crispresso2_info['running_info']['args'].name != '':
        output_title += f"<br/>{crispresso2_info['running_info']['args'].name}"

    make_multi_report(
        run_names,
        failed_runs,
        sub_html_files,
        crispressoBatch_report_file,
        batch_folder,
        _ROOT,
        output_title,
        'batch',
        summary_plots={
            'names': summary_plot_names,
            'titles': summary_plot_titles,
            'labels': summary_plot_labels,
            'datas': summary_plot_datas,
        },
        window_nuc_pct_quilts=window_nuc_pct_quilts,
        nuc_pct_quilts=nuc_pct_quilts,
        window_nuc_conv_plots=window_nuc_conv_plots,
        nuc_conv_plots=nuc_conv_plots,
        allele_modification_heatmap_plot=allele_modification_heatmap_plot,
        allele_modification_line_plot=allele_modification_line_plot,
    )


def make_pooled_report_from_folder(crispresso_report_file, crispresso2_info, folder, _ROOT):
    print(crispresso2_info['results']['failed_batch_arr'])
    names_arr = crispresso2_info['results']['good_region_names']
    output_title = 'CRISPResso Pooled Output'
    if crispresso2_info['running_info']['args'].name != '':
        output_title += f"<br/>{crispresso2_info['running_info']['args'].name}"
    make_multi_report_from_folder(crispresso2_info, names_arr, output_title, crispresso_report_file, folder, _ROOT, 'pooled')


def make_compare_report_from_folder(crispresso_report_file, crispresso2_info, folder, _ROOT):
    names_arr = []
    output_title = 'CRISPResso Compare Output'
    if crispresso2_info['running_info']['args'].name != '':
        output_title += "<br/>{crispresso2_info['running_info']['args'].name}"
    make_multi_report_from_folder(crispresso2_info, names_arr, output_title, crispresso_report_file, folder, _ROOT, 'compare')


def make_meta_report_from_folder(crispresso_report_file, crispresso2_info, folder, _ROOT):
    names_arr = crispresso2_info['meta_names_arr']
    input_names = crispresso2_info['meta_input_names']
    output_title = 'CRISPresso Meta Output'
    if crispresso2_info['running_info']['args'].name != '':
        output_title += "<br/>{crispresso2_info['running_info']['args'].name}"
    make_multi_report_from_folder(crispresso2_info, names_arr, output_title, crispresso_report_file, folder, _ROOT, 'meta',
                                  display_names=input_names)


def make_wgs_report_from_folder(crispresso_report_file, crispresso2_info, folder, _ROOT):
    names_arr = crispresso2_info['results']['good_region_names']
    output_title = 'CRISPResso WGS Output'
    if crispresso2_info['running_info']['args'].name != '':
        output_title += "<br/>{crispresso2_info['running_info']['args'].name}"
    make_multi_report_from_folder(crispresso2_info, names_arr, output_title, crispresso_report_file, folder, _ROOT, 'wgs')


def make_multi_report_from_folder(crispresso2_info, names_arr, report_name, crispresso_report_file, folder, _ROOT, crispresso_tool,
                                  display_names=None):
    """
    Prepares information to make a report of multiple CRISPResso runs - like CRISPRessoWGS or CRISPRessoPooled

    Parameters:
    crispresso2_info (dict): information from the crispresso multi run
    names_arr (arr of strings): Names of the crispresso runs
    report_name (string): text to be shown at top of report
    crispresso_report_file (string): path to write report to
    folder (string): folder containing crispresso runs
    _ROOT (string): location of crispresso assets (images, templates, etc)
    display_names (dict): report_name->display_name; Titles to be shown for crispresso runs 
        (if different from names_arr, e.g. if display_names have spaces or bad chars, they won't be the same as names_arr)

    Returns:
    Nothing
    """

    summary_plot_names = []
    if 'summary_plot_names' in crispresso2_info['results']['general_plots']:
        summary_plot_names = crispresso2_info['results']['general_plots']['summary_plot_names']
    summary_plot_titles = {}
    if 'summary_plot_titles' in crispresso2_info['results']['general_plots']:
        summary_plot_titles = crispresso2_info['results']['general_plots']['summary_plot_titles']
    summary_plot_labels = {}
    if 'summary_plot_labels' in crispresso2_info['results']['general_plots']:
        summary_plot_labels = crispresso2_info['results']['general_plots']['summary_plot_labels']
    summary_plot_datas = {}
    if 'summary_plot_datas' in crispresso2_info['results']['general_plots']:
        summary_plot_datas = crispresso2_info['results']['general_plots']['summary_plot_datas']

    run_names = []
    failed_runs = crispresso2_info['results']['failed_batch_arr']
    print("Failed runs:")
    print(failed_runs)
    sub_html_files = {}
    sub_2a_labels = {}
    sub_2a_pdfs = {}

    for name in names_arr:
        display_name = name
        if display_names is not None:
            display_name = display_names[name]

        folder_name = f'CRISPResso_on_{name}'
        sub_folder = os.path.join(folder, folder_name)
        run_data = CRISPRessoShared.load_crispresso_info(sub_folder)
        if 'running_info' not in run_data:
            raise Exception(f'CRISPResso run {sub_folder} has no report. Cannot add to report.')

        run_names.append(display_name)

        this_sub_html_file = os.path.basename(folder_name) + ".html"
        if run_data['running_info']['args'].place_report_in_output_folder:
            this_sub_html_file = os.path.join(os.path.basename(sub_folder), run_data['running_info']['report_filename'])
        sub_html_files[display_name] = this_sub_html_file

        this_sub_2a_labels = []
        this_sub_2a_pdfs = []
        for ref_name in run_data['results']['ref_names']:
            if 'plot_2a_root' in run_data['results']['refs'][ref_name]:
                pdf_file = run_data['results']['refs'][ref_name]['plot_2a_root'] + ".pdf"
                if os.path.exists(pdf_file):
                    this_sub_2a_pdfs.append(run_data['results']['refs'][ref_name]['plot_2a_root'] + ".pdf")
                    this_sub_2a_labels.append("Nucleotide distribution across " + ref_name)

        sub_2a_labels[display_name] = this_sub_2a_labels
        sub_2a_pdfs[display_name] = this_sub_2a_pdfs
    make_multi_report(
        run_names,
        failed_runs,
        sub_html_files,
        crispresso_report_file,
        folder,
        _ROOT,
        report_name,
        crispresso_tool,
        summary_plots={
            'names': summary_plot_names,
            'titles': summary_plot_titles,
            'labels': summary_plot_labels,
            'datas': summary_plot_datas,
        },
        )


def make_multi_report(
    run_names,
    failed_runs,
    sub_html_files,
    crispresso_multi_report_file,
    crispresso_folder,
    _ROOT,
    report_name,
    crispresso_tool,
    window_nuc_pct_quilts=None,
    nuc_pct_quilts=None,
    window_nuc_conv_plots=None,
    nuc_conv_plots=None,
    summary_plots=None,
    compact_plots_to_show=None,
    allele_modification_heatmap_plot=None,
    allele_modification_line_plot=None,
):
    """
    Makes an HTML report for a run containing multiple crispresso runs

    Parameters:
    run_names (arr of strings): names of runs
    sub_html_files (dict): dict of run_name->file_loc
    crispresso_multi_report_file (string): path of file to write to
    report_name (string): description of report type to be shown at top of report
    crispresso_folder (string): absolute path to the crispresso output
    _ROOT (string): absolute path to the crispresso executable
    summary_plots (dict): a dict with the following keys:
        names (list): list of plot names - keys for following dicts
        titles (dict): dict of plot_name->plot_title
        labels (dict): dict of plot_name->plot_label
        datas (dict): dict of plot_name->[(datafile_description, data_filename), ...]
    compact_plots_to_show (dict): name=>{'href': path to target(report) when user clicks on image, 'img': path to png image to show}
    allele_modification_heatmap_plot (dict): a dict with the following keys:
        names (list): list of plot names for heatmaps, keys for dicts below
        htmls (dict): dict of plot_name->HTML for the plot
        titles (dict): dict of plot_name->plot_title
        labels (dict): dict of plot_name->plot_label
        datas (dict): dict of plot_name->[(datafile_description, data_filename), ...]
    """

    def dirname(path):
        return os.path.basename(os.path.dirname(path))

    def fill_default(dictionary, key, default_type=list):
        if key not in dictionary:
            dictionary[key] = default_type()

    j2_env = Environment(
        loader=FileSystemLoader(os.path.join(_ROOT, 'CRISPRessoReports', 'templates')),
    )
    j2_env.filters['dirname'] = dirname
    if crispresso_tool == 'batch':
        template = 'batchReport.html'
    elif crispresso_tool == 'pooled':
        template = 'pooledReport.html'
    elif crispresso_tool == 'wgs':
        template = 'wgsReport.html'
    else:
        template = 'multiReport.html'

    crispresso_data_path = os.path.relpath(
        crispresso_folder, os.path.dirname(crispresso_multi_report_file),
    )
    if crispresso_data_path == ".":
        crispresso_data_path = ""
    else:
        crispresso_data_path += "/"

    if allele_modification_heatmap_plot is None:
        allele_modification_heatmap_plot = {}
    if allele_modification_line_plot is None:
        allele_modification_line_plot = {}
    dictionaries = [
        allele_modification_heatmap_plot, allele_modification_line_plot,
    ]
    keys_and_default_types = [
        ('names', list),
        ('htmls', dict),
        ('titles', list),
        ('labels', dict),
        ('datas', dict),
    ]
    for dictionary in dictionaries:
        for key, default_type in keys_and_default_types:
            fill_default(
                dictionary,
                key,
                default_type,
            )
    if summary_plots is None:
        summary_plots={
            'names': [],
            'titles': [],
            'labels': [],
            'datas': [],
        }

    for html in sub_html_files:
        sub_html_files[html] = crispresso_data_path + sub_html_files[html]
    with open(crispresso_multi_report_file, 'w', encoding="utf-8") as outfile:
        outfile.write(render_template(
            template,
            j2_env,
            window_nuc_pct_quilts=[] if window_nuc_pct_quilts is None else window_nuc_pct_quilts,
            nuc_pct_quilts=[] if nuc_pct_quilts is None else nuc_pct_quilts,
            window_nuc_conv_plots=[] if window_nuc_conv_plots is None else window_nuc_conv_plots,
            nuc_conv_plots=[] if nuc_conv_plots is None else nuc_conv_plots,
            crispresso_data_path=crispresso_data_path,
            report_data={
                'names': summary_plots['names'],
                'titles': summary_plots['titles'],
                'labels': summary_plots['labels'],
                'datas': summary_plots['datas'],
                'htmls': [],
                'crispresso_data_path': crispresso_data_path,
            },
            run_names=run_names,
            failed_runs=failed_runs,
            sub_html_files=sub_html_files,
            report_name=report_name,
            compact_plots_to_show=[] if compact_plots_to_show is None else compact_plots_to_show,
            allele_modification_heatmap_plot_names=allele_modification_heatmap_plot['names'],
            allele_modification_heatmap_plot_htmls=allele_modification_heatmap_plot['htmls'],
            allele_modification_heatmap_plot_titles=allele_modification_heatmap_plot['titles'],
            allele_modification_heatmap_plot_labels=allele_modification_heatmap_plot['labels'],
            allele_modification_heatmap_plot_datas=allele_modification_heatmap_plot['datas'],
            allele_modification_line_plot_names=allele_modification_line_plot['names'],
            allele_modification_line_plot_htmls=allele_modification_line_plot['htmls'],
            allele_modification_line_plot_titles=allele_modification_line_plot['titles'],
            allele_modification_line_plot_labels=allele_modification_line_plot['labels'],
            allele_modification_line_plot_datas=allele_modification_line_plot['datas'],
        ))


def make_aggregate_report(
    crispresso2_info,
    report_name,
    crispresso_report_file,
    crispresso_report_folder,
    _ROOT,
    folder_arr,
    crispresso_html_reports,
    compact_plots_to_show=None,
    display_names=None,
):
    """
    Prepares information to make a report of a CRISPRessoAggregate run

    Parameters:
    crispresso2_info (dict): information from the crispresso aggregate run
    report_name (string): text to be shown at top of report
    crispresso_report_file (string): path to write report to
    crispresso_report_folder (string): path containing aggregated plots, etc.
    _ROOT (string): location of crispresso assets (images, templates, etc)
    folder_arr (arr of strings): paths to the aggregated crispresso folders
    crispresso_html_reports (dict): folder->html_path; Paths to the aggregated crispresso run html reports
    compact_plots_to_show (dict): name=>{'href': path to target(report) when user clicks on image, 'img': path to png image to show}
    display_names (dict): folder->display_name; Titles to be shown for crispresso runs 
        (if different from names_arr, e.g. if display_names have spaces or bad chars, they won't be the same as names_arr)

    Returns:
    Nothing
    """
    summary_plots = {}
    if 'summary_plot_names' in crispresso2_info['results']['general_plots']:
        summary_plots['names'] = crispresso2_info['results']['general_plots']['summary_plot_names']
    else:
        summary_plots['names'] = []
    if 'summary_plot_titles' in crispresso2_info['results']['general_plots']:
        summary_plots['titles'] = crispresso2_info['results']['general_plots']['summary_plot_titles']
    else:
        summary_plots['titles'] = {}
    if 'summary_plot_labels' in crispresso2_info['results']['general_plots']:
        summary_plots['labels'] = crispresso2_info['results']['general_plots']['summary_plot_labels']
    else:
        summary_plots['labels'] = {}
    if 'summary_plot_datas' in crispresso2_info['results']['general_plots']:
        summary_plots['datas'] = crispresso2_info['results']['general_plots']['summary_plot_datas']
    else:
        summary_plots['datas'] = {}

    allele_modification_heatmap_plot = {}
    if 'allele_modification_heatmap_plot_names' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['names'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_names']
    else:
        allele_modification_heatmap_plot['names'] = []
    if 'allele_modification_heatmap_plot_paths' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['paths'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_paths']
    else:
        allele_modification_heatmap_plot['paths'] = {}
    if 'allele_modification_heatmap_plot_titles' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['titles'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_titles']
    else:
        allele_modification_heatmap_plot['titles'] = {}
    if 'allele_modification_heatmap_plot_labels' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['labels'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_labels']
    else:
        allele_modification_heatmap_plot['labels'] = {}
    if 'allele_modification_heatmap_plot_datas' in crispresso2_info['results']['general_plots']:
        allele_modification_heatmap_plot['datas'] = crispresso2_info['results']['general_plots']['allele_modification_heatmap_plot_datas']
    else:
        allele_modification_heatmap_plot['datas'] = {}

    allele_modification_line_plot = {}
    if 'allele_modification_line_plot_names' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['names'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_names']
    else:
        allele_modification_line_plot['names'] = []
    if 'allele_modification_line_plot_paths' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['paths'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_paths']
    else:
        allele_modification_line_plot['paths'] = {}
    if 'allele_modification_line_plot_titles' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['titles'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_titles']
    else:
        allele_modification_line_plot['titles'] = {}
    if 'allele_modification_line_plot_labels' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['labels'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_labels']
    else:
        allele_modification_line_plot['labels'] = {}
    if 'allele_modification_line_plot_datas' in crispresso2_info['results']['general_plots']:
        allele_modification_line_plot['datas'] = crispresso2_info['results']['general_plots']['allele_modification_line_plot_datas']
    else:
        allele_modification_line_plot['datas'] = {}

    window_nuc_pct_quilts = []
    if 'window_nuc_pct_quilt_plot_names' in crispresso2_info['results']['general_plots']:
        window_nuc_pct_quilts = crispresso2_info['results']['general_plots']['window_nuc_pct_quilt_plot_names']
    nuc_pct_quilts = []
    if 'nuc_pct_quilt_plot_names' in crispresso2_info['results']['general_plots']:
        nuc_pct_quilts = crispresso2_info['results']['general_plots']['nuc_pct_quilt_plot_names']

    run_names = []
    sub_html_files = {}

    for folder in folder_arr:
        display_name = folder
        if display_names is not None:
            display_name = display_names[folder]

        run_names.append(display_name)
        sub_html_file = os.path.relpath(crispresso_html_reports[folder], crispresso_report_folder)
        sub_html_files[display_name] = sub_html_file
    if compact_plots_to_show is None:
        compact_plots_to_show = {}
    for compact_plot in compact_plots_to_show:
        old_href = compact_plots_to_show[compact_plot]['href']
        compact_plots_to_show[compact_plot]['href'] = os.path.relpath(old_href, crispresso_report_folder)
        old_img = compact_plots_to_show[compact_plot]['img']
        compact_plots_to_show[compact_plot]['img'] = os.path.relpath(old_img, crispresso_report_folder)

    allele_modification_heatmap_plot['htmls'] = {}
    for heatmap_plot_name, heatmap_plot_path in allele_modification_heatmap_plot['paths'].items():
        with open(heatmap_plot_path, encoding="utf-8") as fh:
            allele_modification_heatmap_plot['htmls'][heatmap_plot_name] = fh.read()

    allele_modification_line_plot['htmls'] = {}
    for line_plot_name, line_plot_path in allele_modification_line_plot['paths'].items():
        with open(line_plot_path, encoding="utf-8") as fh:
            allele_modification_line_plot['htmls'][line_plot_name] = fh.read()

    make_multi_report(
        run_names,
        sub_html_files,
        crispresso_report_file,
        crispresso_report_folder,
        _ROOT,
        report_name,
        window_nuc_pct_quilts=window_nuc_pct_quilts,
        nuc_pct_quilts=nuc_pct_quilts,
        summary_plots=summary_plots,
        compact_plots_to_show=compact_plots_to_show,
        allele_modification_heatmap_plot=allele_modification_heatmap_plot,
        allele_modification_line_plot=allele_modification_line_plot,
    )
