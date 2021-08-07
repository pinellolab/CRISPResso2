'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''

import argparse
import os
import numpy as np
import pandas as pd
import zipfile
from CRISPResso2 import CRISPRessoPlot
from CRISPResso2 import CRISPRessoShared

def main():
    parser = argparse.ArgumentParser(description="Plot ambiguous reads")
    parser.add_argument("-f","--CRISPResso2_folder",type=str,help="CRISPResso output folder to plot ambiguous reads for",required=True)
    parser.add_argument("-o","--output_root",type=str,help="Plot output root (should not include '.pdf' or '.png')",required=True)
    parser.add_argument("--min_freq","--min_frequency_alleles_around_cut_to_plot",type=float,help="Minimum frequency of alleles to plot")
    parser.add_argument("--max_rows","--max_rows_alleles_around_cut_to_plot",type=float,help="Maximum number of rows to plot")
    parser.add_argument("--plot_cut_point",help="If set, a line at the cut point will be plotted.",action="store_true")
    parser.add_argument("--save_png",help="If set, pngs will also be produced (as well as pdfs).",action="store_true")

    args = parser.parse_args()

    plot_ambiguous_alleles_tables_from_folder(args.CRISPResso2_folder,args.output_root,MIN_FREQUENCY=args.min_freq,MAX_N_ROWS=args.max_rows,SAVE_ALSO_PNG=args.save_png,plot_cut_point=args.plot_cut_point)

def arrStr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]

def plot_ambiguous_alleles_tables_from_folder(crispresso_output_folder,fig_filename_root,MIN_FREQUENCY=None,MAX_N_ROWS=None,SAVE_ALSO_PNG=False,custom_colors=None,plot_cut_point=True,sgRNA_intervals=None,sgRNA_names=None,sgRNA_mismatches=None):
    """
    Plots an allele table plot of ambiguous alleles from a completed CRISPResso run
    This function is only used for one-off plotting purposes and not for the general CRISPResso analysis
    Important: The run must have been run with the --write_detailed_allele_table parameter
    Ambiguous reads align to multiple reference amplicons with the same score
    In this function, ambiguous reads are filtered from the allele tables and the allele plots for these ambiguous reads are plotted
    Note that each ambiguous read is assigned to a reference (usually the first one) and mutations/indels are plotted in relation to this reference sequence.
    crispresso_output_folder: completed analysis crispresso2 output folder
    fig_filename_root: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    plot_cut_point: if false, won't draw 'predicted cleavage' line
    example:
    """
    crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_output_folder)

    if not crispresso2_info['running_info']['args'].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    if MIN_FREQUENCY is None:
        MIN_FREQUENCY = crispresso2_info['running_info']['args'].min_frequency_alleles_around_cut_to_plot
    if MAX_N_ROWS is None:
        MAX_N_ROWS = crispresso2_info['running_info']['args'].max_rows_alleles_around_cut_to_plot

    plot_count = 0

    z = zipfile.ZipFile(os.path.join(crispresso_output_folder, crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
    zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
    df_alleles = pd.read_csv(zf,sep="\t")
    full_len = df_alleles['#Reads'].sum()
    df_alleles['ref_positions'] = df_alleles['ref_positions'].apply(arrStr_to_arr)

    #pd.set_option('display.max_columns', None)
    #print(df_alleles.head())
    df_ambiguous = df_alleles[df_alleles['Reference_Name'].str.contains('AMBIGUOUS')]
    ambig_len = df_ambiguous['#Reads'].sum()

    print("Filtered to " + str(ambig_len) + "/" + str(full_len) + " ambiguous reads")

    ref_names = crispresso2_info['results']['ref_names']
    refs = crispresso2_info['results']['refs']
    print("Ambiguous alleles will be plotted against to the sequence of the first reference sequence ("+ref_names[0]+")")
    for ref_name in ref_names:
        sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
        sgRNA_cut_points = refs[ref_name]['sgRNA_cut_points']
        sgRNA_plot_cut_points = refs[ref_name]['sgRNA_plot_cut_points']
        sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
        sgRNA_names = refs[ref_name]['sgRNA_names']
        sgRNA_mismatches = refs[ref_name]['sgRNA_mismatches']
        sgRNA_plot_idxs = refs[ref_name]['sgRNA_plot_idxs']

        reference_seq = refs[ref_name]['sequence']

        for ind,sgRNA in enumerate(sgRNA_sequences):
            sgRNA_label = sgRNA # for file names
            if sgRNA_names[ind] != "":
                sgRNA_label = sgRNA_names[ind]

            cut_point = sgRNA_cut_points[ind]
            plot_cut_point = sgRNA_plot_cut_points[ind]
            plot_idxs = sgRNA_plot_idxs[ind]
            plot_half_window = max(1,crispresso2_info['running_info']['args'].plot_window_size)
            ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_half_window+1:cut_point+plot_half_window+1]

            ambiguous_ref_name = "AMBIGUOUS_"+ref_name
            df_alleles_around_cut=CRISPRessoShared.get_dataframe_around_cut(df_alleles.loc[df_alleles['Reference_Name'] == ambiguous_ref_name],cut_point,plot_half_window)
            this_ambig_allele_count = len(df_alleles_around_cut.index)
            if this_ambig_allele_count < 1:
                print('No ambiguous reads found for ' + ref_name)
                continue
            this_ambig_count = df_alleles_around_cut['#Reads'].sum()
            print('Plotting ' + str(this_ambig_count) + ' ambiguous reads for ' + ref_name)

            new_sgRNA_intervals = []
            #adjust coordinates of sgRNAs
            new_sel_cols_start = cut_point - plot_half_window
            for (int_start, int_end) in refs[ref_name]['sgRNA_intervals']:
                new_sgRNA_intervals += [(int_start - new_sel_cols_start - 1,int_end - new_sel_cols_start - 1)]
            fig_filename_root = fig_filename_root+"_"+ref_name+"_"+sgRNA_label
            CRISPRessoPlot.plot_alleles_table(ref_seq_around_cut,df_alleles=df_alleles_around_cut,fig_filename_root=fig_filename_root, MIN_FREQUENCY=MIN_FREQUENCY,MAX_N_ROWS=MAX_N_ROWS,SAVE_ALSO_PNG=SAVE_ALSO_PNG,plot_cut_point=plot_cut_point,sgRNA_intervals=new_sgRNA_intervals,sgRNA_names=sgRNA_names,sgRNA_mismatches=sgRNA_mismatches,annotate_wildtype_allele=crispresso2_info['running_info']['args'].annotate_wildtype_allele)

            plot_count += 1
    print('Plotted ' + str(plot_count) + ' plots')

if __name__ == "__main__":
    # execute only if run as a script
    main()
