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
from CRISPResso2 import CRISPRessoShared
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors as colors_mpl
import matplotlib.patches as patches

def main():
    parser = argparse.ArgumentParser(description="Plot custom allele plots")
    parser.add_argument("-f","--CRISPResso2_folder",type=str,help="CRISPResso output folder containing finished analysis",required=True)
    parser.add_argument("-o","--output_root",type=str,help="Plot output root (should not include '.pdf' or '.png')",required=True)
    parser.add_argument("--min_freq","--min_frequency_alleles_around_cut_to_plot",type=float,help="Minimum frequency of alleles to plot")
    parser.add_argument("--max_rows","--max_rows_alleles_around_cut_to_plot",type=int,help="Maximum number of rows to plot")
    parser.add_argument("--plot_cut_point",help="If set, a line at the cut point will be plotted.",action="store_true")
    parser.add_argument("--save_png",help="If set, pngs will also be produced (as well as pdfs).",action="store_true")
    parser.add_argument("--plot_left",help="Number of bases to plot to the left of the cut site",type=int,default=20)
    parser.add_argument("--plot_right",help="Number of bases to plot to the right of the cut site",type=int,default=20)
    parser.add_argument("--plot_center",help="Center of plot. If set, plots for guide RNAs will not be generated -- only a plot centered at this position will be plotted.",type=int,default=None)
    
    # CRISPRessoPro params
    parser.add_argument('--use_matplotlib', action='store_true',
                        help='Use matplotlib for plotting instead of plotly/d3 when CRISPRessoPro is installed')

    args = parser.parse_args()
    if args.use_matplotlib or not CRISPRessoShared.is_C2Pro_installed():
        from CRISPResso2 import CRISPRessoPlot
    else:
        from CRISPRessoPro import plot as CRISPRessoPlot

    plot_alleles_tables_from_folder(args.CRISPResso2_folder,args.output_root,MIN_FREQUENCY=args.min_freq,MAX_N_ROWS=args.max_rows,SAVE_ALSO_PNG=args.save_png,plot_cut_point=args.plot_cut_point,plot_left=args.plot_left,plot_right=args.plot_right,plot_center=args.plot_center)

def arrStr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]

def get_row_around_cut_asymmetrical(row,cut_point,plot_left,plot_right):
    cut_idx=row['ref_positions'].index(cut_point)
    return row['Aligned_Sequence'][cut_idx-plot_left+1:cut_idx+plot_right+1],row['Reference_Sequence'][cut_idx-plot_left+1:cut_idx+plot_right+1],row['Read_Status']=='UNMODIFIED',row['n_deleted'],row['n_inserted'],row['n_mutated'],row['#Reads'], row['%Reads']

def get_dataframe_around_cut_asymmetrical(df_alleles, cut_point,plot_left,plot_right,collapse_by_sequence=True):
    if df_alleles.shape[0] == 0:
        return df_alleles
    ref1 = df_alleles['Reference_Sequence'].iloc[0]
    ref1 = ref1.replace('-','')
    if (cut_point + plot_right + 1 > len(ref1)):
        raise(CRISPRessoShared.BadParameterException('The plotting window cannot extend past the end of the amplicon. Amplicon length is ' + str(len(ref1)) + ' but plot extends to ' + str(cut_point+plot_right+1)))

    df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut_asymmetrical(row,cut_point,plot_left,plot_right),axis=1).values),
                    columns=['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated','#Reads','%Reads'])

    df_alleles_around_cut=df_alleles_around_cut.groupby(['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by=['#Reads', 'Aligned_Sequence', 'Reference_Sequence'], inplace=True, ascending=[False, True, True])
    df_alleles_around_cut['Unedited']=df_alleles_around_cut['Unedited']>0
    return df_alleles_around_cut

def plot_alleles_tables_from_folder(crispresso_output_folder,fig_filename_root,plot_left=20,plot_right=20,plot_center=None,MIN_FREQUENCY=None,MAX_N_ROWS=None,SAVE_ALSO_PNG=False,custom_colors=None,plot_cut_point=True,sgRNA_intervals=None,sgRNA_names=None,sgRNA_mismatches=None):
    """
    Plots an allele table plot from a completed CRISPResso run but plots a specified number of bases left and right from the cut site
    This function is only used for one-off plotting purposes and not for the general CRISPResso analysis
    Important: The run must have been run with the --write_detailed_allele_table parameter
    crispresso_output_folder: completed analysis crispresso2 output folder
    fig_filename_root: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    plot_cut_point: if false, won't draw 'predicted cleavage' line
    plot_left: number of bases left to plot from cut point
    plot_right: number of bases right to plot from cut point
    """
    crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_output_folder)

    if not crispresso2_info['running_info']['args'].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')


    if MIN_FREQUENCY is None:
        MIN_FREQUENCY = crispresso2_info['running_info']['args'].min_frequency_alleles_around_cut_to_plot
    if MAX_N_ROWS is None:
        MAX_N_ROWS = crispresso2_info['running_info']['args'].max_rows_alleles_around_cut_to_plot

    plot_count = 0

    z = zipfile.ZipFile(os.path.join(crispresso_output_folder,crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
    zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
    df_alleles = pd.read_csv(zf,sep="\t")
    full_len = df_alleles['#Reads'].sum()
    df_alleles['ref_positions'] = df_alleles['ref_positions'].apply(arrStr_to_arr)

    ref_names = crispresso2_info['results']['ref_names']
    refs = crispresso2_info['results']['refs']
    for ref_name in ref_names:
        sgRNA_sequences = refs[ref_name]['sgRNA_sequences']
        sgRNA_cut_points = refs[ref_name]['sgRNA_cut_points']
        sgRNA_plot_cut_points = refs[ref_name]['sgRNA_plot_cut_points']
        sgRNA_intervals = refs[ref_name]['sgRNA_intervals']
        sgRNA_names = refs[ref_name]['sgRNA_names']
        sgRNA_mismatches = refs[ref_name]['sgRNA_mismatches']
        sgRNA_plot_idxs = refs[ref_name]['sgRNA_plot_idxs']

        reference_seq = refs[ref_name]['sequence']

        if plot_center is not None:
            sgRNA_label = 'custom'

            cut_point = plot_center
            plot_cut_point = plot_center
            ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_left+1:cut_point+plot_right+1]

            df_alleles_around_cut=get_dataframe_around_cut_asymmetrical(df_alleles, cut_point, plot_left, plot_right)
            this_allele_count = len(df_alleles_around_cut.index)
            if this_allele_count < 1:
                print('No reads found for ' + ref_name)
                continue
            this_reads_count = df_alleles_around_cut['#Reads'].sum()
            print('Plotting ' + str(this_reads_count) + ' reads for ' + ref_name)

            new_sgRNA_intervals = []
            #adjust coordinates of sgRNAs
            new_sel_cols_start = cut_point - plot_left
            for (int_start, int_end) in refs[ref_name]['sgRNA_intervals']:
                new_sgRNA_intervals += [(int_start - new_sel_cols_start - 1,int_end - new_sel_cols_start - 1)]

            fig_filename_root = fig_filename_root+"_"+ref_name+"_"+sgRNA_label
            CRISPRessoPlot.plot_alleles_table(ref_seq_around_cut,
                                              df_alleles=df_alleles_around_cut,
                                              fig_filename_root=fig_filename_root,
                                              cut_point_ind=cut_point-new_sel_cols_start,
                                              custom_colors=custom_colors,
                                              MIN_FREQUENCY=MIN_FREQUENCY,
                                              MAX_N_ROWS=MAX_N_ROWS,
                                              SAVE_ALSO_PNG=SAVE_ALSO_PNG,
                                              plot_cut_point=plot_cut_point,
                                              sgRNA_intervals=new_sgRNA_intervals,
                                              sgRNA_names=sgRNA_names,
                                              sgRNA_mismatches=sgRNA_mismatches,
                                              annotate_wildtype_allele=crispresso2_info['running_info']['args'].annotate_wildtype_allele)

            plot_count += 1
        else:
            for ind,sgRNA in enumerate(sgRNA_sequences):
                sgRNA_label = sgRNA # for file names
                if sgRNA_names[ind] != "":
                    sgRNA_label = sgRNA_names[ind]

                cut_point = sgRNA_cut_points[ind]
                plot_cut_point = sgRNA_plot_cut_points[ind]
                plot_idxs = sgRNA_plot_idxs[ind]
                ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_left+1:cut_point+plot_right+1]

                df_alleles_around_cut=get_dataframe_around_cut_asymmetrical(df_alleles, cut_point, plot_left, plot_right)
                this_allele_count = len(df_alleles_around_cut.index)
                if this_allele_count < 1:
                    print('No reads found for ' + ref_name)
                    continue
                this_reads_count = df_alleles_around_cut['#Reads'].sum()
                print('Plotting ' + str(this_reads_count) + ' reads for ' + ref_name)

                new_sgRNA_intervals = []
                #adjust coordinates of sgRNAs
                new_sel_cols_start = cut_point - plot_left
                for (int_start, int_end) in refs[ref_name]['sgRNA_intervals']:
                    new_sgRNA_intervals += [(int_start - new_sel_cols_start - 1,int_end - new_sel_cols_start - 1)]

                fig_filename_root = fig_filename_root+"_"+ref_name+"_"+sgRNA_label
                CRISPRessoPlot.plot_alleles_table(ref_seq_around_cut,
                                                  df_alleles=df_alleles_around_cut,
                                                  fig_filename_root=fig_filename_root,
                                                  cut_point_ind=cut_point-new_sel_cols_start,
                                                  custom_colors=custom_colors,
                                                  MIN_FREQUENCY=MIN_FREQUENCY,
                                                  MAX_N_ROWS=MAX_N_ROWS,
                                                  SAVE_ALSO_PNG=SAVE_ALSO_PNG,
                                                  plot_cut_point=plot_cut_point,
                                                  sgRNA_intervals=new_sgRNA_intervals,
                                                  sgRNA_names=sgRNA_names,
                                                  sgRNA_mismatches=sgRNA_mismatches,
                                                  annotate_wildtype_allele=crispresso2_info['running_info']['args'].annotate_wildtype_allele)

                plot_count += 1
    print('Plotted ' + str(plot_count) + ' plots')

if __name__ == "__main__":
    # execute only if run as a script
    main()

