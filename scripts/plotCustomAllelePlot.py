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
    parser.add_argument("--max_rows","--max_rows_alleles_around_cut_to_plot",type=float,help="Maximum number of rows to plot")
    parser.add_argument("--plot_cut_point",help="If set, a line at the cut point will be plotted.",action="store_true")
    parser.add_argument("--save_png",help="If set, pngs will also be produced (as well as pdfs).",action="store_true")
    parser.add_argument("--plot_left",help="Number of bases to plot to the left of the cut site",type=int,default=20)
    parser.add_argument("--plot_right",help="Number of bases to plot to the right of the cut site",type=int,default=20)

    args = parser.parse_args()

    plot_alleles_tables_from_folder(args.CRISPResso2_folder,args.output_root,MIN_FREQUENCY=args.min_freq,MAX_N_ROWS=args.max_rows,SAVE_ALSO_PNG=args.save_png,plot_cut_point=args.plot_cut_point,plot_left=args.plot_left,plot_right=args.plot_right)

def arrStr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]

def get_row_around_cut_assymetrical(row,cut_point,plot_left,plot_right):
    cut_idx=row['ref_positions'].index(cut_point)
    return row['Aligned_Sequence'][cut_idx-plot_left+1:cut_idx+plot_right+1],row['Reference_Sequence'][cut_idx-plot_left+1:cut_idx+plot_right+1],row['Read_Status']=='UNMODIFIED',row['n_deleted'],row['n_inserted'],row['n_mutated'],row['#Reads'], row['%Reads']

def get_dataframe_around_cut_assymetrical(df_alleles, cut_point,plot_left,plot_right,collapse_by_sequence=True):
    if df_alleles.shape[0] == 0:
        return df_alleles
    ref1 = df_alleles['Reference_Sequence'].iloc[0]
    ref1 = ref1.replace('-','')
    if (cut_point + plot_right + 1 > len(ref1)):
        raise(BadParameterException('The plotting window cannot extend past the end of the amplicon. Amplicon length is ' + str(len(ref1)) + ' but plot extends to ' + str(cut_point+plot_right+1)))

    df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut_assymetrical(row,cut_point,plot_left,plot_right),axis=1).values),
                    columns=['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated','#Reads','%Reads'])

    df_alleles_around_cut=df_alleles_around_cut.groupby(['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by='%Reads',inplace=True,ascending=False)
    df_alleles_around_cut['Unedited']=df_alleles_around_cut['Unedited']>0
    return df_alleles_around_cut

def plot_alleles_tables_from_folder(crispresso_output_folder,fig_filename_root,plot_left=20,plot_right=20,MIN_FREQUENCY=None,MAX_N_ROWS=None,SAVE_ALSO_PNG=False,custom_colors=None,plot_cut_point=True,sgRNA_intervals=None,sgRNA_names=None,sgRNA_mismatches=None):
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

    ref_names = crispresso2_info['ref_names']
    refs = crispresso2_info['refs']
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
            ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_left+1:cut_point+plot_right+1]

            df_alleles_around_cut=get_dataframe_around_cut_assymetrical(df_alleles, cut_point, plot_left, plot_right)
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
            plot_alleles_table(ref_seq_around_cut,df_alleles=df_alleles_around_cut,fig_filename_root=fig_filename_root,cut_point_ind=cut_point-new_sel_cols_start, MIN_FREQUENCY=MIN_FREQUENCY,MAX_N_ROWS=MAX_N_ROWS,SAVE_ALSO_PNG=SAVE_ALSO_PNG,plot_cut_point=plot_cut_point,sgRNA_intervals=new_sgRNA_intervals,sgRNA_names=sgRNA_names,sgRNA_mismatches=sgRNA_mismatches,annotate_wildtype_allele=crispresso2_info['running_info']['args'].annotate_wildtype_allele)

            plot_count += 1
    print('Plotted ' + str(plot_count) + ' plots')

def plot_alleles_table(reference_seq,df_alleles,fig_filename_root,MIN_FREQUENCY=0.5,MAX_N_ROWS=100,SAVE_ALSO_PNG=False,plot_cut_point=True,cut_point_ind=None,sgRNA_intervals=None,sgRNA_names=None,sgRNA_mismatches=None,custom_colors=None,annotate_wildtype_allele='****',):
    """
    plots an allele table for a dataframe with allele frequencies
    input:
    reference_seq: the reference amplicon sequence to plot
    df_alleles: merged dataframe (should include columns "#Reads','%Reads')
    fig_filename: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    plot_cut_point: if false, won't draw 'predicted cleavage' line
    cut_point_ind: index to plot cut point at
    sgRNA_intervals: locations where sgRNA is located
    sgRNA_mismatches: array (for each sgRNA_interval) of locations in sgRNA where there are mismatches
    sgRNA_names: array (for each sgRNA_interval) of names of sgRNAs (otherwise empty)
    custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    annotate_wildtype_allele: string to add to the end of the wildtype allele (e.g. ** or '')
    """
    X,annot,y_labels,insertion_dict,per_element_annot_kws,is_reference = CRISPRessoPlot.prep_alleles_table(df_alleles,reference_seq,MAX_N_ROWS,MIN_FREQUENCY)
    if annotate_wildtype_allele != '':
        for ix, is_ref in enumerate(is_reference):
            if is_ref:
                y_labels[ix] += annotate_wildtype_allele
    plot_alleles_heatmap(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,SAVE_ALSO_PNG,plot_cut_point,cut_point_ind,sgRNA_intervals,sgRNA_names,sgRNA_mismatches,custom_colors)

def plot_alleles_heatmap(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,SAVE_ALSO_PNG=False,plot_cut_point=True,cut_point_ind=None,sgRNA_intervals=None,sgRNA_names=None,sgRNA_mismatches=None,custom_colors=None):
    """
    Plots alleles in a heatmap (nucleotides color-coded for easy visualization)
    input:
    -reference_seq: sequence of reference allele to plot
    -fig_filename: figure filename to plot (not including '.pdf' or '.png')
    -X: list of numbers representing nucleotides of the allele
    -annot: list of nucleotides (letters) of the allele
    -y_labels: list of labels for each row/allele
    -insertion_dict: locations of insertions -- red squares will be drawn around these
    -per_element_annot_kws: annotations for each cell (e.g. bold for substitutions, etc.)
    -SAVE_ALSO_PNG: whether to write png file as well
    -plot_cut_point: if false, won't draw 'predicted cleavage' line
    -cut_point_ind: index to plot cut point at
    -sgRNA_intervals: locations where sgRNA is located
    -sgRNA_mismatches: array (for each sgRNA_interval) of locations in sgRNA where there are mismatches
    -sgRNA_names: array (for each sgRNA_interval) of names of sgRNAs (otherwise empty)
    -custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    """
    plot_nuc_len=len(reference_seq)

    # make a color map of fixed colors
    alpha=0.4
    A_color=CRISPRessoPlot.get_nuc_color('A',alpha)
    T_color=CRISPRessoPlot.get_nuc_color('T',alpha)
    C_color=CRISPRessoPlot.get_nuc_color('C',alpha)
    G_color=CRISPRessoPlot.get_nuc_color('G',alpha)
    INDEL_color = CRISPRessoPlot.get_nuc_color('N',alpha)

    if custom_colors is not None:
        if 'A' in custom_colors:
            A_color = custom_colors['A']
        if 'T' in custom_colors:
            T_color = custom_colors['T']
        if 'C' in custom_colors:
            C_color = custom_colors['C']
        if 'G' in custom_colors:
            G_color = custom_colors['G']
        if 'N' in custom_colors:
            INDEL_color = custom_colors['N']

    dna_to_numbers={'-':0,'A':1,'T':2,'C':3,'G':4,'N':5}
    seq_to_numbers= lambda seq: [dna_to_numbers[x] for x in seq]

    cmap = colors_mpl.ListedColormap([INDEL_color, A_color,T_color,C_color,G_color,INDEL_color])

    #ref_seq_around_cut=reference_seq[max(0,cut_point-plot_nuc_len/2+1):min(len(reference_seq),cut_point+plot_nuc_len/2+1)]

#    print('per element anoot kws: ' + per_element_annot_kws)
    if len(per_element_annot_kws) > 1:
        per_element_annot_kws=np.vstack(per_element_annot_kws[::-1])
    else:
        per_element_annot_kws=np.array(per_element_annot_kws)
    ref_seq_hm=np.expand_dims(seq_to_numbers(reference_seq),1).T
    ref_seq_annot_hm=np.expand_dims(list(reference_seq),1).T

    NEW_SEABORN=np.sum(np.array(map(int,sns.__version__.split('.')))*(100,10,1))>= 80

    if NEW_SEABORN:
        annot=annot[::-1]
        X=X[::-1]

    N_ROWS=len(X)
    N_COLUMNS=plot_nuc_len

    if N_ROWS < 1:
        fig=plt.figure()
        ax = fig.add_subplot(111)
        plt.text(0.5, 0.5,'No Alleles',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        ax.set_clip_on(False)

        plt.savefig(fig_filename_root+'.pdf',bbox_inches='tight')
        if SAVE_ALSO_PNG:
            plt.savefig(fig_filename_root+'.png',bbox_inches='tight')
        plt.close()
        return

    sgRNA_rows = []
    num_sgRNA_rows = 0

    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        sgRNA_rows = CRISPRessoPlot.get_rows_for_sgRNA_annotation(sgRNA_intervals,plot_nuc_len)
        num_sgRNA_rows = max(sgRNA_rows) + 1
        fig=plt.figure(figsize=(plot_nuc_len*0.3,(N_ROWS+1 + num_sgRNA_rows)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+2,N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+2,N_COLUMNS)
        #ax_hm_ref heatmap for the reference
        ax_hm_ref=plt.subplot(gs1[0:1, :])
        ax_hm=plt.subplot(gs2[2:, :])
    else:
        fig=plt.figure(figsize=(plot_nuc_len*0.3,(N_ROWS+1)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+1,N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+1,N_COLUMNS)
        #ax_hm_ref heatmap for the reference
        ax_hm_ref=plt.subplot(gs1[0, :])
        ax_hm=plt.subplot(gs2[1:, :])


    CRISPRessoPlot.custom_heatmap(ref_seq_hm,annot=ref_seq_annot_hm,annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm_ref,vmin=0,vmax=5,square=True)
    CRISPRessoPlot.custom_heatmap(X,annot=np.array(annot),annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm,vmin=0,vmax=5,square=True, per_element_annot_kws=per_element_annot_kws)

    ax_hm.yaxis.tick_right()
    ax_hm.yaxis.set_ticklabels(y_labels[::-1],rotation=True,va='center')
    ax_hm.xaxis.set_ticks([])

    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        this_sgRNA_y_start = -1*num_sgRNA_rows
        this_sgRNA_y_height = num_sgRNA_rows - 0.3
        CRISPRessoPlot.add_sgRNA_to_ax(ax_hm_ref,sgRNA_intervals,sgRNA_y_start=this_sgRNA_y_start,sgRNA_y_height=this_sgRNA_y_height,amp_len=plot_nuc_len,font_size='small',clip_on=False,sgRNA_names=sgRNA_names,sgRNA_mismatches=sgRNA_mismatches,x_offset=0,label_at_zero=True,sgRNA_rows=sgRNA_rows)

# todo -- add sgRNAs below reference plot
#    if sgRNA_intervals:
#        ax_hm_anno=plt.subplot(gs3[2, :])
#        sgRNA_y_start = 0.3
##        sgRNA_y_height = 0.1
#        sgRNA_y_height = 10
#        min_sgRNA_x = None
#        for idx,sgRNA_int in enumerate(sgRNA_intervals):
#            ax_hm_anno.add_patch(
#                patches.Rectangle((2+sgRNA_int[0], sgRNA_y_start), 1+sgRNA_int[1]-sgRNA_int[0], sgRNA_y_height,facecolor=(0,0,0,0.15))
#                )
#            #set left-most sgrna start
#            if not min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#            if sgRNA_int[0] < min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#        ax_hm_anno.text(2+min_sgRNA_x,sgRNA_y_start + sgRNA_y_height/2,'sgRNA ',horizontalalignment='right',verticalalignment='center')

    #print lines


    #create boxes for ins
    for idx,lss in insertion_dict.iteritems():
        for ls in lss:
            ax_hm.add_patch(patches.Rectangle((ls[0],N_ROWS-idx-1),ls[1]-ls[0],1,linewidth=3,edgecolor='r',fill=False))

    #cut point vertical line
    if plot_cut_point:
        if cut_point_ind is None:
            ax_hm.vlines([plot_nuc_len/2],*ax_hm.get_ylim(),linestyles='dashed')
        else:
            ax_hm.vlines(cut_point_ind,*ax_hm.get_ylim(),linestyles='dashed')

    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(['Reference'],rotation=True,va='center')



    gs2.update(left=0,right=1, hspace=0.05,wspace=0,top=1*(((N_ROWS)*1.13))/(N_ROWS))
    gs1.update(left=0,right=1, hspace=0.05,wspace=0,)

    sns.set_context(rc={'axes.facecolor':'white','lines.markeredgewidth': 1,'mathtext.fontset' : 'stix','text.usetex':True,'text.latex.unicode':True} )

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                    mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'),ms=18),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='r', marker='s',ms=8,markeredgewidth=2.5),
              matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='black', marker='_',ms=2,)]
    descriptions=['Substitutions','Insertions','Deletions']

    if plot_cut_point:
        proxies.append(
              matplotlib.lines.Line2D([0], [1], linestyle='--',c='black',ms=6))
        descriptions.append('Predicted cleavage position')

    #ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)
    lgd = ax_hm.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='upper center', bbox_to_anchor=(0.5, 0),ncol=1,fancybox=True,shadow=False)

    plt.savefig(fig_filename_root+'.pdf',bbox_inches='tight',bbox_extra_artists=(lgd,))
    if SAVE_ALSO_PNG:
        plt.savefig(fig_filename_root+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,))
    plt.close()


if __name__ == "__main__":
    # execute only if run as a script
    main()
