'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from collections import defaultdict
import re
from matplotlib import colors as colors_mpl
import seaborn as sns
sns.set_style("white")
sns.set_context('poster')
sns.set(font_scale=2.2)




def setMatplotlibDefaults():
    font = {'size'   : 22}
    matplotlib.rc('font', **font)
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams['axes.facecolor'] = 'white'

setMatplotlibDefaults()


def get_nuc_color(nuc,alpha):
    get_color=lambda x,y,z: (x/255.0,y/255.0,z/255.0,alpha)
    if nuc == "A":
        return get_color(127,201,127)
    elif nuc == "T":
        return get_color(190,174,212)
    elif nuc == "C":
        return get_color(253,192,134)
    elif nuc == "G":
        return get_color(255,255,153)
    elif nuc == "N":
        return get_color(200,200,200)
    elif nuc == "INS":
#        return get_color(185,219,253)
#        return get_color(177,125,76)
        return get_color(193,129,114)
    elif nuc == "DEL":
        #return get_color(177,125,76)
#        return get_color(202,109,87)
        return get_color(193,129,114)
    elif nuc == "-":
        #return get_color(177,125,76)
#        return get_color(202,109,87)
        return get_color(30,30,30)
    else: #return a random color (that is based on the nucleotide given)
        charSum = 0
        for char in nuc.upper():
            thisval = ord(char) - 65 #'A' is 65
            if thisval < 0 or thisval > 90:
                thisval = 0
            charSum += thisval
        charSum = (charSum/len(nuc))/90.0

        return (charSum,(1-charSum),(2*charSum*(1-charSum)))

def get_color_lookup(nucs,alpha):
    colorLookup = {}
    for nuc in nucs:
        colorLookup[nuc] = get_nuc_color(nuc,alpha)
    return colorLookup

def plot_nucleotide_quilt(nuc_pct_df,mod_pct_df,fig_filename_root,save_also_png=False,sgRNA_intervals=None,min_text_pct=0.5,max_text_pct=0.95,quantification_window_idxs=None):
    """Plots a nucleotide quilt with each square showing the percentage of each base at that position in the reference
    nuc_pct_df: dataframe with percents of each base (ACTGN-) at each position
    mod_pct_df: dataframe with percents of modifications at each position (this function uses 'Insertions_Left' to plot insertions)
    fig_filename_root: filename root (will add .pdf or .png)
    save_also_png: whether png should also be saved
    sgRNA_intervals: ranges for sgRNA annotation on plot
    quantification_window_idxs: indices for quantification window annotation on plot
    min_text_pct: add text annotation if the percent is greater than this number
    max_text_pct: add text annotation if the percent is less than this number
    """

    plotPct = 0.9 #percent of vertical space to plot in (the rest will be white)
    min_plot_pct = 0.01 #if value is less than this, it won't plot the rectangle (with white boundary)

    if nuc_pct_df.iloc[1,2] > 1:
        raise Exception('Expecting nucleotide percentage. Instead, got numbers in nuc_pct_df: ' + str(nuc_pct_df.iloc[1,5]))

    nrows = nuc_pct_df.shape[0]
    amp_len = nuc_pct_df.shape[1] - 2 #Batch, Nucleotide, nuc1, nuc2, nuc3 ...
    nucs = nuc_pct_df.Nucleotide.unique()
    nNucs = len(nucs)
    nSamples = nrows/nNucs
    samplesList = []
    for i in range(nSamples): #iterate over all samples
        sample_row_start = nNucs * i
        samplesList.append(nuc_pct_df.iloc[sample_row_start,0])

    # make a color map of fixed colors
    color_lookup = get_color_lookup(['A','T','C','G','N','INS','-'],alpha=1)

    #fig = plt.figure(figsize=(amp_len/2.0,nSamples*2))
    #fig = plt.figure(figsize=(amp_len,nSamples))
    #fig = plt.figure(figsize=(amp_len,nSamples*2))
    #fig = plt.figure(figsize=(amp_len,(nSamples+1)*2))
    fig = plt.figure(figsize=((amp_len+10)/2.0,(nSamples+1)*2))
    ax = fig.add_subplot(111)

    #remove box around plot
    for spine in plt.gca().spines.values():
        spine.set_visible(False)


    for pos_ind in range(2,amp_len+2): #iterate over all nucleotide positions in the sequence (0=Batch, 1=Nucleotide, so start at 2)
        x_start = pos_ind
        x_end = pos_ind + 1
        for i in range(nSamples): #iterate over all samples
            sample_row_start = nNucs * i
            y_start = nSamples - i
            sumPct = 0
            for nuc_ind in range(nNucs): #iterate over each nucleotide at this position in this sample
                pct = float(nuc_pct_df.iloc[sample_row_start + nuc_ind, pos_ind])
                sumPct += pct
                if pct > min_plot_pct:
                    obs_pct = pct * plotPct
                    curr_nuc = nuc_pct_df.iloc[sample_row_start + nuc_ind,1]
                    ax.add_patch(
                        patches.Rectangle((x_start, y_start), x_end-x_start, obs_pct,facecolor=color_lookup[curr_nuc],edgecolor='w')
                        )
                    if pct > min_text_pct and pct < max_text_pct:
                        ax.text(x_start+0.55,y_start + obs_pct/2.0,format(pct*100,'.1f'),horizontalalignment='center',verticalalignment='center',rotation=90)
                    y_start += obs_pct

    for pos_ind in range(2,amp_len+1): #iterate over all nucleotide positions in the sequence (0=Batch, 1=Nucleotide, so start at 2)
        x_start = pos_ind + 0.7
        x_end = pos_ind + 1.3
        for i in range(nSamples): #iterate over all samples
            sampleName = samplesList[i]

            sample_row_start = nNucs * i
            y_start = nSamples - i


            ins_pct = float(mod_pct_df.loc[(mod_pct_df['Batch'] == sampleName) &
                    (mod_pct_df['Modification'] == "Insertions_Left")].iloc[:,pos_ind])
            if ins_pct > min_plot_pct:
                obs_pct = ins_pct * plotPct
                ax.add_patch(
                    patches.Rectangle((x_start, y_start), x_end-x_start, obs_pct,facecolor=color_lookup['INS'],edgecolor='w')
                    )
                if ins_pct > min_text_pct and ins_pct < max_text_pct:
                    ax.text(x_start+0.15,y_start + obs_pct/2.0,format(ins_pct*100,'.1f'),horizontalalignment='center',verticalalignment='center',rotation=90)

    #draw black box around each sample
    for i in range(nSamples):
        y_start = nSamples - i
        ax.add_patch(
            patches.Rectangle((2, y_start), amp_len, plotPct,facecolor='none',edgecolor='black')
            )

    #draw reference sequence
    ref_y_start = 0.5
    ref_y_height = 0.4
    ref_seq = nuc_pct_df.columns.values
    for pos_ind in range(2,amp_len+2): #iterate over all nucleotide positions in the sequence (0=Batch, 1=Nucleotide, so start at 2)
        ax.add_patch(
            patches.Rectangle((pos_ind, ref_y_start), 1, ref_y_height,facecolor=color_lookup[ref_seq[pos_ind]],edgecolor='w')
            )
        ax.text(pos_ind+0.5,ref_y_start + ref_y_height/2.3,ref_seq[pos_ind],horizontalalignment='center',verticalalignment='center')

    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=False)

    ax.set_yticks([ref_y_start + ref_y_height/2.0]+[x+0.5 for x in range(1,nSamples+1)])
#    sampleLabs = list(nuc_pct_df.iloc[[((nSamples-1)-x)*nNucs for x in range(0,nSamples)],0]))
#    print(mod_pct_df)
#    sampleReadCounts = list(nuc_pct_df.iloc[[((nSamples-1)-x)*nNucs for x in range(0,nSamples)],0]))
    ax.set_yticklabels(['Reference'] + list(nuc_pct_df.iloc[[((nSamples-1)-x)*nNucs for x in range(0,nSamples)],0]))


    plot_y_start = ref_y_start

    if sgRNA_intervals:
        plot_y_start = 0
        sgRNA_y_start = 0.1
        sgRNA_y_height = 0.2
        min_sgRNA_x = None
        for idx,sgRNA_int in enumerate(sgRNA_intervals):
            this_sgRNA_start = sgRNA_int[0]
            this_sgRNA_end = sgRNA_int[1]
#            print('this sgRNA_start is ' + str(this_sgRNA_start))
#            print('this sgRNA_end is ' + str(this_sgRNA_end))
#            print(nuc_pct_df)
            ax.add_patch(
                patches.Rectangle((2+sgRNA_int[0], sgRNA_y_start), 1+sgRNA_int[1]-sgRNA_int[0], sgRNA_y_height,facecolor=(0,0,0,0.15))
                )
            #set left-most sgrna start
            if not min_sgRNA_x:
                min_sgRNA_x = sgRNA_int[0]
            if sgRNA_int[0] < min_sgRNA_x:
                min_sgRNA_x = sgRNA_int[0]
        ax.text(2+min_sgRNA_x,sgRNA_y_start + sgRNA_y_height/2,'sgRNA ',horizontalalignment='right',verticalalignment='center')

    if quantification_window_idxs is not None:
        q_win_y_start = 0
        q_win_y_height = nSamples+1

        q_list = sorted(list(quantification_window_idxs))

        lastStart = q_list[0]
        lastIdx = q_list[0]
        for idx in range(1,len(q_list)):
            if q_list[idx] == lastIdx + 1:
                lastIdx = q_list[idx]
            else:
                ax.add_patch(
                    patches.Rectangle((2+lastStart, q_win_y_start), 1+(lastIdx-lastStart), q_win_y_height,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                    )
                lastStart = q_list[idx]
                lastIdx = q_list[idx]
        ax.add_patch(
            patches.Rectangle((2+lastStart, q_win_y_start), 1+(lastIdx-lastStart), q_win_y_height,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
            )

    ax.set_xlim([2,amp_len+3])
    ax.set_ylim([plot_y_start,nSamples+1.2])

    legend_patches = []
    for nuc in nucs:
        if nuc == "-":
            continue
        if nuc == "N":
            continue
        patch = patches.Patch(color=color_lookup[nuc], label=nuc)
        legend_patches.append(patch)

    n_tab = nuc_pct_df.loc[nuc_pct_df.Nucleotide == 'N']
    n_sum = n_tab.iloc[:,2:len(n_tab.columns)].sum().sum()
    ins_tab = mod_pct_df.loc[mod_pct_df.Modification == 'Insertions_Left']
    ins_sum = ins_tab.iloc[:,2:len(ins_tab.columns)].sum().sum()
    del_tab = mod_pct_df.loc[mod_pct_df.Modification == 'Deletions']
    del_sum = del_tab.iloc[:,2:len(del_tab.columns)].sum().sum()
    if n_sum > 0:
        n_patch = patches.Patch(color=color_lookup['N'], label='N')
        legend_patches.append(n_patch)
    if ins_sum > 0:
        ins_patch = patches.Patch(color=color_lookup['INS'], label='Insertion')
        legend_patches.append(ins_patch)
    if del_sum > 0:
        del_patch = patches.Patch(color=color_lookup['-'], label='Deletion')
        legend_patches.append(del_patch)
    if quantification_window_idxs:
        q_win_patch = patches.Patch(fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
        legend_patches.append(q_win_patch)


    plt.legend(handles=legend_patches,loc='center left',ncol=1,bbox_to_anchor=(1,0.5))

    fig.savefig(fig_filename_root+'.pdf',bbox_inches='tight')
    if save_also_png:
        fig.savefig(fig_filename_root+'.png',bbox_inches='tight',pad=1)
    plt.close()

def plot_conversion_map(nuc_pct_df,fig_filename_root,conversion_nuc_from,conversion_nuc_to,save_also_png,plotPct = 0.9,min_text_pct=0.3,max_text_pct=0.9,conversion_scale_max=None,sgRNA_intervals=None):
    """
    Plots a heatmap of conversion across several sequences
    :param nuc_pct_df combined df of multiple batches
    :param plotPct: percent of vertical space to plot in (the rest will be white)
    :param min_text_pct: add text annotation if the percent is greater than this number
    :param max_text_pct: add text annotation if the percent is less than this number
    """

    from_nuc = conversion_nuc_from
    to_nuc = conversion_nuc_to
    nucs = nuc_pct_df.Nucleotide.unique()
    nNucs = len(nucs)

    nrows = nuc_pct_df.shape[0]
    amp_len = nuc_pct_df.shape[1] - 2 #batch, nucleotide, nuc1, nuc2, nuc3 ...
    nSamples = nrows/nNucs

    nuc_pct_conversion = []
    for i in range(nSamples):
        newRow = [nuc_pct_df.iloc[i*nNucs,0], from_nuc + ">" + to_nuc]
        sampleRows = nuc_pct_df.iloc[i*nNucs:i*nNucs+nNucs,:]
        sub1 = sampleRows.ix[sampleRows['Nucleotide'] == to_nuc,2:amp_len+2]
        sub2 = pd.DataFrame(sampleRows.ix[sampleRows['Nucleotide'].isin([from_nuc,to_nuc]),2:amp_len+2].sum(axis=0)).transpose()
        #conversion_pcts = sampleRows.ix[sampleRows['Nucleotide'] == to_nuc,2:amp_len+2].div(sampleRows.ix[sampleRows['Nucleotide'].isin([from_nuc,to_nuc]),2:amp_len+2].sum(axis=0))
        conversion_pcts = sub1.div(sub2.values)
        newRow.extend(conversion_pcts.values.tolist()[0])
        nuc_pct_conversion.append(newRow)

    nuc_pct_conversion_df = pd.DataFrame(data=nuc_pct_conversion, columns=nuc_pct_df.columns)

    #get max pct conversion for the color bar
    from_nuc_indices = nuc_pct_conversion_df.columns == from_nuc
    if sum(from_nuc_indices) == 0:
        print('Skipping conversion plot. No ' + from_nuc + ' in sequence')
        return()


    max_pct_conversion = 1 # default
    if sum(from_nuc_indices) == 1: #if only one row where nuc_pct_conversion is from_nuc...
        max_pct_conversion = float(nuc_pct_conversion_df.iloc[:,from_nuc_indices].max())
    elif sum(from_nuc_indices) > 1: #if multiple rows
        max_pct_conversion = nuc_pct_conversion_df.iloc[:,from_nuc_indices].max().max() #one max returns column-based max, second takes max of those
    if (max_pct_conversion < 0.01): # the min conversion perctent is 0.01. The legend axis are rounded to the nearest 0.01, so if this lower limit isn't set, the max will appear as 0.00%
        max_pct_conversion = 0.02
    if conversion_scale_max:
        max_pct_conversion = conversion_scale_max

    #set up color map
    color_map = cm.Blues
    color_map_normalization = matplotlib.colors.Normalize(vmin=0, vmax=max_pct_conversion)

    # Use contourf to create colorbar info, then clear the figure
    plt_color_bar = plt.contourf([[0,0],[0,0]],[x/100.0 for x in range(101)], cmap=color_map)#,normalize=color_map_normalization)
    plt.clf()

    # make a color map of fixed colors (for coloring reference in this example)
    color_lookup = get_color_lookup(['A','T','C','G'],alpha=1)

    fig = plt.figure(figsize=(amp_len/2.0,nSamples*2))
    ax = fig.add_subplot(111)

    #remove box around plot
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    #draw gray background behind each sample
    for i in range(nSamples):
        y_start = nSamples - i
        ax.add_patch(
            patches.Rectangle((2, y_start), amp_len, plotPct,facecolor='0.95',edgecolor='black')
            )

    ref_seq = nuc_pct_df.columns.values

    #new stuff here
    for pos_ind in range(2,amp_len+2): #iterate over all nucleotide positions in the sequence (0=batch, 1=nucleotide, so start at 2)
        ref_nuc = ref_seq[pos_ind]
        if ref_nuc != from_nuc:
            continue
        x_start = pos_ind
        x_end = pos_ind + 1
        for i in range(nSamples): #iterate over all samples
            sample_row_start = nNucs * i
            y_start = nSamples - i
            conversion_pct = nuc_pct_conversion_df.iloc[i,pos_ind]
            ax.add_patch(
                        patches.Rectangle((x_start, y_start), x_end-x_start, plotPct,facecolor=color_map(color_map_normalization(conversion_pct)),edgecolor='w')
                        )
            if conversion_pct > min_text_pct and conversion_pct < max_text_pct:
                textCol = 'k'
                if (conversion_pct/max_pct_conversion > 0.75):
                    textCol = 'w'
                ax.text(x_start+0.6,y_start + plotPct/2.0,format(conversion_pct*100,'.1f'),horizontalalignment='center',verticalalignment='center',rotation=90,color=textCol)

    #draw black box around each sample
    for i in range(nSamples):
        y_start = nSamples - i
        ax.add_patch(
            patches.Rectangle((2, y_start), amp_len, plotPct,facecolor='None',edgecolor='black')
            )

    #draw reference sequence
    ref_y_start = 0.5
    ref_y_height = 0.4
    for pos_ind in range(2,amp_len+2): #iterate over all nucleotide positions in the sequence (0=nucName, 1=batch, so start at 2)
        ax.add_patch(
            patches.Rectangle((pos_ind, ref_y_start), 1, ref_y_height,facecolor=color_lookup[ref_seq[pos_ind]],edgecolor='w')
            )
        ax.text(pos_ind+0.5,ref_y_start + ref_y_height/2.0,ref_seq[pos_ind],horizontalalignment='center',verticalalignment='center')

    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=False)
#    plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='off')

    ax.set_yticks([ref_y_start + ref_y_height/2.0]+[x+0.5 for x in range(1,nSamples+1)])
    ax.set_yticklabels(['Reference'] + list(nuc_pct_df.iloc[[((nSamples-1)-x)*nNucs for x in range(0,nSamples)],0]))

    ax.set_xlim([2,amp_len+3])
    ax.set_ylim([ref_y_height - 0.2,nSamples+1.2])

    if sgRNA_intervals:
        sgRNA_y_start = 0.3
        sgRNA_y_height = 0.1
        min_sgRNA_x = None
        for idx,sgRNA_int in enumerate(sgRNA_intervals):
            ax.add_patch(
                patches.Rectangle((2+sgRNA_int[0], sgRNA_y_start), 1+sgRNA_int[1]-sgRNA_int[0], sgRNA_y_height,facecolor=(0,0,0,0.15))
                )
            #set left-most sgrna start
            if not min_sgRNA_x:
                min_sgRNA_x = sgRNA_int[0]
            if sgRNA_int[0] < min_sgRNA_x:
                min_sgRNA_x = sgRNA_int[0]
        ax.text(2+min_sgRNA_x,sgRNA_y_start + sgRNA_y_height/2,'sgRNA ',horizontalalignment='right',verticalalignment='center')

    #legend
    #cbar_pad = max(0.05,nSamples * 0.01) #so the cbar label doesn't get written on the data for large numbers of samples
    cbar_pad = 0.05
#    cbar_shrink = min(1.0,3/float(nSamples))
    #cbar = plt.colorbar(plt_color_bar,ax=ax,pad=cbar_pad) # using the colorbar info I got from contourf above
    cbar = plt.colorbar(plt_color_bar,ax=ax,fraction=0.046, pad=0.04) # using the colorbar info I got from contourf above
    ticks_at = [x/100.0 for x in range(0,101,25)]
    cbar.set_ticks(ticks_at)
    cbar.set_ticklabels([ "{:0.2f}".format(x*max_pct_conversion*100) for x in ticks_at ])
    #cbar.ax.text(-1,0.5,'Percentage %s to %s conversion'%(from_nuc,to_nuc), rotation=90,verticalalignment='center',horizontalalignment='right')
    cbar.ax.text(0.02,0.5,'Percentage %s to %s conversion'%(from_nuc,to_nuc), rotation=90,verticalalignment='center',horizontalalignment='right')

    #title
    plt.title('%s to %s Conversion Percent'%(from_nuc,to_nuc))

    #throws error here in tight_layout... if the selected reference is too small  (e.g. 5bp)
#    plt.tight_layout()
    fig.savefig(fig_filename_root+'.pdf',bbox_inches='tight')
    if save_also_png:
        fig.savefig(fig_filename_root+'.png',bbox_inches='tight',pad=1)
    plt.close()


def plot_subs_across_ref(ref_len, ref_seq, ref_name, ref_count, all_substitution_base_vectors, plot_title, fig_filename_root, save_also_png, quantification_window_idxs=None):
    """
    Plots substitutions across the reference sequece - each position on the x axis reprsents a nucleotide in the reference
    bars at each x posion show the number of times the reference nucleotide was substituted for another reference
    """

    fig, ax = plt.subplots(figsize=(16,8))
    ind = np.arange(ref_len)

    alph = ['A','C','G','T','N']
    color_lookup = get_color_lookup(alph,alpha=1)

    pA = ax.bar(ind,all_substitution_base_vectors[ref_name+"_A"], color=color_lookup['A'])
    pC = ax.bar(ind,all_substitution_base_vectors[ref_name+"_C"], color=color_lookup['C'],bottom=all_substitution_base_vectors[ref_name+"_A"])
    pG = ax.bar(ind,all_substitution_base_vectors[ref_name+"_G"], color=color_lookup['G'],bottom=all_substitution_base_vectors[ref_name+"_A"]+all_substitution_base_vectors[ref_name+"_C"])
    pT = ax.bar(ind,all_substitution_base_vectors[ref_name+"_T"], color=color_lookup['T'],bottom=all_substitution_base_vectors[ref_name+"_A"]+all_substitution_base_vectors[ref_name+"_C"]+all_substitution_base_vectors[ref_name+"_G"])
    pN = ax.bar(ind,all_substitution_base_vectors[ref_name+"_N"], color=color_lookup['N'],bottom=all_substitution_base_vectors[ref_name+"_A"]+all_substitution_base_vectors[ref_name+"_C"]+all_substitution_base_vectors[ref_name+"_G"]+all_substitution_base_vectors[ref_name+"_T"])
    tots = all_substitution_base_vectors[ref_name+"_N"]+all_substitution_base_vectors[ref_name+"_A"]+all_substitution_base_vectors[ref_name+"_C"]+all_substitution_base_vectors[ref_name+"_G"]+all_substitution_base_vectors[ref_name+"_T"]
    y_max = (max(max(tots),1))*1.1 #max to avoid ylim of 0,0 which makes python freak out
    ax.set_ylim(0,y_max) #max to avoid ylim of 0,0 which makes python freak out
    ax.set_xlim([0,ref_len])

    legend_patches = [pA[0],pC[0],pG[0],pT[0],pN[0]]
    legend_labels = ['A','C','T','G','N']

    if quantification_window_idxs and len(quantification_window_idxs) > 1:
        include_idxs_list = sorted(list(quantification_window_idxs))
        lastStart = include_idxs_list[0]
        lastIdx = include_idxs_list[0]
        for idx in range(1,len(include_idxs_list)):
            if include_idxs_list[idx] == lastIdx + 1:
                lastIdx = include_idxs_list[idx]
            else:
                p = matplotlib.patches.Rectangle((lastStart, 0), 1+(lastIdx-lastStart), y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                plt.gca().add_patch(p) #gca = get current axis
                lastStart = include_idxs_list[idx]
                lastIdx = include_idxs_list[idx]
        p = matplotlib.patches.Rectangle((lastStart, 0), 1+(lastIdx-lastStart), y_max,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
        plt.gca().add_patch(p)
        q_win_patch = patches.Patch(fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2,label='Quantification window')
        legend_patches.append(q_win_patch)
        legend_labels.append('Quantification window')

    ax.set_title(plot_title)
    ax.set_ylabel('% of total bases (Number of substitutions)')
    ax.set_xlabel('Reference position')
#            ax.set_facecolor('lightgray')
#    lgd=ax.legend((pA[0],pC[0],pG[0],pT[0],pN[0]),('A','C','G','T','N'),ncol=1)
#    lgd=ax.legend(legend_patches,legend_labels,ncol=1)
    lgd=ax.legend(handles=legend_patches, labels=legend_labels,loc='upper center', bbox_to_anchor=(0.3, -0.15),ncol=2, fancybox=True, shadow=True)

    y_label_values= np.round(np.linspace(0, min(max(tots),max(ax.get_yticks())),6))# np.arange(0,y_max,y_max/6.0)
    plt.yticks(y_label_values,['%.1f%% (%d)' % (n_reads/ref_count*100,n_reads) for n_reads in y_label_values])

    plt.savefig(fig_filename_root + '.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
    if save_also_png:
        plt.savefig(fig_filename_root + '.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close()

def plot_sub_freqs(alt_nuc_counts,plot_title,fig_filename_root,save_also_png):
    """
    Plots histogram of substitution frequencies for each nucleotide (from nuc X to nuc Y)
    input:
    alt_nuc_counts: dict of number of changes
    """

    #plot all substitution rates
    plt.figure(figsize=(8.3,8))

    alph = ['A','C','G','T','N']
    color_lookup = get_color_lookup(alph,alpha=1)

    plt.bar([1,2,3],[alt_nuc_counts['A']['C'],alt_nuc_counts['A']['G'],alt_nuc_counts['A']['T']],color=color_lookup['A'])
    plt.bar([5,6,7],[alt_nuc_counts['C']['A'],alt_nuc_counts['C']['G'],alt_nuc_counts['C']['T']],color=color_lookup['C'])
    plt.bar([9,10,11],[alt_nuc_counts['G']['A'],alt_nuc_counts['G']['C'],alt_nuc_counts['G']['T']],color=color_lookup['G'])
    plt.bar([13,14,15],[alt_nuc_counts['T']['A'],alt_nuc_counts['T']['C'],alt_nuc_counts['T']['G']],color=color_lookup['T'])
    plt.title(plot_title)
    plt.ylabel('Number of substitutions')
    plt.xticks([1,2,3,5,6,7,9,10,11,13,14,15],['A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G'],rotation='vertical')

    plt.savefig(fig_filename_root + '.pdf',bbox_inches='tight')
    if save_also_png:
        plt.savefig(fig_filename_root + '.png',bbox_inches='tight')
    plt.close()

def plot_nuc_freqs(df_nuc_freq,tot_aln_reads,plot_title,fig_filename_root,save_also_png):
    """
    Plots a heatmap of the percentage of reads that had each nucletide at each base in the reference
    Positions in the reference that have more than one allele can be spotted using this plot
    """

    fig, ax = plt.subplots(figsize=(len(df_nuc_freq.columns),4))
#                sns_plot = sns.heatmap(df_nuc_freq,vmin=0,vmax=tot_aln_reads,cmap="YlGnBu",square=True).get_figure()
    sns.heatmap(df_nuc_freq,vmin=0,vmax=tot_aln_reads,cmap="YlGnBu",square=True,ax=ax)
    ax.set_title(plot_title)

    plt.savefig(fig_filename_root + '.pdf',bbox_inches='tight')
    if save_also_png:
        plt.savefig(fig_filename_root + '.png',bbox_inches='tight')
    plt.close()

def plot_log_nuc_freqs(df_nuc_freq,tot_aln_reads,plot_title,fig_filename_root,save_also_png,quantification_window_idxs=None):
    """
    Plots a heatmap of the percentage of reads that had each nucletide at each base in the reference
    Positions in the reference that have more than one allele can be spotted using this plot
    """

    fig, ax = plt.subplots(figsize=(len(df_nuc_freq.columns),4))
#                sns_plot = sns.heatmap(df_nuc_freq,vmin=0,vmax=tot_aln_reads,cmap="YlGnBu",square=True).get_figure()
    sns.heatmap(np.log2(df_nuc_freq+1),vmin=0,vmax=np.log2(tot_aln_reads+1),cmap="YlGnBu",ax=ax)#,xticklabels=1)
    ax.set_title(plot_title)

    if quantification_window_idxs is not None:
        q_win_y_start = -0.1
        q_win_y_height = 6.1

        q_list = sorted(list(quantification_window_idxs))

        lastStart = q_list[0]
        lastIdx = q_list[0]
        for idx in range(1,len(q_list)):
            if q_list[idx] == lastIdx + 1:
                lastIdx = q_list[idx]
            else:
                ax.add_patch(
                    patches.Rectangle((2+lastStart, q_win_y_start), 1+(lastIdx-lastStart), q_win_y_height,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
                    )
                lastStart = q_list[idx]
                lastIdx = q_list[idx]
        ax.add_patch(
            patches.Rectangle((2+lastStart, q_win_y_start), 1+(lastIdx-lastStart), q_win_y_height,fill=None,edgecolor=(0,0,0,0.25),linestyle=(0,(5,2)),linewidth=2)
            )

    plt.savefig(fig_filename_root + '.pdf',bbox_inches='tight')
    if save_also_png:
        plt.savefig(fig_filename_root + '.png',bbox_inches='tight')
    plt.close()


def plot_conversion_at_sel_nucs(df_subs,ref_name,ref_sequence,plot_title,conversion_nuc_from,fig_filename_root,save_also_png):
    '''
    Plots the conversion at selected nucleotides
    Looks for the 'conversion_nuc_from' in the ref_sequence and sets those as 'selected nucleotides'
    At selected nucleotides, the proportion of each base is shown as a barplot
    '''
    nucs = list(df_subs.index)
    color_lookup = get_color_lookup(nucs,alpha=1)
    amp_len = len(ref_sequence)

    fig = plt.figure(figsize=(amp_len,6))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])

    from_nuc_indices = [pos for pos, char in enumerate (ref_sequence) if char == conversion_nuc_from]
    ax = plt.subplot(gs[0])
    ind = np.arange(len(from_nuc_indices))

    #pandas was the cause of much blood, sweat, and tears which were shed here
    #bottom_so_far keeps track of the bottom of each barplot for each base
    bottom_so_far = np.zeros(len(from_nuc_indices))
    for n in nucs:
        vals = df_subs.iloc[:,from_nuc_indices].loc[n] * 100
        ax.bar(ind,vals,color=color_lookup[n],bottom=bottom_so_far)
        bottom_so_far += vals

    ax.set_xticks(ind)
    ax.set_xticklabels([conversion_nuc_from + str(x+1) for x in from_nuc_indices])
    ax.set_xlim([-0.5,len(ind)-0.5])
    ax.set_ylabel('Nucleotide frequency',fontsize='small')

    #plot legend
    legend_patches = []
    for nuc in nucs:
        patch = patches.Patch(color=color_lookup[nuc], label=nuc)
        legend_patches.append(patch)

    plt.legend(handles=legend_patches,loc='center left',ncol=1,bbox_to_anchor=(1,0.5))

    ax = plt.subplot(gs[1])
    #draw reference sequence
    ref_y_start = 0
    ref_y_height = 1
    for pos_ind in range(amp_len):
        ax.add_patch(
            patches.Rectangle((pos_ind, ref_y_start), 1, ref_y_height,facecolor=color_lookup[ref_sequence[pos_ind]],edgecolor='w')
            )
        ax.text(pos_ind+0.5,ref_y_start + ref_y_height/2.0,ref_sequence[pos_ind],horizontalalignment='center',verticalalignment='center')

    ax.set_xticks([x + 0.5 for x in from_nuc_indices])
    ax.set_xticklabels([conversion_nuc_from + str(x+1) for x in from_nuc_indices])
    ax.set_xlim([0,amp_len])

    ax.set_yticks([0.5])
    ax.set_yticklabels(['Reference'])


    plt.tight_layout()
    fig.savefig(fig_filename_root+'.pdf',bbox_inches='tight')
    if save_also_png:
        fig.savefig(fig_filename_root+'.png',bbox_inches='tight',pad=1)
    plt.close()

def plot_conversion_at_sel_nucs_not_include_ref(df_subs,ref_name,ref_sequence,plot_title,conversion_nuc_from,fig_filename_root,save_also_png):
    '''
    Plots the conversion at selected nucleotides but ignores non-substitutions (for example at nucs that are 'C' in the reference, bars show the proportion of A T G (not C))
    Looks for the 'conversion_nuc_from' in the ref_sequence and sets those as 'selected nucleotides'
    At selected nucleotides, the proportion of each substitution is shown as a barplot
    '''
    nucs = list(df_subs.index)
    color_lookup = get_color_lookup(nucs,alpha=1)
    amp_len = len(ref_sequence)

    fig = plt.figure(figsize=(amp_len,6))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])

    from_nuc_indices = [pos for pos, char in enumerate (ref_sequence) if char == conversion_nuc_from]
    ax = plt.subplot(gs[0])
    ind = np.arange(len(from_nuc_indices))


    nucs_only_sub = nucs
    nucs_only_sub.remove(conversion_nuc_from)
    sub_freq = df_subs.iloc[:,from_nuc_indices].loc[nucs_only_sub].sum(axis=0)
    bottom_so_far = np.zeros(len(from_nuc_indices))
    for n in nucs_only_sub:
        vals = df_subs.iloc[:,from_nuc_indices].loc[n]
        pcts = (vals/sub_freq.mask(sub_freq == 0)).fillna(0) * 100
        ax.bar(ind,pcts,color=color_lookup[n],bottom=bottom_so_far)
        bottom_so_far += pcts

    #bottom x axis
    ax.set_xticks(ind)
    ax.set_xticklabels([conversion_nuc_from + str(x+1) for x in from_nuc_indices])
    ax.set_xlim([-0.5,len(ind)-0.5])
    ax.set_ylabel('Percentage non-reference',fontsize='small')

    #top x axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.append(ind,len(from_nuc_indices)-0.5))
    ax2.set_xticklabels(["%0.2f%%"%(x*100)  for x in sub_freq] +["         $\it{\% Non-" + conversion_nuc_from + "}$"])

    #plot legend
    legend_patches = []
    for nuc in nucs_only_sub:
        patch = patches.Patch(color=color_lookup[nuc], label=nuc)
        legend_patches.append(patch)

    plt.legend(handles=legend_patches,loc='center left',ncol=1,bbox_to_anchor=(1,0.5))

    ax = plt.subplot(gs[1])
    #draw reference sequence
    ref_y_start = 0
    ref_y_height = 1
    for pos_ind in range(amp_len):
        ax.add_patch(
            patches.Rectangle((pos_ind, ref_y_start), 1, ref_y_height,facecolor=color_lookup[ref_sequence[pos_ind]],edgecolor='w')
            )
        ax.text(pos_ind+0.5,ref_y_start + ref_y_height/2.0,ref_sequence[pos_ind],horizontalalignment='center',verticalalignment='center')

    ax.set_xticks([x + 0.5 for x in from_nuc_indices])
    ax.set_xticklabels([conversion_nuc_from + str(x+1) for x in from_nuc_indices])
    ax.set_xlim([0,amp_len])

    ax.set_yticks([0.5])
    ax.set_yticklabels(['Reference'])


    plt.tight_layout()
    fig.savefig(fig_filename_root+'.pdf',bbox_inches='tight')
    if save_also_png:
        fig.savefig(fig_filename_root+'.png',bbox_inches='tight',pad=1)
    plt.close()

def plot_conversion_at_sel_nucs_not_include_ref_scaled(df_subs,ref_name,ref_sequence,plot_title,conversion_nuc_from,fig_filename_root,save_also_png):
    '''
    Plots the conversion at selected nucleotides not including reference base, scaled by number of events
    Looks for the 'conversion_nuc_from' in the ref_sequence and sets those as 'selected nucleotides'
    At selected nucleotides, the count of each base is shown as a barplot
    '''
    nucs = list(df_subs.index)
    color_lookup = get_color_lookup(nucs,alpha=1)
    nucs.remove(conversion_nuc_from)
    amp_len = len(ref_sequence)

    fig = plt.figure(figsize=(amp_len,6))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])

    from_nuc_indices = [pos for pos, char in enumerate (ref_sequence) if char == conversion_nuc_from]
    ax = plt.subplot(gs[0])
    ind = np.arange(len(from_nuc_indices))

    #pandas was the cause of much blood, sweat, and tears which were shed here
    #bottom_so_far keeps track of the bottom of each barplot for each base
    bottom_so_far = np.zeros(len(from_nuc_indices))
    for n in nucs:
        vals = df_subs.iloc[:,from_nuc_indices].loc[n] * 100
        ax.bar(ind,vals,color=color_lookup[n],bottom=bottom_so_far)
        bottom_so_far += vals

    ax.set_xticks(ind)
    ax.set_xticklabels([conversion_nuc_from + str(x+1) for x in from_nuc_indices])
    ax.set_xlim([-0.5,len(ind)-0.5])
    ax.set_ylabel('Nucleotide frequency',fontsize='small')

    #plot legend
    legend_patches = []
    for nuc in nucs:
        patch = patches.Patch(color=color_lookup[nuc], label=nuc)
        legend_patches.append(patch)

    plt.legend(handles=legend_patches,loc='center left',ncol=1,bbox_to_anchor=(1,0.5))

    ax = plt.subplot(gs[1])
    #draw reference sequence
    ref_y_start = 0
    ref_y_height = 1
    for pos_ind in range(amp_len):
        ax.add_patch(
            patches.Rectangle((pos_ind, ref_y_start), 1, ref_y_height,facecolor=color_lookup[ref_sequence[pos_ind]],edgecolor='w')
            )
        ax.text(pos_ind+0.5,ref_y_start + ref_y_height/2.0,ref_sequence[pos_ind],horizontalalignment='center',verticalalignment='center')

    ax.set_xticks([x + 0.5 for x in from_nuc_indices])
    ax.set_xticklabels([conversion_nuc_from + str(x+1) for x in from_nuc_indices])
    ax.set_xlim([0,amp_len])

    ax.set_yticks([0.5])
    ax.set_yticklabels(['Reference'])


    plt.tight_layout()
    fig.savefig(fig_filename_root+'.pdf',bbox_inches='tight')
    if save_also_png:
        fig.savefig(fig_filename_root+'.png',bbox_inches='tight',pad=1)
    plt.close()

### Allele plot
#We need to customize the seaborn heatmap class and function
class Custom_HeatMapper(sns.matrix._HeatMapper):

    def __init__(self, data, vmin, vmax, cmap, center, robust, annot, fmt,
                 annot_kws,per_element_annot_kws,cbar, cbar_kws,
                 xticklabels=True, yticklabels=True, mask=None):

        super(Custom_HeatMapper, self).__init__(data, vmin, vmax, cmap, center, robust, annot, fmt,
                 annot_kws, cbar, cbar_kws,
                 xticklabels, yticklabels, mask)


        if annot is not None:
            if per_element_annot_kws is None:
                self.per_element_annot_kws=np.empty_like(annot,dtype=np.object)
                self.per_element_annot_kws[:]=dict()
            else:
                self.per_element_annot_kws=per_element_annot_kws

    #add per element dict to syle the annotatiin
    def _annotate_heatmap(self, ax, mesh):
        """Add textual labels with the value in each cell."""
        mesh.update_scalarmappable()
        xpos, ypos = np.meshgrid(ax.get_xticks(), ax.get_yticks())


        for x, y, m, color, val,per_element_dict  in zip(xpos.flat, ypos.flat,
                                       mesh.get_array(), mesh.get_facecolors(),
                                       self.annot_data.flat,self.per_element_annot_kws.flat):
            #print per_element_dict
            if m is not np.ma.masked:
                l = sns.utils.relative_luminance(color)
                text_color = ".15" if l > .408 else "w"
                annotation = ("{:" + self.fmt + "}").format(val)
                text_kwargs = dict(color=text_color, ha="center", va="center")
                text_kwargs.update(self.annot_kws)
                text_kwargs.update(per_element_dict)

                ax.text(x, y, annotation, **text_kwargs)


    #removed the colobar
    def plot(self, ax, cax, kws):
        """Draw the heatmap on the provided Axes."""
        # Remove all the Axes spines
        sns.utils.despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        mesh = ax.pcolormesh(self.plot_data, vmin=self.vmin, vmax=self.vmax,
                             cmap=self.cmap, **kws)

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Add row and column labels
        ax.set(xticks=self.xticks, yticks=self.yticks)
        xtl = ax.set_xticklabels(self.xticklabels)
        ytl = ax.set_yticklabels(self.yticklabels, rotation="vertical")

        # Possibly rotate them if they overlap
        plt.draw()
        if sns.utils.axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if sns.utils.axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, mesh)

def custom_heatmap(data, vmin=None, vmax=None, cmap=None, center=None, robust=False,
            annot=None, fmt=".2g", annot_kws=None,per_element_annot_kws=None,
            linewidths=0, linecolor="white",
            cbar=True, cbar_kws=None, cbar_ax=None,
            square=False, ax=None, xticklabels=True, yticklabels=True,
            mask=None,
            **kwargs):

    # Initialize the plotter object
    plotter = Custom_HeatMapper(data, vmin, vmax, cmap, center, robust, annot, fmt,
                          annot_kws, per_element_annot_kws,cbar, cbar_kws, xticklabels,
                          yticklabels, mask)

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax

def prep_alleles_table(df_alleles,MAX_N_ROWS,MIN_FREQUENCY):
    """
    Prepares a df of alleles for Plotting
    input:
    -df_alleles: pandas dataframe of alleles to plot
    -MAX_N_ROWS: max number of rows to plot
    -MIN_FREQUENCY: min frequency for a row to be plotted
    returns:
    -X: list of numbers representing nucleotides of the allele
    -annot: list of nucleotides (letters) of the allele
    -y_labels: list of labels for each row/allele
    -insertion_dict: locations of insertions -- red squares will be drawn around these
    -per_element_annot_kws: annotations for each cell (e.g. bold for substitutions, etc.)
    """
    dna_to_numbers={'-':0,'A':1,'T':2,'C':3,'G':4,'N':5}
    seq_to_numbers= lambda seq: [dna_to_numbers[x] for x in seq]

    X=[]
    annot=[]
    y_labels=[]
    insertion_dict=defaultdict(list)
    per_element_annot_kws=[]

    re_find_indels=re.compile("(-*-)")
    idx_row=0
    for idx,row in df_alleles.ix[df_alleles['%Reads']>=MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(seq_to_numbers(str.upper(idx)))
        annot.append(list(idx))
        y_labels.append('%.2f%% (%d reads)' % (row['%Reads'],row['#Reads']))


        for p in re_find_indels.finditer(row['Reference_Sequence']):
            insertion_dict[idx_row].append((p.start(),p.end()))

        idx_row+=1


        idxs_sub= [i_sub for i_sub in range(len(idx)) if \
                   (row['Reference_Sequence'][i_sub]!=idx[i_sub]) and \
                   (row['Reference_Sequence'][i_sub]!='-') and\
                   (idx[i_sub]!='-')]
        to_append=np.array([{}]*len(idx),dtype=np.object)
        to_append[ idxs_sub]={'weight':'bold', 'color':'black','size':16}
        per_element_annot_kws.append(to_append)

    return X,annot,y_labels,insertion_dict,per_element_annot_kws

def prep_alleles_table_compare(df_alleles,sample_name_1,sample_name_2,MAX_N_ROWS,MIN_FREQUENCY):
    """
    Prepares a df of alleles for Plotting
    takes a merged allele table, and sets labels to read percents and counts from each sample
    input:
    -df_alleles: merged pandas dataframe of alleles to plot
    -sample_name_1: sample name 1
    -sample_name_2: sample name 2
    --- y_labels will be determined using the columns named like: '#Reads_s1' (where s1 is the sample 1 name)
    -MAX_N_ROWS: max number of rows to plot
    -MIN_FREQUENCY: min frequency for a row to be plotted
    returns:
    -X: list of numbers representing nucleotides of the allele
    -annot: list of nucleotides (letters) of the allele
    -y_labels: list of labels for each row/allele
    -insertion_dict: locations of insertions -- red squares will be drawn around these
    -per_element_annot_kws: annotations for each cell (e.g. bold for substitutions, etc.)
    """
    dna_to_numbers={'-':0,'A':1,'T':2,'C':3,'G':4,'N':5}
    seq_to_numbers= lambda seq: [dna_to_numbers[x] for x in seq]

    X=[]
    annot=[]
    y_labels=[]
    insertion_dict=defaultdict(list)
    per_element_annot_kws=[]

    re_find_indels=re.compile("(-*-)")
    idx_row=0
    for idx,row in df_alleles.ix[df_alleles['%Reads_'+sample_name_1] + df_alleles['%Reads_'+sample_name_2]>=MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(seq_to_numbers(str.upper(idx)))
        annot.append(list(idx))
        y_labels.append('%.2f%% (%d reads) %.2f%% (%d reads) ' % (row['%Reads_'+sample_name_1],row['#Reads_'+sample_name_1],
                                                    row['%Reads_'+sample_name_2],row['#Reads_'+sample_name_2]))


        for p in re_find_indels.finditer(row['Reference_Sequence']):
            insertion_dict[idx_row].append((p.start(),p.end()))

        idx_row+=1


        idxs_sub= [i_sub for i_sub in range(len(idx)) if \
                   (row['Reference_Sequence'][i_sub]!=idx[i_sub]) and \
                   (row['Reference_Sequence'][i_sub]!='-') and\
                   (idx[i_sub]!='-')]
        to_append=np.array([{}]*len(idx),dtype=np.object)
        to_append[ idxs_sub]={'weight':'bold', 'color':'black','size':16}
        per_element_annot_kws.append(to_append)

    return X,annot,y_labels,insertion_dict,per_element_annot_kws

def plot_alleles_heatmap(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,SAVE_ALSO_PNG=False,base_editor_output=False,sgRNA_intervals=None):
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
    -base_editor_output: if true, won't draw 'predicted cleavage' line
    """
    plot_nuc_len=len(reference_seq)

    # make a color map of fixed colors
    alpha=0.4
    A_color=get_nuc_color('A',alpha)
    T_color=get_nuc_color('T',alpha)
    C_color=get_nuc_color('C',alpha)
    G_color=get_nuc_color('G',alpha)
    INDEL_color = get_nuc_color('N',alpha)

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

    fig=plt.figure(figsize=(plot_nuc_len*0.3,(N_ROWS+1)*0.6))
    gs1 = gridspec.GridSpec(N_ROWS+1,N_COLUMNS)
    gs2 = gridspec.GridSpec(N_ROWS+1,N_COLUMNS)

    #ax_hm_ref heatmap for the reference
    ax_hm_ref=plt.subplot(gs1[0, :])
    ax_hm=plt.subplot(gs2[1:, :])

    custom_heatmap(ref_seq_hm,annot=ref_seq_annot_hm,annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm_ref,vmin=0,vmax=5,square=True)
    custom_heatmap(X,annot=np.array(annot),annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm,vmin=0,vmax=5,square=True, per_element_annot_kws=per_element_annot_kws)

    ax_hm.yaxis.tick_right()
    ax_hm.yaxis.set_ticklabels(y_labels[::-1],rotation=True),
    ax_hm.xaxis.set_ticks([])

# todo -- add sgRNAs below reference plot
#    if sgRNA_intervals:
#        sgRNA_y_start = 0.3
#        sgRNA_y_height = 0.1
#        min_sgRNA_x = None
#        for idx,sgRNA_int in enumerate(sgRNA_intervals):
#            ax_hm_ref.add_patch(
#                patches.Rectangle((2+sgRNA_int[0], sgRNA_y_start), 1+sgRNA_int[1]-sgRNA_int[0], sgRNA_y_height,facecolor=(0,0,0,0.15))
#                )
#            #set left-most sgrna start
#            if not min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#            if sgRNA_int[0] < min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#        ax_hm_ref.text(2+min_sgRNA_x,sgRNA_y_start + sgRNA_y_height/2,'sgRNA ',horizontalalignment='right',verticalalignment='center')

    #print lines


    #create boxes for ins
    for idx,lss in insertion_dict.iteritems():
        for ls in lss:
            for l in ls:
                ax_hm.vlines([l],N_ROWS-idx-1,N_ROWS-idx,color='red',linewidths=3)

            ax_hm.hlines(N_ROWS-idx-1,ls[0],ls[1],color='red',linewidths=3)
            ax_hm.hlines(N_ROWS-idx,ls[0],ls[1],color='red',linewidths=3)

    #cut point vertical line
    if base_editor_output is False:
        ax_hm.vlines([plot_nuc_len/2],*ax_hm.get_ylim(),linestyles='dashed')


    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(['Reference'],rotation=True)



    gs2.update(left=0,right=1, hspace=0.05,wspace=0,top=1*(((N_ROWS)*1.13))/(N_ROWS))
    gs1.update(left=0,right=1, hspace=0.05,wspace=0,)

    sns.set_context(rc={'axes.facecolor':'white','lines.markeredgewidth': 1,'mathtext.fontset' : 'stix','text.usetex':True,'text.latex.unicode':True} )

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                    mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'),ms=18),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='red', marker='s',ms=8,markeredgewidth=2.5),
              matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='black', marker='_',ms=2,)]
    descriptions=['Substitutions','Insertions','Deletions']

    if base_editor_output is False:
        proxies.append(
              matplotlib.lines.Line2D([0], [1], linestyle='--',c='black',ms=6))
        descriptions.append('Predicted cleavage position')

    #ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)
    lgd = ax_hm.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='upper center', bbox_to_anchor=(0.5, 0),ncol=1,fancybox=True,shadow=False)

    plt.savefig(fig_filename_root+'.pdf',bbox_inches='tight',bbox_extra_artists=(lgd,))
    if SAVE_ALSO_PNG:
        plt.savefig(fig_filename_root+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,))
    plt.close()

def plot_alleles_heatmap_hist(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,count_values,SAVE_ALSO_PNG=False,base_editor_output=False,sgRNA_intervals=None):
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
    -base_editor_output: if true, won't draw 'predicted cleavage' line
    """
    plot_nuc_len=len(reference_seq)

    # make a color map of fixed colors
    alpha=0.4
    A_color=get_nuc_color('A',alpha)
    T_color=get_nuc_color('T',alpha)
    C_color=get_nuc_color('C',alpha)
    G_color=get_nuc_color('G',alpha)
    INDEL_color = get_nuc_color('N',alpha)

    dna_to_numbers={'-':0,'A':1,'T':2,'C':3,'G':4,'N':5}
    seq_to_numbers= lambda seq: [dna_to_numbers[x] for x in seq]

    cmap = colors_mpl.ListedColormap([INDEL_color, A_color,T_color,C_color,G_color,INDEL_color])

    #ref_seq_around_cut=reference_seq[max(0,cut_point-plot_nuc_len/2+1):min(len(reference_seq),cut_point+plot_nuc_len/2+1)]

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

    fig=plt.figure(figsize=(plot_nuc_len*0.3,(N_ROWS+1)*0.6))
    gs = gridspec.GridSpec(2,2,width_ratios=[1,N_ROWS],height_ratios=[N_COLUMNS,1])

    #ax_hm_ref heatmap for the reference
    ax_hm_ref = plt.subplot(gs[0])
    ax_hm = plt.subplot(gs[2])
    ax_bar = plt.subplot(gs[3])

    custom_heatmap(ref_seq_hm,annot=ref_seq_annot_hm,annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm_ref,vmin=0,vmax=5,square=True)
    custom_heatmap(X,annot=np.array(annot),annot_kws={'size':16},cmap=cmap,fmt='s',ax=ax_hm,vmin=0,vmax=5,square=True, per_element_annot_kws=per_element_annot_kws)

    ax_hm.yaxis.tick_right()
    ax_hm.yaxis.set_ticklabels(y_labels[::-1],rotation=True),
    ax_hm.xaxis.set_ticks([])

# todo -- add sgRNAs below reference plot
#    if sgRNA_intervals:
#        sgRNA_y_start = 0.3
#        sgRNA_y_height = 0.1
#        min_sgRNA_x = None
#        for idx,sgRNA_int in enumerate(sgRNA_intervals):
#            ax_hm_ref.add_patch(
#                patches.Rectangle((2+sgRNA_int[0], sgRNA_y_start), 1+sgRNA_int[1]-sgRNA_int[0], sgRNA_y_height,facecolor=(0,0,0,0.15))
#                )
#            #set left-most sgrna start
#            if not min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#            if sgRNA_int[0] < min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#        ax_hm_ref.text(2+min_sgRNA_x,sgRNA_y_start + sgRNA_y_height/2,'sgRNA ',horizontalalignment='right',verticalalignment='center')

    #print lines

    #cut point vertical line
    if base_editor_output is False:
        ax_hm.vlines([plot_nuc_len/2],*ax_hm.get_ylim(),linestyles='dashed')

    #create boxes for ins
    for idx,lss in insertion_dict.iteritems():
        for ls in lss:
            for l in ls:
                ax_hm.vlines([l],N_ROWS-idx-1,N_ROWS-idx,color='red',linewidths=3)

            ax_hm.hlines(N_ROWS-idx-1,ls[0],ls[1],color='red',linewidths=3)
            ax_hm.hlines(N_ROWS-idx,ls[0],ls[1],color='red',linewidths=3)


    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(['Reference'],rotation=True)



#    gs2.update(left=0,right=1, hspace=0.05,wspace=0,top=1*(((N_ROWS)*1.13))/(N_ROWS))
#    gs1.update(left=0,right=1, hspace=0.05,wspace=0,)

    sns.set_context(rc={'axes.facecolor':'white','lines.markeredgewidth': 1,'mathtext.fontset' : 'stix','text.usetex':True,'text.latex.unicode':True} )

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                    mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'),ms=18),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='red', marker='s',ms=8,markeredgewidth=2.5),
              matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='black', marker='_',ms=2,)]
    descriptions=['Substitutions','Insertions','Deletions']

    if base_editor_output is False:
        proxies.append(
              matplotlib.lines.Line2D([0], [1], linestyle='--',c='black',ms=6))
        descriptions.append('Predicted cleavage position')

    #ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)
    lgd = ax_hm.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='upper center', bbox_to_anchor=(0.5, 0),ncol=1,fancybox=True,shadow=False)

    plt.savefig(fig_filename_root+'.pdf',bbox_inches='tight',bbox_extra_artists=(lgd,))
    if SAVE_ALSO_PNG:
        plt.savefig(fig_filename_root+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,),pad=1)
    plt.close()

def plot_alleles_table(reference_seq,df_alleles,fig_filename_root,MIN_FREQUENCY=0.5,MAX_N_ROWS=100,SAVE_ALSO_PNG=False,base_editor_output=False,sgRNA_intervals=None):
    """
    plots an allele table for a dataframe with allele frequencies
    input:
    reference_seq: the reference amplicon sequence to plot
    df_alleles: merged dataframe (should include columns "#Reads','%Reads')
    fig_filename: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    sgRNA_intervals: locations where sgRNA is located
    """
    X,annot,y_labels,insertion_dict,per_element_annot_kws = prep_alleles_table(df_alleles,MAX_N_ROWS,MIN_FREQUENCY)
    plot_alleles_heatmap(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,SAVE_ALSO_PNG,base_editor_output,sgRNA_intervals)

def plot_alleles_table_from_file(alleles_file_name,fig_filename_root,MIN_FREQUENCY=0.5,MAX_N_ROWS=100,SAVE_ALSO_PNG=False,base_editor_output=False,sgRNA_intervals=None):
    """
    plots an allele table for a dataframe with allele frequencies
    infers the reference sequence by finding reference sequences without gaps (-)
    This function is only used for one-off plotting purposes and not for the general CRISPResso analysis because it dies if the ref can't be found

    input:
    alleles_file_name: file name with alleles (should include columns "#Reads','%Reads')
    fig_filename: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    sgRNA_intervals: locations where sgRNA is located
    """

    df_alleles = pd.read_table(alleles_file_name)
    df_alleles = df_alleles.reset_index().set_index('Aligned_Sequence')

    rows_include_reference_seq = df_alleles.loc[df_alleles['Reference_Sequence'].str.contains('-')==False]
    if len(rows_include_reference_seq) > 0:
        reference_seq = rows_include_reference_seq['Reference_Sequence'].iloc[0]
    else:
        raise Exception('Could not infer reference sequence from allele table')

    X,annot,y_labels,insertion_dict,per_element_annot_kws = prep_alleles_table(df_alleles,MAX_N_ROWS,MIN_FREQUENCY)
    plot_alleles_heatmap(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,SAVE_ALSO_PNG,base_editor_output,sgRNA_intervals)

def plot_alleles_table_compare(reference_seq,df_alleles,sample_name_1,sample_name_2,fig_filename_root,MIN_FREQUENCY=0.5,MAX_N_ROWS=100,SAVE_ALSO_PNG=False,base_editor_output=False,sgRNA_intervals=None):
    """
    plots an allele table for a dataframe with allele frequencies from two CRISPResso runs
    input:
    reference_seq: the reference amplicon sequence to plot
    df_alleles: merged dataframe (should include columns "#Reads_s1','%Reads_s1','#Reads_s2','%Reads_s2','each_LFC'" where s1 and s2 are the sample names)
    sample_name_1, sample_name_2: names of the samples
    fig_filename: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    """
    X,annot,y_labels,insertion_dict,per_element_annot_kws = prep_alleles_table_compare(df_alleles,sample_name_1,sample_name_2,MAX_N_ROWS,MIN_FREQUENCY)
    plot_alleles_heatmap(reference_seq,fig_filename_root,X,annot,y_labels,insertion_dict,per_element_annot_kws,SAVE_ALSO_PNG,base_editor_output,sgRNA_intervals)
