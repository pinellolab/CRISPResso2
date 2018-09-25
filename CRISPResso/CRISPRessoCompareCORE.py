# -*- coding: utf-8 -*-
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''
import os
import errno
import sys
import traceback
import argparse
import re
import cPickle as cp
from CRISPResso import CRISPRessoShared
from CRISPResso import CRISPRessoPlot



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
                error('You need to install %s module to use CRISPRessoCompare!' % library_name)
                sys.exit(1)


def get_amplicon_output(amplicon_name,output_folder):
    profile_file=os.path.join(output_folder,amplicon_name+'.effect_vector_combined.txt')
    if os.path.exists(quantification_file) and profile_file:
        return quantification_file,profile_file
    else:
        raise CRISPRessoShared.OutputFolderIncompleteException('The folder %s  is not a valid CRISPResso2 output folder. Cannot find profile file %s for amplicon %s.' % (output_folder,profile_file,amplicon_name))

def parse_profile(profile_file):
    return np.loadtxt(profile_file,skiprows=1)


def load_cut_points_sgRNA_intervals(output_folder,amplicon_name):
    cut_points_file=os.path.join(output_folder,amplicon_name+'.cut_points.pickle')
    sgRNA_intervals_file=os.path.join(output_folder,amplicon_name+'.sgRNA_intervals.pickle')
    if os.path.exists(cut_points_file):
        cut_points=cp.load(open(cut_points_file))
    else:
        cut_points=[]


    if os.path.exists(sgRNA_intervals_file):
        sgRNA_intervals=cp.load(open(sgRNA_intervals_file))
    else:
        sgRNA_intervals=[]

    return  cut_points,sgRNA_intervals


###EXCEPTIONS############################

class MixedRunningModeException(Exception):
    pass

class DifferentAmpliconLengthException(Exception):
    pass
############################



matplotlib=check_library('matplotlib')
from matplotlib import font_manager as fm
font = {'size'   : 20}
matplotlib.rc('font', **font)
matplotlib.use('Agg')

plt=check_library('pylab')
np=check_library('numpy')
pd=check_library('pandas')
#scipy=check_library('scipy.stats')
import scipy.stats as stats


_ROOT = os.path.abspath(os.path.dirname(__file__))



def main():
    try:
        description = ['~~~CRISPRessoCompare~~~','-Comparison of two CRISPResso analyses-']
        compare_header = r'''
 ___________________________
| __ __      __      __  __ |
|/  /  \|\/||__) /\ |__)|_  |
|\__\__/|  ||   /--\| \ |__ |
|___________________________|
        '''
        compare_header = CRISPRessoShared.get_crispresso_header(description,compare_header)
        print(compare_header)

        parser = argparse.ArgumentParser(description='CRISPRessoCompare Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('crispresso_output_folder_1', type=str,  help='First output folder with CRISPResso analysis')
        parser.add_argument('crispresso_output_folder_2', type=str,  help='Second output folder with CRISPResso analysis')

        #OPTIONALS
        parser.add_argument('-n','--name',  help='Output name', default='')
        parser.add_argument('-n1','--sample_1_name',  help='Sample 1 name', default='Sample_1')
        parser.add_argument('-n2','--sample_2_name',  help='Sample 2 name', default='Sample_2')
        parser.add_argument('-o','--output_folder',  help='', default='')
        parser.add_argument('--min_frequency_alleles_around_cut_to_plot', type=float, help='Minimum %% reads required to report an allele in the alleles table plot.', default=0.2)
        parser.add_argument('--max_rows_alleles_around_cut_to_plot',  type=int, help='Maximum number of rows to report in the alleles table plot. ', default=50)
        parser.add_argument('--save_also_png',help='Save also .png images additionally to .pdf files',action='store_true')
        parser.add_argument('--debug', help='Show debug messages', action='store_true')

        args = parser.parse_args()
        debug_flag = args.debug


        #check that the CRISPResso output is present and fill amplicon_info
        quantification_file_1,amplicon_names_1,amplicon_info_1=CRISPRessoShared.check_output_folder(args.crispresso_output_folder_1)
        quantification_file_2,amplicon_names_2,amplicon_info_2=CRISPRessoShared.check_output_folder(args.crispresso_output_folder_2)

        get_name_from_folder=lambda x: os.path.basename(os.path.abspath(x)).replace('CRISPResso_on_','')

        if not args.name:
                 database_id='%s_VS_%s' % (get_name_from_folder(args.crispresso_output_folder_1),get_name_from_folder(args.crispresso_output_folder_2))
        else:
                 database_id=args.name


        OUTPUT_DIRECTORY='CRISPRessoCompare_on_%s' % database_id

        if args.output_folder:
                 OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)

        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoCompare_RUNNING_LOG.txt')


        try:
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 os.makedirs(OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)

        log_filename=_jp('CRISPRessoCompare_RUNNING_LOG.txt')
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        with open(log_filename,'w+') as outfile:
                  outfile.write('[Command used]:\nCRISPRessoCompare %s\n\n[Execution log]:\n' % ' '.join(sys.argv))

        #LOAD DATA
        amplicon_names_in_both = [amplicon_name for amplicon_name in amplicon_names_1 if amplicon_name in amplicon_names_2]
        n_refs = len(amplicon_names_in_both)
        def get_plot_title_with_ref_name(plotTitle,refName):
            if n_refs > 1:
                return (plotTitle + ": " + refName)
            return plotTitle

        for amplicon_name in amplicon_names_in_both:
            profile_1=parse_profile(amplicon_info_1[amplicon_name]['quantification_file'])
            profile_2=parse_profile(amplicon_info_2[amplicon_name]['quantification_file'])

            try:
                assert np.all(profile_1[:,0]==profile_2[:,0])
            except:
                raise DifferentAmpliconLengthException('Different amplicon lengths for the two amplicons.')
            len_amplicon=profile_1.shape[0]
            effect_vector_any_1=profile_1[:,1]
            effect_vector_any_2=profile_2[:,1]
            cut_points,sgRNA_intervals=load_cut_points_sgRNA_intervals(args.crispresso_output_folder_1,amplicon_name)


            #Quantification comparison barchart
            fig=plt.figure(figsize=(30,15))
            n_groups = 2

            N_TOTAL_1 = float(amplicon_info_1[amplicon_name]['Total'])
            N_UNMODIFIED_1 = float(amplicon_info_1[amplicon_name]['Unmodified'])
            N_MODIFIED_1 = float(amplicon_info_1[amplicon_name]['Modified'])

            N_TOTAL_2 = float(amplicon_info_2[amplicon_name]['Total'])
            N_UNMODIFIED_2 = float(amplicon_info_2[amplicon_name]['Unmodified'])
            N_MODIFIED_2 = float(amplicon_info_2[amplicon_name]['Modified'])


            means_sample_1= np.array([N_UNMODIFIED_1,N_MODIFIED_1])/N_TOTAL_1*100
            means_sample_2 = np.array([N_UNMODIFIED_2,N_MODIFIED_2])/N_TOTAL_2*100

            ax1=fig.add_subplot(1,2,1)

            index = np.arange(n_groups)
            bar_width = 0.35

            opacity = 0.4
            error_config = {'ecolor': '0.3'}

            rects1 = ax1.bar(index, means_sample_1, bar_width,
                             alpha=opacity,
                             color=(0,0,1,0.4),
                             label=args.sample_1_name)

            rects2 = ax1.bar(index + bar_width, means_sample_2, bar_width,
                             alpha=opacity,
                             color=(1,0,0,0.4),
                             label=args.sample_2_name)

            plt.ylabel('% Sequences')
            plt.title(get_plot_title_with_ref_name('%s VS %s' % (args.sample_1_name,args.sample_2_name),amplicon_name))
            plt.xticks(index + bar_width/2.0, ('Unmodified', 'Modified'))
            plt.legend()
#            plt.xlim(index[0]-0.2,(index + bar_width)[-1]+bar_width+0.2)
            plt.tight_layout()

            ax2=fig.add_subplot(1,2,2)
            ax2.bar(index, means_sample_1- means_sample_2, bar_width+0.35,
                             alpha=opacity,
                             color=(0,1,1,0.4),
                             label='')


            plt.ylabel('% Sequences Difference')
            plt.title(get_plot_title_with_ref_name('%s - %s' % (args.sample_1_name,args.sample_2_name),amplicon_name))
            plt.xticks(index,['Unmodified', 'Modified'])


#            plt.xlim(index[0]-bar_width/2, (index+bar_width)[-1]+2*bar_width)
            plt.tight_layout()
            plt.savefig(_jp('1.'+amplicon_name+'.Comparison_Efficiency.pdf'), bbox_inches='tight')
            if args.save_also_png:
                plt.savefig(_jp('1.'+amplicon_name+'.Comparison_Efficiency.png'), bbox_inches='tight')


            #profile comparion
            fig=plt.figure(figsize=(20,10))

            ax1=fig.add_subplot(1,2,1)
            plt.title(get_plot_title_with_ref_name('Mutation position distribution',amplicon_name))
            y_max=max(effect_vector_any_1.max(),effect_vector_any_2.max())*1.2

            plt.plot(effect_vector_any_1,color=(0,0,1,0.3),lw=4,label='%s combined mutations' % args.sample_1_name)
#            plt.hold(True)
            plt.plot(effect_vector_any_2,color=(1,0,0,0.3),lw=4,label='%s combined mutations' % args.sample_2_name)

            if cut_points:
                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                            plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2,label='Predicted cleavage position')
                    else:
                            plt.plot([cut_point,cut_point],[0,y_max],'--k',lw=2,label='_nolegend_')


                for idx,sgRNA_int in enumerate(sgRNA_intervals):
                    if idx==0:
                       plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='sgRNA')
                    else:
                       plt.plot([sgRNA_int[0],sgRNA_int[1]],[0,0],lw=10,c=(0,0,0,0.15),label='_nolegend_')


            lgd=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.3),ncol=1, fancybox=True, shadow=False)


            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
            plt.xlabel('Reference amplicon position (bp)')
            plt.ylabel('Sequences %')
            plt.ylim(0,max(1,y_max))
            plt.xlim(xmax=len_amplicon-1)

            ax2=fig.add_subplot(1,2,2)

            effect_vector_any_diff=effect_vector_any_1-effect_vector_any_2

            y_max=effect_vector_any_diff.max()*1.2
            y_min=effect_vector_any_diff.min()*1.2

            plt.title(get_plot_title_with_ref_name('%s - %s' % (args.sample_1_name,args.sample_2_name),amplicon_name))
            plt.plot(effect_vector_any_diff,color=(0,1,0,0.4),lw=3,label='Difference' )


            if cut_points:
                for idx,cut_point in enumerate(cut_points):
                    if idx==0:
                            plt.plot([cut_point,cut_point],[min(-1,y_min),max(1,y_max)],'--k',lw=2,label='Predicted cleavage position')
                    else:
                            plt.plot([cut_point,cut_point],[min(-1,y_min),max(1,y_max)],'--k',lw=2,label='_nolegend_')


                for idx,sgRNA_int in enumerate(sgRNA_intervals):
                    if idx==0:
                       plt.plot([sgRNA_int[0],sgRNA_int[1]],[min(-1,y_min),min(-1,y_min)],lw=10,c=(0,0,0,0.15),label='sgRNA')
                    else:
                       plt.plot([sgRNA_int[0],sgRNA_int[1]],[min(-1,y_min),min(-1,y_min)],lw=10,c=(0,0,0,0.15),label='_nolegend_')

            lgd2=plt.legend(loc='center', bbox_to_anchor=(0.5, -0.2),ncol=1, fancybox=True, shadow=False)
            plt.xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int) )
            plt.xlabel('Reference amplicon position (bp)')
            plt.ylabel('Sequences Difference %')
            plt.xlim(xmax=len_amplicon-1)

            plt.ylim(min(-1,y_min),max(1,y_max))

            plt.savefig(_jp('2.'+amplicon_name+'.Comparison_Combined_Insertion_Deletion_Substitution_Locations.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
            if args.save_also_png:
                    plt.savefig(_jp('2.'+amplicon_name+'.Comparison_Insertion_Deletion_Substitution_Locations.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')

            mod_file_1 = amplicon_info_1[amplicon_name]['modification_count_file']
            amp_seq_1,mod_freqs_1 = CRISPRessoShared.parse_count_file(mod_file_1)
            mod_file_2 = amplicon_info_2[amplicon_name]['modification_count_file']
            amp_seq_2,mod_freqs_2 = CRISPRessoShared.parse_count_file(mod_file_2)
            consensus_sequence = amp_seq_1
            if amp_seq_2 != consensus_sequence:
                raise DifferentAmpliconLengthException('Different amplicon lengths for the two amplicons.')


            for mod in ['Insertions','Deletions','Substitutions','All_modifications']:

                mod_counts_1 = np.array(mod_freqs_1[mod], dtype=float)
                tot_counts_1 = np.array(mod_freqs_1['Total'],dtype=float)
                unmod_counts_1 = tot_counts_1 - mod_counts_1

                mod_counts_2 = np.array(mod_freqs_2[mod], dtype=float)
                tot_counts_2 = np.array(mod_freqs_2['Total'],dtype=float)
                unmod_counts_2 = tot_counts_2 - mod_counts_2

                fisher_results = [stats.fisher_exact([[z[0],z[1]],[z[2],z[3]]]) if max(z) > 0 else [nan,1.0] for z in zip(mod_counts_1,unmod_counts_1,mod_counts_2,unmod_counts_2)]
                oddsratios,pvalues = [ a for a,b in fisher_results ], [ b for a,b in fisher_results ]

                mod_df = []
                row = [args.sample_1_name+'_'+mod]
                row.extend(mod_counts_1)
                mod_df.append(row)

                row = [args.sample_1_name+'_total']
                row.extend(tot_counts_1)
                mod_df.append(row)

                row = [args.sample_2_name+'_'+mod]
                row.extend(mod_counts_2)
                mod_df.append(row)

                row = [args.sample_2_name+'_total']
                row.extend(tot_counts_2)
                mod_df.append(row)

                row = ['odds_ratios']
                row.extend(oddsratios)
                mod_df.append(row)

                row = ['pvalues']
                row.extend(pvalues)
                mod_df.append(row)

                colnames = ['Reference']
                colnames.extend(list(consensus_sequence))
                mod_df = pd.DataFrame(mod_df,columns=colnames)
#                mod_df = pd.concat([mod_df.iloc[:,0:2], mod_df.iloc[:,2:].apply(pd.to_numeric)],axis=1)
                #write to file
                mod_df.to_csv(_jp(amplicon_name + '.' + mod + '_quantification.txt'),sep='\t',index=None)

                #plot
                fig=plt.figure(figsize=(20,10))
                ax1 = fig.add_subplot(2,1,1)

                diff = np.divide(mod_counts_1,tot_counts_1) - np.divide(mod_counts_2,tot_counts_2)
                diff_plot = ax1.plot(diff,color=(0,1,0,0.4),lw=3,label='Difference' )
                ax1.set_title(get_plot_title_with_ref_name('%s: %s - %s' % (mod,args.sample_1_name,args.sample_2_name),amplicon_name))
                ax1.set_xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int))
                ax1.set_ylabel('Sequences Difference %')
                ax1.set_xlim(xmin=0,xmax=len_amplicon-1)

                pvalues = np.array(pvalues)
                min_nonzero = np.min(pvalues[np.nonzero(pvalues)])
                pvalues[pvalues == 0] = min_nonzero
                #ax2 = ax1.twinx()
                ax2 = fig.add_subplot(2,1,2)
                pval_plot = ax2.plot(-1*np.log10(pvalues),color=(1,0,0,0.4),lw=2,label='-log10 P-value')
                ax2.set_ylabel('-log10 P-value')
                ax2.set_xlim(xmin=0,xmax=len_amplicon-1)
                ax2.set_xticks(np.arange(0,len_amplicon,max(3,(len_amplicon/6) - (len_amplicon/6)%5)).astype(int))
                ax2.set_xlabel('Reference amplicon position (bp)')


                #bonferroni correction
                corrected_p = -1*np.log10(0.01/float(len(consensus_sequence)))
                cutoff_plot = ax2.plot([0,len(consensus_sequence)],[corrected_p,corrected_p],color='k',dashes=(5,10),label='Bonferronni corrected cutoff')

                plots = diff_plot + pval_plot + cutoff_plot

                diff_y_min,diff_y_max = ax1.get_ylim()
                p_y_min,p_y_max = ax2.get_ylim()
                if cut_points:
                    for idx,cut_point in enumerate(cut_points):
                        if idx==0:
                                plot_cleavage = ax1.plot([cut_point,cut_point],[diff_y_min,diff_y_max],'--k',lw=2,label='Predicted cleavage position')
                                ax2.plot([cut_point,cut_point],[p_y_min,p_y_max],'--k',lw=2,label='Predicted cleavage position')
                                plots = plots + plot_cleavage
                        else:
                                ax1.plot([cut_point,cut_point],[diff_y_min,diff_y_max],'--k',lw=2,label='_nolegend_')
                                ax2.plot([cut_point,cut_point],[diff_y_min,diff_y_max],'--k',lw=2,label='_nolegend_')


                    for idx,sgRNA_int in enumerate(sgRNA_intervals):
                        if idx==0:
                           p2 = ax1.plot([sgRNA_int[0],sgRNA_int[1]],[diff_y_min,diff_y_min],lw=10,c=(0,0,0,0.15),label='sgRNA')
                           ax2.plot([sgRNA_int[0],sgRNA_int[1]],[p_y_min,p_y_min],lw=10,c=(0,0,0,0.15),label='sgRNA')
                           plots = plots + p2
                        else:
                           ax1.plot([sgRNA_int[0],sgRNA_int[1]],[diff_y_min,diff_y_min],lw=10,c=(0,0,0,0.15),label='_nolegend_')
                           ax2.plot([sgRNA_int[0],sgRNA_int[1]],[p_y_min,p_y_min],lw=10,c=(0,0,0,0.15),label='_nolegend_')


                labs = [p.get_label() for p in plots]
                lgd=plt.legend(plots,labs,loc='upper center', bbox_to_anchor=(0.5, -0.2),ncol=1, fancybox=True, shadow=False)

                plt.savefig(_jp('2.'+amplicon_name+'.'+mod+'.quantification.pdf'),bbox_extra_artists=(lgd,), bbox_inches='tight')
                if args.save_also_png:
                    plt.savefig(_jp('2.'+amplicon_name+'.'+mod+'.quantification.png'),bbox_extra_artists=(lgd,), bbox_inches='tight')


            #create merged heatmaps for each cut site
            allele_files_1 = amplicon_info_1[amplicon_name]['allele_files']
            allele_files_2 = amplicon_info_2[amplicon_name]['allele_files']
            for allele_file_1 in allele_files_1:
                allele_file_1_name = os.path.split(allele_file_1)[1] #get file part of path
                for allele_file_2 in allele_files_2:
                    allele_file_2_name = os.path.split(allele_file_2)[1] #get file part of path
                    #if files are the same (same amplicon, cut site, guide), run comparison
                    if allele_file_1_name == allele_file_2_name:
                        df1 = pd.read_csv(allele_file_1,sep="\t")
                        df2 = pd.read_csv(allele_file_2,sep="\t")

                        #find unmodified reference for comparison (if it exists)
                        ref_seq_around_cut = ""
                        if len(df1.loc[df1['Reference_Sequence'].str.contains('-')==False]) > 0:
                            ref_seq_around_cut = df1.loc[df1['Reference_Sequence'].str.contains('-')==False]['Reference_Sequence'].iloc[0]
                        #otherwise figure out which sgRNA was used for this comparison
                        elif len(df2.loc[df2['Reference_Sequence'].str.contains('-')==False]) > 0:
                            ref_seq_around_cut = df2.loc[df2['Reference_Sequence'].str.contains('-')==False]['Reference_Sequence'].iloc[0]
                        else:
                            seq_len = df2[df2['Unedited']==True]['Reference_Sequence'].iloc[0]
                            for sgRNA_interval,cut_point in zip(sgRNA_intervals,cut_points):
                                sgRNA_seq = consensus_sequence[sgRNA_interval[0]:sgRNA_interval[1]]
                                if sgRNA_seq in allele_file_1_name:
                                    this_sgRNA_seq = sgRNA_seq
                                    this_cut_point = cut_point
                                    ref_seq_around_cut=consensus_sequence[max(0,this_cut_point-args.offset_around_cut_to_plot+1):min(len(reference_seq),cut_point+args.offset_around_cut_to_plot+1)]
                                    break

                        merged = pd.merge(df1, df2, on = ['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated'],suffixes=('_' + args.sample_1_name,'_'+args.sample_2_name),how='outer')
                        quant_cols = ['#Reads_'+args.sample_1_name,'%Reads_'+args.sample_1_name,'#Reads_'+args.sample_2_name,'%Reads_'+args.sample_2_name]
                        merged[quant_cols] = merged[quant_cols].fillna(0)
                        lfc_error =0.1
                        merged['each_LFC'] = np.log2(((merged['%Reads_'+args.sample_1_name]+lfc_error)/(merged['%Reads_'+args.sample_2_name]+lfc_error)).astype(float)).replace([np.inf,np.NaN],0)
                        merged = merged.reset_index().set_index('Aligned_Sequence')
                        output_root = allele_file_1_name.replace(".txt","")
                        merged.to_csv(_jp(output_root+".txt"),sep="\t",index=None)
                        CRISPRessoPlot.plot_alleles_table_compare(ref_seq_around_cut,merged.sort_values(['each_LFC'],ascending=True),args.sample_1_name,args.sample_2_name,_jp('3.'+output_root+"_top"),
                                    MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot,MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot,SAVE_ALSO_PNG=args.save_also_png)
                        CRISPRessoPlot.plot_alleles_table_compare(ref_seq_around_cut,merged.sort_values(['each_LFC'],ascending=False),args.sample_1_name,args.sample_2_name,_jp('3.'+output_root+"_bottom"),
                                    MIN_FREQUENCY=args.min_frequency_alleles_around_cut_to_plot,MAX_N_ROWS=args.max_rows_alleles_around_cut_to_plot,SAVE_ALSO_PNG=args.save_also_png)


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
