#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Aleksei Mironov
# Company: Mihaela Zavolan, Biozentrum, Basel
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numpy as np
import itertools

import seaborn as sns; sns.set()

import matplotlib
matplotlib.style.use('seaborn')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['backend'] = "Qt5Agg"
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def main():
    """ Get genomic segmentation out of the .gtf file"""

    __doc__ = "Get genomic segmentation out of the .gtf file"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--start_samples",
                        dest="start_samples_path",
                        help="Path to the input start samples file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--sample",
                        dest="sample",
                        help="sample id",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--binned_genome_dir",
                        dest="binned_genome_dir",
                        help="path for the binned genome mapping",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--selected_TSSs_over_chrs_dir",
                        dest="selected_TSSs_over_chrs_dir",
                        help="path for the directory with chr-separated sel TSSs",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--metaplot_dir",
                        dest="metaplot_dir",
                        help="path for the output directory",
                        required=True,
                        metavar="FILE")
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    def q25(x):
        return np.quantile(x,q=0.25)
    def q75(x):
        return np.quantile(x,q=0.75)
        
    start_samples = pd.read_csv(options.start_samples_path,delimiter="\t",index_col=None,header=0)
    sample = options.sample

    stats = [np.mean,q25,np.median,q75]
    stats_labels = ['mean','q25','median','q75']

    binned_genome_dir = options.binned_genome_dir+'/'
    metaplot_dir = options.metaplot_dir+'/'

    selected_TSSs_over_chrs_dir = options.selected_TSSs_over_chrs_dir+'/'
    
    chr_list = (start_samples.loc[start_samples['name']==sample].iloc[0]['chromosomes'].split(' '))

    command = 'mkdir -p '+metaplot_dir
    out = subprocess.check_output(command, shell=True)

    CC = 1
    for count_col in ['uniquely_mapped;0']:
        bin_counts_over_chrs_list = []
        for chr in chr_list:
            
            input_file_path = binned_genome_dir+sample+'.'+chr+'.wd_counts.tsv'
            if not (os.path.isfile(input_file_path) and (os.stat(input_file_path).st_size > 0)):
                print(chr+': no wd_counts file')
                continue
            
            selected_TSSs_over_chrs_file = selected_TSSs_over_chrs_dir+chr+'.bed'
            if not (os.path.isfile(selected_TSSs_over_chrs_file) and (os.stat(selected_TSSs_over_chrs_file).st_size > 0)):
                print(chr+': no selected_TSSs file')
                continue

            bin_counts = pd.read_csv(input_file_path,delimiter="\t",index_col=None,header=0,dtype={"chr": "category","strand":"category","segment_name":"category"},nrows=5)
            
            i=1
            col_number = -1
            for elem in bin_counts.columns:
                if elem==count_col:
                    col_number = i
                i=i+1
            if col_number==-1:
                score_val = "0"
            else:
                score_val = "$"+str(col_number)
            
            nonsorted_temp_file_path =metaplot_dir+sample+'.'+chr+'.nonsorted.bed'
            sorted_temp_file_path = metaplot_dir+sample+'.'+chr+'.sorted.bed'
            sel_file_path = metaplot_dir+sample+'.'+chr+'.selected.bed'
            
            command = """awk \'(NR>1) {print $1"\\t"$2"\\t"$3"\\t""t""\\t\""""+score_val+"""\"\\t"$4}' """+input_file_path+' > '+nonsorted_temp_file_path
            out = subprocess.check_output(command, shell=True)
    
            command = 'bedtools sort -i '+nonsorted_temp_file_path+' > '+sorted_temp_file_path
            out = subprocess.check_output(command, shell=True)   
            command = 'bedtools intersect -sorted -a '+sorted_temp_file_path+' -b '+selected_TSSs_over_chrs_file+' -wa -wb -s'+""" | awk \'{print $1"\\t"$2"\\t"$3"\\t"$10"\\t"$5"\\t"$6}\' > """+sel_file_path
            out = subprocess.check_output(command, shell=True)
    
            command = 'rm '+sorted_temp_file_path+' '+nonsorted_temp_file_path
            out = subprocess.check_output(command, shell=True)
            
            bin_counts = pd.read_csv(sel_file_path,delimiter="\t",index_col=None,header=None)
            bin_counts.columns = ['chr','start','end','name','r','strand']
            bin_counts['TSS'] = bin_counts['name'].astype('str').str.split(';',expand=True)[1].astype('int')
    
            bin_counts['pos_rel_to_TSS'] = ((bin_counts['end']+bin_counts['start'])*0.5-bin_counts['TSS']).astype('int')*((bin_counts['strand']=='+').astype('int')*2-1)
            bin_counts['pos_rel_to_TSS'] = (np.round((bin_counts['pos_rel_to_TSS']/10),0)*10).astype('int')
            bin_counts = bin_counts.drop_duplicates(['chr','start','end','strand']).reset_index(drop=True)
            bin_counts_over_chrs_list.append(bin_counts)
            
            command = 'rm '+sel_file_path
            out = subprocess.check_output(command, shell=True)            
            
            print(chr+' done')
        bin_counts_over_chrs_df = pd.concat(bin_counts_over_chrs_list).reset_index(drop=True)
        
        # remove the maximum and minimum positions - they are always enriched in metaplots
        min_val,max_val = bin_counts_over_chrs_df['pos_rel_to_TSS'].min(),bin_counts_over_chrs_df['pos_rel_to_TSS'].max()
        bin_counts_over_chrs_df = bin_counts_over_chrs_df.loc[~bin_counts_over_chrs_df['pos_rel_to_TSS'].isin([min_val,max_val])].reset_index(drop=True)
        
        if 'r_sum' in list(bin_counts_over_chrs_df.columns):
            bin_counts_over_chrs_df = bin_counts_over_chrs_df.drop('r_sum',1)
        gr = bin_counts_over_chrs_df.groupby('TSS').agg({'r':np.sum}).reset_index().rename(columns={'r':'r_sum'})
        bin_counts_over_chrs_df = pd.merge(bin_counts_over_chrs_df,gr,how='left',on='TSS')
        bin_counts_over_chrs_df['log2_r_sum'] = np.log2(bin_counts_over_chrs_df['r_sum']+1)
        
        # start plotting
        
        data = bin_counts_over_chrs_df
        
        expr_thresholds = [0,2,3,4,5]
        if max(expr_thresholds)<data['log2_r_sum'].max():
            expr_thresholds = expr_thresholds+[data['log2_r_sum'].max()]
        data['expr_cat'] = pd.cut(data['log2_r_sum'],bins=expr_thresholds,include_lowest=True)
        data['frac'] = (data['r'])/(data['r_sum']+1)

        data['r_pseudo'] = data['r']+1
        data['cpm'] = data['r_pseudo']*10**5/data['r_pseudo'].sum() # bin length is 10, therefore multiply by 10**5, not 10**6
        data['log2_cpm'] = np.log2(data['cpm'])

        # for density
        sns.set(font_scale=1)
        sns.set_style("white")
        fig, axes = plt.subplots(1+len(stats), 1, sharey=False, sharex=False, figsize=(1+2*(len(stats)),12),gridspec_kw={'height_ratios': [1]+[2]*len(stats)})
        
        gr = data[['TSS','r_sum','log2_r_sum']].drop_duplicates().reset_index()
        ax = sns.histplot(ax=axes[0],data=gr['log2_r_sum'],bins=40,stat='probability')
        ylims = (ax.get_ylim()[0],ax.get_ylim()[1])
        for expr_threshold in expr_thresholds:
            ax.vlines(x=expr_threshold,ymin=ylims[0],ymax=ylims[1],color='grey',linestyles='--',linewidth=1)
        ax.set(title='total # of reads + 1, $log_2$, '+count_col+', PDF',xlabel='')
        
        i=1
        for statistic in stats:
            x_feature,y_feature = 'pos_rel_to_TSS','frac'
            l = list(data['expr_cat'].unique().sort_values())
            palette = sns.color_palette("rocket",len(l))
            k=0
            for expr_cat in l:
                data_sel = data.loc[(data['expr_cat']==expr_cat)]
                data_sel['t']=1
                data_sel_gr = data_sel.groupby(x_feature).agg({y_feature:statistic,'t':np.sum}).reset_index()
                ax1 = sns.lineplot(ax=axes[i],data = data_sel_gr,x=x_feature,y=y_feature,color=palette[k],label=str(expr_cat)+', # of TSSs ='+str(int(data_sel_gr.iloc[0]['t'])))
                ax1.set(title=stats_labels[i-1]+' over TSSs, pos relative to TSS',xlabel='',ylabel='density of reads')
                k=k+1
            ylims = (ax1.get_ylim()[0],ax1.get_ylim()[1])
            ax1.vlines(x=0,ymin=ylims[0],ymax=ylims[1],color='grey',linestyles='--',linewidth=1)
            ax1.legend(bbox_to_anchor=(1, 0, 1, 1),loc='upper left',borderaxespad=0,title='total # of reads + 1, $log_2$',markerscale=1.2,ncols=1,fontsize=9)

            gr = data.groupby(['expr_cat',x_feature]).agg({y_feature:statistic}).reset_index().rename(columns={y_feature:stats_labels[i-1]})
            if i==1:
                table_results = gr.copy()
            else:
                table_results = pd.merge(table_results,gr,how='outer',on=['expr_cat',x_feature])
            i=i+1
                
        table_results.to_csv(metaplot_dir+sample+'.'+str(CC)+'.'+count_col.replace(';','_')+'.metaplot_density.tsv', sep=str('\t'),header=True,index=None)
        
        fig.tight_layout(pad=0.5)
        fig.savefig(metaplot_dir+sample+'.'+str(CC)+'.'+count_col.replace(';','_')+'.metaplot_density.pdf',bbox_inches='tight',dpi=300)
        
        # for CPM
        sns.set(font_scale=1)
        sns.set_style("white")
        fig, axes = plt.subplots(len(stats), 1, sharey=False, sharex=False, figsize=(2*(len(stats)),12),)

        i=0
        for statistic in stats:
            x_feature,y_feature = 'pos_rel_to_TSS','log2_cpm'
            data_sel_gr = data.groupby(x_feature).agg({y_feature:statistic}).reset_index()
            ax1 = sns.lineplot(ax=axes[i],data = data_sel_gr,x=x_feature,y=y_feature,color='black')
            ax1.set(title=stats_labels[i]+' over TSSs, pos relative to TSS',xlabel='',ylabel=r'$log_2CPM$')
            ylims = (ax1.get_ylim()[0],ax1.get_ylim()[1])
            ax1.vlines(x=0,ymin=ylims[0],ymax=ylims[1],color='grey',linestyles='--',linewidth=1)

            gr = data_sel_gr.rename(columns={y_feature:stats_labels[i]})
            if i==0:
                table_results = gr.copy()
            else:
                table_results = pd.merge(table_results,gr,how='outer',on=[x_feature])
            i=i+1
        table_results.to_csv(metaplot_dir+sample+'.'+str(CC)+'.'+count_col.replace(';','_')+'.metaplot_cpm.tsv', sep=str('\t'),header=True,index=None)
        fig.tight_layout(pad=0.5)
        fig.savefig(metaplot_dir+sample+'.'+str(CC)+'.'+count_col.replace(';','_')+'.metaplot_cpm.pdf',bbox_inches='tight',dpi=300)
        CC = CC+1
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)