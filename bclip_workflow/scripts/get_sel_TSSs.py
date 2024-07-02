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

def main():
    """ Get chr-separated .bed files out of the .gtf file"""

    __doc__ = "Get chr-separated .bed files out of the .gtf file"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_gtf",
                        dest="input_gtf_file_path",
                        help="Path to the input gtf annotation file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--output_dir",
                        dest="output_dir",
                        help="Path to the input genome fai index file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--borders",
                        dest="borders",
                        help="string of 6 border coordinates, signal, conrol_up, control_down",
                        required=True,
                        metavar="FILE")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    # create output directory if it doesn't exist
    output_dir = options.output_dir+'/'
    out = subprocess.check_output('mkdir -p '+output_dir, shell=True)
    
    gtf = pd.read_csv(options.input_gtf_file_path,delimiter="\t",index_col=None,header=None)
    # expected format is 0,150,-150,0,150,300
    borders = [int(elem) for elem in (options.borders).split(',')]
    signal_borders,control_up_borders,control_down_borders = borders[:2],borders[2:4],borders[4:6]
    
    min_distance_btw_TSSs = max(control_down_borders)-min(control_up_borders)

    # get TSSs from annotation (the most distal TSS)
    transcripts = gtf.loc[gtf[2]=='transcript'][[0,3,4,6,8]]
    transcripts['gene_id'] = transcripts[8].str.split('gene_id "',expand=True)[1].str.split('";',expand=True)[0].astype('int32')
    transcripts['gene_name'] = transcripts[8].str.split('gene_name "',expand=True)[1].str.split('";',expand=True)[0]
    transcripts['gene_type'] = transcripts[8].str.split('gene_biotype "|gene_type "',expand=True)[1].str.split('";',expand=True)[0]
    transcripts = transcripts.rename(columns = {0:'chr'})
    transcripts = transcripts.drop([8],1)
    transcripts['TSS'] = transcripts.apply(lambda x:x[3] if x[6]=='+' else x[4],1).astype('int')
    transcripts['len'] = transcripts[4]-transcripts[3]+1
    
    transcripts[6] = transcripts[6].astype('category')
    transcripts = transcripts.rename(columns={6:'strand'})
    
    transcripts['chr'] = transcripts['chr'].astype('category')
    transcripts['gene_type'] = transcripts['gene_type'].astype('category')
    transcripts = transcripts.loc[(transcripts['gene_type']=='protein_coding')].reset_index(drop=True)
    
    updated_transcripts = []
    for strand in ['+','-']:
        tmp = transcripts.loc[transcripts['strand']==strand].sort_values(['chr','TSS'],ascending=[True,(True if strand=='+' else False)]).reset_index(drop=True)
        coef = (1 if strand=='+' else -1)
        # prioritize the most distal TSS if there are alternatives
        tmp = tmp.drop_duplicates('gene_name').reset_index(drop=True)
        a = []
        chr_prev,TSS_prev = '',-1
        for v in tmp[['chr','TSS']].values:
            if v[1]*coef>TSS_prev*coef+min_distance_btw_TSSs or v[0]!=chr_prev:
                a.append(1)
            else:
                a.append(0)
            chr_prev,TSS_prev = v[0],v[1]
        tmp['TSS_suites'] = a
        updated_transcripts.append(tmp)
    transcripts = pd.concat(updated_transcripts).reset_index(drop=True)
    
    sel_TSSs_bed = transcripts.loc[(transcripts['TSS_suites']==1)&(~transcripts['gene_name'].isna())][['chr','strand','TSS','gene_id','gene_name']].drop_duplicates().reset_index(drop=True)
    
    sel_TSSs_bed_signal = sel_TSSs_bed.copy()
    sel_TSSs_bed_signal['start'] = sel_TSSs_bed_signal.apply(lambda x:x['TSS']+(min(signal_borders) if x['strand']=='+' else max(signal_borders)*(-1)),1)
    sel_TSSs_bed_signal['end'] = sel_TSSs_bed_signal.apply(lambda x:x['TSS']+(max(signal_borders) if x['strand']=='+' else min(signal_borders)*(-1)),1)
    sel_TSSs_bed_signal['name'] = 'signal'+';'+sel_TSSs_bed_signal['TSS'].astype('str')+';'+sel_TSSs_bed_signal['gene_id'].astype('str')+';'+sel_TSSs_bed_signal['gene_name']
    
    sel_TSSs_bed_control_up = sel_TSSs_bed.copy()
    sel_TSSs_bed_control_up['start'] = sel_TSSs_bed_control_up.apply(lambda x:x['TSS']+(min(control_up_borders) if x['strand']=='+' else max(control_up_borders)*(-1)),1)
    sel_TSSs_bed_control_up['end'] = sel_TSSs_bed_control_up.apply(lambda x:x['TSS']+(max(control_up_borders) if x['strand']=='+' else min(control_up_borders)*(-1)),1)
    sel_TSSs_bed_control_up['name'] = 'control_up'+';'+sel_TSSs_bed_control_up['TSS'].astype('str')+';'+sel_TSSs_bed_control_up['gene_id'].astype('str')+';'+sel_TSSs_bed_control_up['gene_name']
    
    sel_TSSs_bed_control_down = sel_TSSs_bed.copy()
    sel_TSSs_bed_control_down['start'] = sel_TSSs_bed_control_down.apply(lambda x:x['TSS']+(min(control_down_borders) if x['strand']=='+' else max(control_down_borders)*(-1)),1)
    sel_TSSs_bed_control_down['end'] = sel_TSSs_bed_control_down.apply(lambda x:x['TSS']+(max(control_down_borders) if x['strand']=='+' else min(control_down_borders)*(-1)),1)
    sel_TSSs_bed_control_down['name'] = 'control_down'+';'+sel_TSSs_bed_control_down['TSS'].astype('str')+';'+sel_TSSs_bed_control_down['gene_id'].astype('str')+';'+sel_TSSs_bed_control_down['gene_name']
    
    sel_TSSs_bed = pd.concat([sel_TSSs_bed_signal,sel_TSSs_bed_control_up,sel_TSSs_bed_control_down]).reset_index(drop=True)
    
    sel_TSSs_bed['score'] = 0
    for chr in list(sel_TSSs_bed['chr'].unique()):    
        sel_TSSs_bed.loc[sel_TSSs_bed['chr']==chr][['chr','start','end','name','score','strand']].to_csv(output_dir+chr+'.nonsorted.bed', sep=str('\t'),header=False,index=None)
        command = 'bedtools sort -i '+output_dir+chr+'.nonsorted.bed > '+output_dir+chr+'.bed && rm '+output_dir+chr+'.nonsorted.bed'
        out = subprocess.check_output(command, shell=True)
    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)