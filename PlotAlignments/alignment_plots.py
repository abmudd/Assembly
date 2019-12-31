#!/usr/bin/env python

# This script is derived from a Jupyter notebook written by Sofia Medina Ruiz.

# Import modules
# ============================================================

import argparse, os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script writes out alignments of three species chromosomes with centromere locations and repeat histograms.')
parser.add_argument("-o", "--output", metavar='STR', help="output directory [pwd]", type=str, default=os.getcwd()+'/')
parser.add_argument("-x", "--xtropicalis", help="flag if mid is Xtropicalis", action='store_true')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("-c", "--centromere", metavar='STR', help="concat centromere positions", type=str, required=True)
required.add_argument("-f", "--fai", metavar='STR', help="concat fai", type=str, required=True)
required.add_argument("-l", "--links", metavar='STR', help="links bedpe", type=str, required=True)
required.add_argument("-r", "--repeatbins", metavar='STR', help="concat repeat bins output from gff2bins.py", type=str, required=True)
required.add_argument("-s", "--seqconfig", metavar='STR', help="config with all included sequence", type=str, required=True)
args = parser.parse_args()

# Set variables
# ============================================================

Mb = 1000000
int_ = 0.5;
step = 0.25;
fig = plt.figure(figsize=(8,3))
extra_space = step /4
ymin_ = -10;
ymax_ = 12
qtile_ = 0.1
xtr_chr = {
    "Xtropicalis.Chr1": "89200000",
    "Xtropicalis.Chr2": "67500000",
    "Xtropicalis.Chr3": "16700000",
    "Xtropicalis.Chr4": "46600000",
    "Xtropicalis.Chr5": "62000000",
    "Xtropicalis.Chr6": "73100000",
    "Xtropicalis.Chr7": "60300000",
    "Xtropicalis.Chr8": "21500000",
    "Xtropicalis.Chr9": "42100000",
    "Xtropicalis.Chr10": "21200000"
}

# Define functions
# ============================================================

def str2bool(v):
    return v.lower() in ("true")

def Obtain_Chr_size(Fai_lengths, Query_Chr):
    chr_size = int(Fai_lengths[Fai_lengths.Src==Query_Chr].Lengths.tolist()[0])
    return(chr_size)

def place_10Mb_ticks(chr_length, add_extraMbs, position_y, Rev, factor, len_max):
    if len_max > 200*Mb:
        ticks_10Mb = np.arange(0, chr_length, 20*Mb)
        for i in ticks_10Mb:
            if i < chr_length/20*Mb:
                x_pos = i/Mb
                if Rev == True: 
                    x_pos = (chr_length- i)/Mb
                ax5.plot([(x_pos+add_extraMbs)*factor, (x_pos+add_extraMbs)*factor], [position_y-.2, position_y], color='black', linewidth=1, alpha=0.5)
                ax5.text((x_pos+add_extraMbs)*factor, position_y-.6, str(int(i/Mb)), fontsize=8, horizontalalignment='center')
    else:
        ticks_10Mb = np.arange(0, chr_length, 10*Mb)
        for i in ticks_10Mb:
            if i < chr_length/10*Mb:
                x_pos = i/Mb
                if Rev == True: 
                    x_pos = (chr_length- i)/Mb
                ax5.plot([(x_pos+add_extraMbs)*factor, (x_pos+add_extraMbs)*factor], [position_y-.2, position_y], color='black', linewidth=1, alpha=0.5)
                ax5.text((x_pos+add_extraMbs)*factor, position_y-.6, str(int(i/Mb)), fontsize=8, horizontalalignment='center')
    return(ax5)

def Show_chr_backbone(Chr, chr_length, incolor, add_extraMbs, position_y, Rev, factor, second, len_max):
    Chr_name = Chr
    if Rev == True:
        Chr_name = ' '.join((Chr,'(Rev)'))
    if second == 0:
        ax5.text((add_extraMbs-1)*factor, position_y, Chr_name, fontsize=12, color=incolor, horizontalalignment='right')
    else:
        ax5.text((add_extraMbs)*factor, position_y+second, Chr_name, fontsize=12, color=incolor, horizontalalignment='left')
    ax5.plot([(add_extraMbs)*factor, (add_extraMbs+float(chr_length)/Mb)*factor], [position_y,position_y], color='gray', alpha=.5, linewidth=2)
    place_10Mb_ticks(chr_length, add_extraMbs, position_y, Rev, factor, len_max)
    return(ax5)

def show_repeat(Repeat_bins, Query_Chr, incolor, add_extraMbs, position_y, chr_length, Rev, factor):
    Rep_pos = Repeat_bins[Repeat_bins.Src==Query_Chr].Position
    if Rev==True:
        Rep_pos = (chr_length - Rep_pos)/Mb
    else:
        Rep_pos = Rep_pos/Mb
    Repeat_distribution = Repeat_bins[Repeat_bins.Src==Query_Chr]['Repeats']/Repeat_bins[Repeat_bins.Src==Query_Chr]['Repeats'].max()
    ax5.fill_between((add_extraMbs+Rep_pos)*factor, position_y, position_y+2*Repeat_distribution, color=incolor, alpha=0.4, linewidth=0)
    return(ax5)

def draw_centromere(Qry_Chr, Chromosome_size, add_extraMbs, position_y, Rev, factor, incolor):
    Centromere = int(Centromere_postions[Qry_Chr])
    if Rev == True: 
        Centromere = Chromosome_size - Centromere
    ax5.plot((add_extraMbs+Centromere/Mb)*factor, position_y, '*', color=incolor, alpha=0.85, markersize=12)
    return(ax5)

def draw_Xtr_centromere(Qry_Chr, Chromosome_size, add_extraMbs, position_y, Rev, factor, incolor):
    Centromere = int(xtr_chr[Qry_Chr])
    if Rev == True: 
        Centromere = Chromosome_size - Centromere
    ax5.plot((add_extraMbs+Centromere/Mb)*factor, position_y, '*', color=incolor, alpha=0.85, markersize=12)
    return(ax5)

def draw_links(Qry_Tar_DF, Query_Chr, Target_Chr, Query_chr_length, incolor, add_extraMbs, position_y, Rev, Ref_add_extraMbs, factor):
    Qry_starts = Qry_Tar_DF[(Qry_Tar_DF.Q_Src.isin([Query_Chr])) & (Qry_Tar_DF.T_Src.isin([Target_Chr]))]['Q_mid']
    Tar_starts = Qry_Tar_DF[(Qry_Tar_DF.Q_Src.isin([Query_Chr])) & (Qry_Tar_DF.T_Src.isin([Target_Chr]))]['T_mid']
    if Rev ==True: 
        Qry_starts = Query_chr_length - Qry_starts
    coordinates = zip(list(Qry_starts/Mb),list(Tar_starts/Mb))
    
    for i in range(0, len(coordinates)):
        a,b = coordinates[i]
        ax5.plot([(a+add_extraMbs)*factor,b+Ref_add_extraMbs], [position_y,8], alpha=0.1, color=incolor, linewidth=0.5)
    return(ax5)

def figure_aspects_():
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['left'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    return(ax5)

# Import files
# ============================================================
Centromere_postions = {} # import centromeres
for line in open(args.centromere):
    line = line.rstrip().split('\t')
    Centromere_postions[line[0]] = line[1]
Fai_lengths = pd.read_csv(args.fai, sep='\t', names=['Src','Lengths','NA1','NA2','NA3']) # import fai
Repeat_bins = pd.read_csv(args.repeatbins, sep='\t', names=['Src','Position','Repeats']) # import repeats
# import and parse links bedpe
A_B_link = pd.read_csv(args.links, sep='\t', names=['Q_Src','Q_Start','Q_End','T_Src','T_Start','T_End'])
A_B_link['Q_strand'] = '+'
A_B_link['T_strand'] = '+'
A_B_link['Q_mid'] = (A_B_link.Q_End + A_B_link.Q_Start)/2
A_B_link['T_mid'] = (A_B_link.T_End + A_B_link.T_Start)/2

# Import and parse config
# ============================================================ 
all_chr = []
high = {}
high_count = 0
high_length = 0
mid = False
mid_count = 0
mid_length = 0
low = {}
low_count = 0
low_length = 0
for line in open(args.seqconfig,'r'):
    if not line.startswith('#'):
        line = line.rstrip().split('\t')
        all_chr.append(line[0])
        if line[3].lower() == 'high':
            high_count += 1
            high[str(high_count)] = ';'.join(line[0:3])
            high_length += int(line[2])*Mb + Obtain_Chr_size(Fai_lengths, line[0])
        elif line[3].lower() == 'mid':
            mid_count += 1
            mid = ';'.join(line[0:3])
            mid_length += int(line[2])*Mb + Obtain_Chr_size(Fai_lengths, line[0])
        elif line[3].lower() == 'low':
            low_count += 1
            low[str(low_count)] = ';'.join(line[0:3])
            low_length += int(line[2])*Mb + Obtain_Chr_size(Fai_lengths, line[0])

# Plot mid
# ============================================================ 
if mid_count > 1:
    sys.err.write('ERROR: Unable to handle more than 1 middle chromosome at this time.\n')
    sys.exit(2)
elif mid_count == 0:
    sys.err.write('ERROR: Must include 1 middle chromosome.\n')
    sys.exit(2)
else:
    Ref_chr = mid.split(';')[0]; Rev = str2bool(mid.split(';')[1].lower());
    Ref_add_extraMbs = int(mid.split(';')[2]);
    position_y = 8; Ref_chr_size = Obtain_Chr_size(Fai_lengths, Ref_chr);
    ax5 = fig.add_axes([0.1, int_, 0.8, step*5],  xlim=(-10,Ref_chr_size/Mb)) 
    Show_chr_backbone(Ref_chr, Ref_chr_size, 'blue', Ref_add_extraMbs, position_y, Rev, 1, 0, mid_length)
    show_repeat(Repeat_bins, Ref_chr, 'orange', Ref_add_extraMbs, position_y, Ref_chr_size, Rev, 1)
    if args.xtropicalis:
        draw_Xtr_centromere(Ref_chr, Ref_chr_size, Ref_add_extraMbs, position_y, Rev, 1, 'red')
    draw_centromere(Ref_chr, Ref_chr_size, Ref_add_extraMbs, position_y, Rev, 1, 'black')

# Plot high
# ============================================================
high_extraMbs = 0
for i in range(1,high_count+1):
    Query_chr = high[str(i)].split(';')[0]; Rev = str2bool(high[str(i)].split(';')[1].lower());
    add_extraMbs = int(high[str(i)].split(';')[2]) + high_extraMbs
    position_y = 16; Query_chr_size = Obtain_Chr_size(Fai_lengths, Query_chr);
    if i == 1:
        Show_chr_backbone(Query_chr, Query_chr_size, 'green', add_extraMbs, position_y, Rev, float(mid_length)/high_length, 0, high_length)
    else:
        Show_chr_backbone(Query_chr, Query_chr_size, 'green', add_extraMbs, position_y, Rev, float(mid_length)/high_length, 1.6, high_length)
    show_repeat(Repeat_bins, Query_chr, 'orange', add_extraMbs, position_y, Query_chr_size, Rev, float(mid_length)/high_length)
    draw_centromere(Query_chr, Query_chr_size, add_extraMbs, position_y, Rev, float(mid_length)/high_length, 'black')
    draw_links(A_B_link, Query_chr, Ref_chr, Query_chr_size, 'gray', add_extraMbs, position_y, Rev, Ref_add_extraMbs, float(mid_length)/high_length)
    high_extraMbs += int(high[str(i)].split(';')[2]) + float(Query_chr_size)/Mb

# Plot low
# ============================================================ 
low_extraMbs = 0
for i in range(1,low_count+1):
    Query_chr = low[str(i)].split(';')[0]; Rev = str2bool(low[str(i)].split(';')[1].lower());
    add_extraMbs = int(low[str(i)].split(';')[2]) + low_extraMbs
    position_y = 0; Query_chr_size = Obtain_Chr_size(Fai_lengths, Query_chr);
    if i == 1:
        Show_chr_backbone(Query_chr, Query_chr_size, 'purple', add_extraMbs, position_y, Rev, float(mid_length)/low_length, 0, low_length)
    else:
        Show_chr_backbone(Query_chr, Query_chr_size, 'purple', add_extraMbs, position_y, Rev, float(mid_length)/low_length, -1.6, low_length)
    show_repeat(Repeat_bins, Query_chr, 'orange', add_extraMbs, position_y, Query_chr_size, Rev, float(mid_length)/low_length)
    draw_centromere(Query_chr, Query_chr_size, add_extraMbs, position_y, Rev, float(mid_length)/low_length, 'black')
    draw_links(A_B_link, Query_chr, Ref_chr, Query_chr_size, 'gray', add_extraMbs, position_y, Rev, Ref_add_extraMbs, float(mid_length)/low_length)
    low_extraMbs += int(low[str(i)].split(';')[2]) + float(Query_chr_size)/Mb

# Save plot
# ============================================================ 
figure_aspects_()
plt.savefig('-'.join(all_chr)+'.png', bbox_inches='tight')
