#!/usr/bin/env python

###############################################################################
# Setup
###############################################################################

# Import modules
# ============================================================

import argparse, gzip, numpy, os, pysam, sys
from collections import defaultdict
from scipy.stats import norm


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='LoReM (Long Read Misassemblies) bins and identifies regions of the genome assembly with low spanning coverage and high number of read terminals (starts and ends) using 10X linked reads, PacBio reads, or Nanopore reads.', epilog='For 10X linked reads, the input bam is the sorted and indexed bam with a BX optional field, such as the output from longranger align. For PacBio or Nanopore reads, the input bam is the sorted and indexed bam from your choice of aligner.')
parser.add_argument("-b", "--binsize", metavar='INT', help="bin size for analysis of sequence [1000]", type=int, default=1000)
parser.add_argument("-l", "--length", metavar='INT', help="minimum 10X molecule or PacBio/Nanopore read length [1000]", type=int, default=1000)
parser.add_argument("-o", "--output", metavar='STR', help="output file prefix [output]", type=str, default='output')
parser.add_argument("-s", "--spanning", metavar='FLOAT', help="spanning coverage threshold: calculates from p-value if <1; otherwise sets cutoff value [0.05]", type=float, default=0.05)
parser.add_argument("-t", "--terminal", metavar='FLOAT', help="minimum terminal (starts/ends) threshold: calculates from p-value if <1; otherwise sets cutoff value [0.05]", type=float, default=0.05)
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 1.0')
tenx = parser.add_argument_group('10X optional arguments')
tenx.add_argument("-d", "--distance", metavar='INT', help="minimum distance between reads with same 10X barcode being called new 10X molecules [50000]", type=int, default=50000)
tenx.add_argument("-r", "--reads", metavar='INT', help="minimum reads required per 10X molecule [4]", type=int, default=4)
required = parser.add_argument_group('required arguments')
required.add_argument("-f", "--fasta", metavar='STR', help="genome fasta", type=str, required=True)
required.add_argument("-i", "--in_bam", metavar='STR', help="sorted and indexed bam", type=str, required=True)
required.add_argument("-y", "--type", metavar='STR', help="type of input: 10X, PacBio, or Nanopore", type=str, required=True)
args = parser.parse_args()


###############################################################################
# Functions
###############################################################################

# Analyze genome
def analyze_fai():
    count = 0
    bases = 0
    bins = {}
    for fasta in pysam.FastxFile(args.fasta):
        len_fasta = len(str(fasta.sequence))
        bases += len_fasta
        bins[fasta.name] = int(float(len_fasta)/args.binsize)
        count += 1
    out_log.write('Extracted '+str(count)+' scaffolds from '+args.fasta+'\n')
    return bases, bins


# Analyze 10X barcode
def analyze_barcode(read, bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density):
    bx = read.get_tag('BX')
    if 'BX:Z:' in bx:
        bx = bx.split('BX:Z:')[1].split(' ')[0]
    if bx not in bx_seen:
        bx_seen[bx] = True
        bx_index[bx] = 0
        bx_int = bx+':'+str(bx_index[bx])
        mi_chrom[bx_int] = read.reference_name
        mi_start[bx_int] = read.reference_start
        mi_end[bx_int] = read.reference_end
        mi_reads[bx_int] = 1
        mi_coverage[bx_int] = read.reference_end - read.reference_start
        mi_density[bx_int].append(read.reference_start)
    elif ((read.reference_start - mi_end[bx+':'+str(bx_index[bx])]) >= args.distance and read.reference_name == mi_chrom[bx+':'+str(bx_index[bx])]) or read.reference_name != mi_chrom[bx+':'+str(bx_index[bx])]:
        bx_index[bx] += 1
        bx_int = bx+':'+str(bx_index[bx])
        mi_chrom[bx_int] = read.reference_name
        mi_start[bx_int] = read.reference_start
        mi_end[bx_int] = read.reference_end
        mi_reads[bx_int] = 1
        mi_coverage[bx_int] = read.reference_end - read.reference_start
        mi_density[bx_int].append(read.reference_start)
    else:
        bx_int = bx+':'+str(bx_index[bx])
        mi_end[bx_int] = read.reference_end
        mi_reads[bx_int] += 1
        mi_coverage[bx_int] = mi_coverage[bx_int] + read.reference_end - read.reference_start
        if read.reference_start not in mi_density[bx_int]:
            mi_density[bx_int].append(read.reference_start)
    old_chrom = read.reference_name
    return bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density


# Clear 10X variables
def clear_variables():
    bx_seen = {}
    bx_index = {}
    mi_chrom = {}
    mi_start = {}
    mi_end = {}
    mi_reads = {}
    old_chrom = False
    temp_list = []
    mi_coverage = {}
    mi_density = defaultdict(list)
    return bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, temp_list, mi_coverage, mi_density


# Analyze 10X bam file
def analyze_10X_bam():
    bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, temp_list, mi_coverage, mi_density = clear_variables()
    bed = []
    reads_count = 0
    mi_count = 0
    out_bed = open(args.output+'.full.bed', 'w')

    # Read bam file and extract BX, chrom, start, and end
    for read in pysam.AlignmentFile(args.in_bam, 'rb').fetch():

        if read.has_tag('BX') and read.reference_name and read.reference_start and read.reference_end:
            if not old_chrom or read.reference_name == old_chrom:
                bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density = analyze_barcode(read, bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density)
            else:
                for mi in sorted(mi_chrom, key=mi_end.get):                    
                    if mi_reads[mi] >= args.reads and (mi_end[mi]-mi_start[mi]) >= args.length:
                        temp_list.append(mi)
                for mi in sorted(temp_list, key=mi_start.get):
                    out_bed.write(mi_chrom[mi]+'\t'+str(mi_start[mi])+'\t'+str(mi_end[mi])+'\t'+str(mi)+'\t'+str(float(mi_coverage[mi])/(mi_end[mi]-mi_start[mi]))+','+str(float(len(mi_density[mi]))/(mi_end[mi]-mi_start[mi]))+'\n')
                    bed.append(mi_chrom[mi]+':'+str(mi_start[mi])+':'+str(mi_end[mi]))
                    mi_count += 1
                bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, temp_list, mi_coverage, mi_density = clear_variables()
                bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density = analyze_barcode(read, bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density)
            reads_count += 1
    for mi in sorted(mi_chrom, key=mi_end.get):
        if mi_reads[mi] >= args.reads and (mi_end[mi]-mi_start[mi]) >= args.length:
            temp_list.append(mi)
    for mi in sorted(temp_list, key=mi_start.get):
        out_bed.write(mi_chrom[mi]+'\t'+str(mi_start[mi])+'\t'+str(mi_end[mi])+'\t'+str(mi)+'\t'+str(float(mi_coverage[mi])/(mi_end[mi]-mi_start[mi]))+','+str(float(len(mi_density[mi]))/(mi_end[mi]-mi_start[mi]))+'\n')
        bed.append(mi_chrom[mi]+':'+str(mi_start[mi])+':'+str(mi_end[mi]))
        mi_count += 1
    bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, temp_list, mi_coverage, mi_density = clear_variables()
    out_log.write('Extracted '+str(reads_count)+' reads from '+args.in_bam+'\nIdentified '+str(mi_count)+' 10X molecules with >='+str(args.reads)+' supporting reads and >='+str(args.length)+' length\n')
    out_bed.flush()
    out_bed.close()
    pysam.tabix_index(args.output+'.full.bed', preset='bed')
    return bed


# Analyze PacBio/Nanopore bam file
def analyze_other_bam():

    # Read bam file and extract chrom, start, and end and write reads to bed
    bed = []
    count = 0
    for read in pysam.AlignmentFile(args.in_bam, 'rb').fetch():
        if (read.reference_end - read.reference_start) >= args.length:
            bed.append(read.reference_name+':'+str(read.reference_start)+':'+str(read.reference_end))
            count += 1
    out_log.write('Extracted '+str(count)+' reads from '+args.in_bam+'\n')
    return bed


# Set threshold
def set_threshold(bins_dict, bases, argsvalue):
    if args.spanning >= 1:
        initial_threshold = argsvalue
    else:
        all_sites = []
        p_value_success = False
        for site in bins_dict:
            if bins_dict[site] > 0:
                all_sites.append(bins_dict[site])
        all_sites_mean = numpy.mean(all_sites)
        all_sites_std = numpy.std(all_sites)
        for n in range(1,int(all_sites_mean)):
            if norm.sf(abs(((n-all_sites_mean)/all_sites_std))) < argsvalue:
                initial_threshold = n
                p_value_success = True
                break
        if not p_value_success:
            all_sites_sum = numpy.sum(all_sites)
            for n in range(1,int(numpy.amax(all_sites))):
                if (1 - numpy.exp(-(float(all_sites_sum)/bases)*n)) > argsvalue:
                    initial_threshold = n
                    p_value_success = True
                    break
        if not p_value_success:
            out_log.write('Unable to determine threshold based on p value; setting to 10.\n')
            initial_threshold = 10
    return initial_threshold


# Analyze bed output
def analyze_bed(bed, bases, bins):
    bins_starts = {}
    bins_ends = {}
    bins_coverage = {}

    # Set all start, end, and coverage for bins to 0
    for site in bins:
        for n in range(0, bins[site]+1):
            bins_starts[site+':'+str(n)] = 0
            bins_ends[site+':'+str(n)] = 0
            bins_coverage[site+':'+str(n)] = 0

    # Extract start, end, and coverage for bins
    count_start = 0
    count_end = 0
    count_coverage = 0
    for site in bed:
        site = site.split(':')
        start_bin = int(float(site[1])/args.binsize)
        end_bin = int(float(site[2])/args.binsize)
        if site[0]+':'+str(start_bin) in bins_starts and site[0]+':'+str(end_bin) in bins_ends:
            if (end_bin - start_bin) > 1:
                for n in range(start_bin+1, end_bin):
                    bins_coverage[site[0]+':'+str(n)] += 1
                count_coverage += 1
            bins_starts[site[0]+':'+str(start_bin)] += 1
            bins_ends[site[0]+':'+str(end_bin)] += 1
            count_start += 1
            count_end += 1
    out_log.write('Identified '+str(count_start)+' long read starts, '+str(count_end)+' long read ends, and '+str(count_coverage)+' long reads with spanning coverage\n')


    # Calculate thresholds
    coverage_threshold = set_threshold(bins_coverage, bases, args.spanning)
    start_treshold = set_threshold(bins_starts, bases, args.terminal)
    end_treshold = set_threshold(bins_ends, bases, args.terminal)
    out_log.write('Set thresholds to coverage '+str(coverage_threshold)+', start '+str(start_treshold)+', and end '+str(end_treshold)+'\n')

    # Determine bins that meet thresholds
    out_tsv = gzip.open(args.output+'.tsv.gz', 'wb')
    out_tsv.write('CHROM\tSTART\tSTOP\tCOVERAGE\tSTARTS\tENDS\n')
    out_bed = open(args.output+'.bed', 'w')
    for site in bins:
        for n in range(0, bins[site]+1):
            if bins_coverage[site+':'+str(n)] > 0 or bins_starts[site+':'+str(n)] > 0 or bins_ends[site+':'+str(n)] > 0:
                out_tsv.write(site+'\t'+str(args.binsize*n)+'\t'+str(args.binsize*(1+n))+'\t'+str(bins_coverage[site+':'+str(n)])+'\t'+str(bins_starts[site+':'+str(n)])+'\t'+str(bins_ends[site+':'+str(n)])+'\n')
            if bins_coverage[site+':'+str(n)] < coverage_threshold and bins_starts[site+':'+str(n)] >= start_treshold and bins_ends[site+':'+str(n)] >= end_treshold:
                out_bed.write(site+'\t'+str(args.binsize*n)+'\t'+str(args.binsize*(1+n))+'\tcoverage:'+str(bins_coverage[site+':'+str(n)])+';starts:'+str(bins_starts[site+':'+str(n)])+';ends:'+str(bins_ends[site+':'+str(n)])+'\n')
    out_tsv.flush()
    out_tsv.close()
    out_bed.flush()
    out_bed.close()


###############################################################################
# Run
###############################################################################

# Run functions
# ============================================================

out_log = sys.stderr
bases, bins = analyze_fai()
if args.type.upper() == '10X':
    bed = analyze_10X_bam()
elif args.type.upper() == 'PACBIO' or args.type.upper() == 'NANOPORE':
    bed = analyze_other_bam()
else:
    out_log.write('Unable to accept read type '+args.type+'. Please select 10X, PacBio, or Nanopore.\n')
    sys.exit(2)
analyze_bed(bed, bases, bins)
out_log.flush()
out_log.close()
