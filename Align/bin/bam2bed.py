#!/usr/bin/env python

###############################################################################
# Setup
###############################################################################

# Import modules
# ============================================================

import argparse, os, pysam, sys
from collections import defaultdict


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script converts a bam to bed with 10X, Matepair, Nanopore, or PacBio reads for misassembly detection with Wombat.')
parser.add_argument("-l", "--length", metavar='INT', help="minimum 10X molecule, Matepair insert, or Nanopore/PacBio read length [1000]", type=int, default=1000)
parser.add_argument("-g", "--log", metavar='STR', help="output log [sys.stderr]", type=str, default='sys.stderr')
parser.add_argument("-o", "--output", metavar='STR', help="output file prefix [output]", type=str, default='output')
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.0')
tenx = parser.add_argument_group('10X optional arguments')
tenx.add_argument("-d", "--distance", metavar='INT', help="minimum distance between reads with same 10X barcode being called new 10X molecules [50000]", type=int, default=50000)
tenx.add_argument("-r", "--reads", metavar='INT', help="minimum reads required per 10X molecule [10]", type=int, default=10)
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--in_bam", metavar='STR', help="sorted and indexed bam", type=str, required=True)
required.add_argument("-t", "--type", metavar='STR', help="type of input: 10X, Matepair, Nanopore, or PacBio", type=str, required=True)
args = parser.parse_args()


###############################################################################
# Functions
###############################################################################

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
    reads_count = 0
    mi_count = 0

    # Read bam file and extract BX, chrom, start, and end
    for read in pysam.AlignmentFile(args.in_bam, 'rb').fetch():

        if read.has_tag('BX') and read.reference_name and read.reference_start and read.reference_end:
            if not old_chrom or read.reference_name == old_chrom:
                bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density = analyze_barcode(read, bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density)
            else:
                for mi in sorted(sorted(mi_chrom, key=mi_end.get), key=mi_start.get):
                    if mi_reads[mi] >= args.reads and (mi_end[mi]-mi_start[mi]) >= args.length:
                        out_bed.write(mi_chrom[mi]+'\t'+str(mi_start[mi])+'\t'+str(mi_end[mi])+'\t'+str(mi)+'\t'+str(float(mi_coverage[mi])/(mi_end[mi]-mi_start[mi]))+','+str(float(len(mi_density[mi]))/(mi_end[mi]-mi_start[mi]))+'\n')
                        mi_count += 1
                bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, temp_list, mi_coverage, mi_density = clear_variables()
                bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density = analyze_barcode(read, bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, mi_coverage, mi_density)
            reads_count += 1
    for mi in sorted(sorted(mi_chrom, key=mi_end.get), key=mi_start.get):
        if mi_reads[mi] >= args.reads and (mi_end[mi]-mi_start[mi]) >= args.length:
            out_bed.write(mi_chrom[mi]+'\t'+str(mi_start[mi])+'\t'+str(mi_end[mi])+'\t'+str(mi)+'\t'+str(float(mi_coverage[mi])/(mi_end[mi]-mi_start[mi]))+','+str(float(len(mi_density[mi]))/(mi_end[mi]-mi_start[mi]))+'\n')
            mi_count += 1
    bx_seen, bx_index, mi_chrom, mi_start, mi_end, mi_reads, old_chrom, temp_list, mi_coverage, mi_density = clear_variables()
    out_log.write('Extracted '+str(reads_count)+' reads from '+args.in_bam+'\nIdentified '+str(mi_count)+' 10X molecules with >='+str(args.reads)+' supporting reads and >='+str(args.length)+' length\n')


# Clear Matepair variables
def clear_matepair():
    first_ref = {}
    first_start = {}
    first_end = {}
    max_values = {}
    min_values = {}
    ref_name = {}
    return first_ref, first_start, first_end, max_values, min_values, ref_name


# Analyze Matepair bam file
def analyze_matepair_bam():

    # Read bam file and extract chrom, start, and end and write reads to bed
    count = 0
    old_ref = False
    for read in pysam.AlignmentFile(args.in_bam, 'rb').fetch():
        if not old_ref:
            first_ref, first_start, first_end, max_values, min_values, ref_name = clear_matepair()
        elif old_ref != read.reference_name:
            for query in sorted(sorted(ref_name, key=max_values.get), key=min_values.get):
                if (max_values[query] - min_values[query]) >= args.length:
                    out_bed.write(ref_name[query]+'\t'+str(min_values[query])+'\t'+str(max_values[query])+'\t'+query+'\n')
                    count += 1
            first_ref, first_start, first_end, max_values, min_values, ref_name = clear_matepair()
        old_ref = read.reference_name
        if read.query_name in first_ref:
            if read.reference_name == first_ref[read.query_name]:
                values = [first_start[read.query_name],first_end[read.query_name],read.reference_start,read.reference_end]
                max_values[read.query_name] = max(values)
                min_values[read.query_name] = min(values)
                ref_name[read.query_name] = read.reference_name
            first_ref.pop(read.query_name)
            first_start.pop(read.query_name)
            first_end.pop(read.query_name)
        else:
            first_ref[read.query_name] = read.reference_name
            first_start[read.query_name] = read.reference_start
            first_end[read.query_name] = read.reference_end
    for query in sorted(sorted(ref_name, key=max_values.get), key=min_values.get):
        if (max_values[query] - min_values[query]) >= args.length:
            out_bed.write(ref_name[query]+'\t'+str(min_values[query])+'\t'+str(max_values[query])+'\t'+query+'\n')
            count += 1
    out_log.write('Extracted '+str(count)+' reads >='+str(args.length)+' length from '+args.in_bam+'\n')


# Analyze Nanopore/PacBio bam file
def analyze_other_bam():

    # Read bam file and extract chrom, start, and end and write reads to bed
    count = 0
    for read in pysam.AlignmentFile(args.in_bam, 'rb').fetch():
        if (read.reference_end - read.reference_start) >= args.length:
            out_bed.write(read.reference_name+'\t'+str(read.reference_start)+'\t'+str(read.reference_end)+'\t'+read.query_name+'\n')
            count += 1
    out_log.write('Extracted '+str(count)+' reads >='+str(args.length)+' length from '+args.in_bam+'\n')


###############################################################################
# Run
###############################################################################

# Run functions
# ============================================================

out_bed = open(args.output+'.bed', 'w')
if args.log != 'sys.stderr':
    out_log = open(args.log, 'w')
else:
    out_log = sys.stderr
if args.type.upper() == '10X':
    analyze_10X_bam()
elif args.type.upper() == 'MATEPAIR':
    analyze_matepair_bam()
elif args.type.upper() == 'PACBIO' or args.type.upper() == 'NANOPORE':
    analyze_other_bam()
else:
    out_log.write('Unable to accept read type '+args.type+'. Please select 10X, Matepair, Nanopore, or PacBio.\n')
    sys.exit(2)
if args.log != sys.stderr:
    out_log.flush()
    out_log.close()
out_bed.flush()
out_bed.close()
pysam.tabix_index(args.output+'.bed', force=True, preset='bed')
