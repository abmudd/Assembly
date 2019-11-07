#!/usr/bin/env python

# This script is derived from a script written by Jessen Bredeson <jessenbredeson@berkeley.edu>.

###############################################################################
# Setup
###############################################################################

# Import modules
# ============================================================

import argparse, numpy, os, pysam, sys


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script outputs the read length, median coverage, and average coverage for each contig.')
parser.add_argument("-b","--in_bed", metavar='STR', help="only reads overlapping this bed", type=str, default=False)
required = parser.add_argument_group('required arguments')
required.add_argument("in_bam", metavar='in_bam', help="sorted and indexed input bam", type=str)
required.add_argument("in_fa", metavar='in_fa', help="input fasta", type=str)
args = parser.parse_args()


# Set variables
# ============================================================

bam = pysam.AlignmentFile(args.in_bam)
fasta = pysam.FastaFile(args.in_fa)


###############################################################################
# Run
###############################################################################

# Import fasta file
# ============================================================

length = dict(zip(fasta.references, fasta.lengths))


# Import bed file
# ============================================================

bed_list = {}
if args.in_bed:
    for line in open(args.in_bed, 'r'):
        bed_list[line.split('\t')[0]] = True


# Parse contigs
# ============================================================

for contig in fasta.references:
    if not args.in_bed or contig in bed_list:
        median_depth = 0
        average_depth = 0
        counts = {}
        total = 0
        sum = 0
        for pile in bam.pileup(contig, stepper='nofilter'):
            if pile.nsegments in counts:
                counts[pile.nsegments] += 1
            else:
                counts[pile.nsegments] = 1

            sum += pile.nsegments
            total += 1
            
        if total > 0:
            average_depth = sum / total
        midpoint = total / 2.0
        total = 0
        for depth in sorted(counts):
            if counts[depth]+total < midpoint:
                total += counts[depth]
            else:
                median_depth = depth
                break

        sys.stdout.write("%s\t%d\t%d\t%d\n" % (contig, length[contig], median_depth, average_depth))
