#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("This script identifies runs of co-linearity from the output of maf2stats.py. User must run script with maf2stats.py output as: tail -n +2 MAF2STATS_OUT | sort -k5,5 -k7,7n | awk '{print $0 \"\\t\" NR}' | sort -k1,1 -k3,3n | "+os.path.basename(sys.argv[0])+" - <species1> <species2>.\n")
    sys.exit(1)

species1_chr = False
species1_strand = False
species1_start = False
species1_end = False
species2_chr = False
species2_strand = False
species2_start = False
species2_end = False
alnsize = False
aln100id = False
species2_order = False

for line in sys.stdin:
    line = line.rstrip().split()
    src1 = line[0].split('.')
    src2 = line[4].split('.')
    if src1[0] != sys.argv[2] or src2[0] != sys.argv[3]:
        sys.stderr.write('ERROR: Line is incorrectly ordered or species not found: '+' '.join(line)+'\n')
        sys.exit(2)
    if not species2_order:
        species1_chr = src1[1]
        species1_strand = line[1]
        species1_start = line[2]
        species1_end = line[3]
        species2_chr = src2[1]
        species2_strand = line[5]
        species2_start = line[6]
        species2_end = line[7]
        alnsize = int(line[8])
        aln100id = int(line[9])
        species2_order = int(line[10])
    elif src1[1] == species1_chr and src2[1] == species2_chr and line[1] == species1_strand and line[5] == species2_strand and abs(int(line[10]) - species2_order) < 2:
        species1_start = str(min(int(species1_start),int(line[2])))
        species1_end = str(max(int(species1_end),int(line[3])))
        species2_start = str(min(int(species2_start),int(line[6])))
        species2_end = str(max(int(species2_end),int(line[7])))
        alnsize += int(line[8])
        aln100id += int(line[9])
        species2_order = int(line[10])
    else:
        sys.stdout.write(sys.argv[2]+'.'+species1_chr+'\t'+species1_strand+'\t'+species1_start+'\t'+species1_end+'\t'+sys.argv[3]+'.'+species2_chr+'\t'+species2_strand+'\t'+species2_start+'\t'+species2_end+'\t'+str(alnsize)+'\t'+str(aln100id)+'\n')
        species1_chr = src1[1]
        species1_strand = line[1]
        species1_start = line[2]
        species1_end = line[3]
        species2_chr = src2[1]
        species2_strand = line[5]
        species2_start = line[6]
        species2_end = line[7]
        alnsize = int(line[8])
        aln100id = int(line[9])
        species2_order = int(line[10])
sys.stdout.write(sys.argv[2]+'.'+species1_chr+'\t'+species1_strand+'\t'+species1_start+'\t'+species1_end+'\t'+sys.argv[3]+'.'+species2_chr+'\t'+species2_strand+'\t'+species2_start+'\t'+species2_end+'\t'+str(alnsize)+'\t'+str(aln100id)+'\n')
