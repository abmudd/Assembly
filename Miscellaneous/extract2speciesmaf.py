#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf> <species1> <species2>\n")
    sys.stderr.write("This script filters maf alignments where each block contains exactly one species 1 alignment and one species 2 alignment and requires the alignment to have at least 1 base of overlap. The output is a new pairwise MAF file.\n")
    sys.exit(1)

header = False
species1 = False
species1_count = 0
species2 = False
species2_count = 0

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

for line in in_maf:
    if line.startswith('#'):
        header = True
        sys.stdout.write(line)
    elif header and not line.split():
        header = False
        sys.stdout.write(line)
    elif line.startswith('a'):
        header = False
    elif line.startswith('s'):
        line_split = line.rstrip().split()
        if line_split[1].startswith(sys.argv[2]):
            species1 = line_split
            species1_count += 1
        elif line_split[1].startswith(sys.argv[3]):
            species2 = line_split
            species2_count += 1
    elif not line.split():
        align = False
        seq1 = ''
        seq2 = ''
        if species1_count == 1 and species2_count == 1 and len(species1[6]) > 0 and len(species2[6]) > 0 and len(species1[6]) == len(species2[6]):
            for i in range(0,len(species1[6])):
                if species1[6][i] != '-' and species1[6][i].upper() != 'N' and species2[6][i] != '-' and species2[6][i].upper() != 'N':
                    align = True
                if species1[6][i] != '-' or species2[6][i] != '-':
                   seq1 += species1[6][i]
                   seq2 += species2[6][i]
        if align:
            sys.stdout.write('a\n'+' '.join(species1[0:6])+' '+seq1+'\n'+' '.join(species2[0:6])+' '+seq2+'\n\n')
        species1_count = 0
        species2_count = 0
