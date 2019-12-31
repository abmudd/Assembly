#!/usr/bin/env python

import os, sys

if len(sys.argv) != 5 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf> <species1> <species2> <MAF/BED>\n")
    sys.stderr.write("This script filters maf alignments where each block contains exactly one species 1 alignment and one species 2 alignment and requires the alignment to have at least 1 base of overlap. The output is either a new pairwise MAF file or a BEDPE file.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

header = False
species1 = False
species1_count = 0
species2 = False
species2_count = 0

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

if sys.argv[4].upper() != "MAF" and sys.argv[4].upper() != "BED":
    sys.stderr.write('Unable to determine output format - please use either MAF or BED.\n')
    sys.exit(2)

for line in in_maf:
    if line.startswith('#'):
        header = True
        if sys.argv[4].upper() == "MAF":
            sys.stdout.write(line)
    elif header and not line.split():
        header = False
        if sys.argv[4].upper() == "MAF":
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
            if sys.argv[4].upper() == "MAF":
                sys.stdout.write('a\n'+' '.join(species1[0:6])+' '+seq1+'\n'+' '.join(species2[0:6])+' '+seq2+'\n\n')
            elif sys.argv[4].upper() == "BED":
                if species2[4] != "-":
                    sys.stdout.write(species1[1]+'\t'+species1[2]+'\t'+str(int(species1[2])+int(species1[3]))+'\t'+species2[1]+'\t'+species2[2]+'\t'+str(int(species2[2])+int(species2[3]))+'\n')
                else:
                    sys.stdout.write(species1[1]+'\t'+species1[2]+'\t'+str(int(species1[2])+int(species1[3]))+'\t'+species2[1]+'\t'+str(int(species2[5])-(int(species2[2])+int(species2[3]))-1)+'\t'+str(int(species2[5])-int(species2[2])-1)+'\n')
        species1_count = 0
        species2_count = 0
