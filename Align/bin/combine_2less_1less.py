#!/usr/bin/env python

import os, sys

if len(sys.argv) != 5 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <2orless.bed> <1orless.bed> <in.fasta> <out.prefix>\n")
    sys.stderr.write("This script combines 2orless.bed and 1orless.bed and extracts regions with 2 or less at the ends of scaffolds and 1 or less in between. The input fasta must have a fai file.\n")
    sys.exit(1)

length = {}
for line in open(sys.argv[3]+'.fai', 'r'):
    line = line.rstrip().split('\t')
    length[line[0]] = line[1]

out_trim = open(sys.argv[4]+'.trim.bed', 'w')
out_mask = open(sys.argv[4]+'.mask.bed', 'w')

for line in open(sys.argv[1], 'r'):
    sline = line.rstrip().split('\t')
    if sline[1] == '0':
        out_trim.write(line)
    elif sline[2] == length[sline[0]]:
        out_trim.write(line)

for line in open(sys.argv[2], 'r'):
    sline = line.rstrip().split('\t')
    if sline[1] != '0' and sline[2] != length[sline[0]]:
        out_mask.write(line)

out_trim.close()
out_mask.close()
