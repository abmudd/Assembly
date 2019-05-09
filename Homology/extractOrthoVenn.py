#!/usr/bin/env python

import os, pysam, sys

if len(sys.argv) != 6 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.orthologs> <species1> <species1.pro.aa> <species2> <species2.pro.aa>\n")
    sys.stderr.write("This script extracts amino acid sequences from OrthoVenn 1-to-1 orthologs.\n")
    sys.exit(1)

in_ortho = sys.stdin
if sys.argv[1] != "-":
    in_ortho = open(sys.argv[1],'r')

species1_fai = {}
with pysam.FastxFile(sys.argv[3]) as fg:
    for entry in fg:
        species1_fai[entry.name] = entry.sequence

species2_fai = {}
with pysam.FastxFile(sys.argv[5]) as fh:
    for entry in fh:
        species2_fai[entry.name] = entry.sequence

for line in in_ortho:
    species1 = False
    species2 = False
    for item in line.rstrip().split('\t')[2].split(' '):
        if item.startswith(sys.argv[2]):
            species1 = item
        elif item.startswith(sys.argv[4]):
            species2 = item
    if species1 and species2:
        species1_split = species1.split('|')[1]
        species2_split = species2.split('|')[1]
        sys.stdout.write('>'+species1+'\n'+species1_fai[species1_split]+'\n')
        sys.stdout.write('>'+species2+'\n'+species2_fai[species2_split]+'\n')
#        sys.stderr.write((species1_split+'\t'+species2_split).replace('N','-').replace('P','+').replace('_','\t')+'\n')
