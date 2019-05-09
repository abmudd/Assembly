#!/usr/bin/env python

import os, pysam, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <keep.scaffolds> <in.fa>\n")
    sys.stderr.write("This script keeps scaffolds from the input fasta.\n")
    sys.exit(1)

file = sys.stdin
if sys.argv[1] != '-':
    file = open(sys.argv[1], 'r')

keep = {}
for line in file:
    keep[line.rstrip()] = True

with pysam.FastxFile(sys.argv[2]) as fh:
    for entry in fh:
        if entry.name in keep:
            print('>'+entry.name)
            print(entry.sequence)
