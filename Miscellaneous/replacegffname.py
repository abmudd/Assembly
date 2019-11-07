#!/usr/bin/env python

import os, sys

if len(sys.argv) != 5 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.gff> <names.txt> <old.name.column> <new.name.column>\n")
    sys.stderr.write("This script replaces old scaffold names with new names in the names.txt files.\nColumn names should be zero-based.\n")
    sys.exit(1)

names = {}

for item in open(sys.argv[2], 'r'):
    item = item.rstrip().split('\t')
    names[item[int(sys.argv[3])]] = item[int(sys.argv[4])]

for line in open(sys.argv[1], 'r'):
    if line.rstrip().split('\t')[0] in names:
        line = line.rstrip().split('\t')
        sys.stdout.write(names[line[0]]+'\t'+'\t'.join(line[1:])+'\n')
    else:
        sys.stdout.write(line)
