#!/usr/bin/env python

import os, sys

if len(sys.argv) != 5 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.fasta> <names.txt> <old.name.column> <new.name.column>\n")
    sys.stderr.write("This script replaces old scaffold names with new names in the names.txt files.\nColumn names should be zero-based.\n")
    sys.exit(1)

names = {}

for item in open(sys.argv[2], 'r'):
    item = item.rstrip().split('\t')
    names[item[int(sys.argv[3])]] = item[int(sys.argv[4])]

for line in open(sys.argv[1], 'r'):
    if line.startswith('>'):
        query = line.rstrip().replace('>','').split(' ')[0]
        if query in names:
            sys.stdout.write('>'+str(names[query])+'\n')
        else:
            sys.stdout.write(line)
    else:
        sys.stdout.write(line)
