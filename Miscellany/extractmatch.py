#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4 or sys.argv[1] == "--help":
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <full.txt> <subset.txt> <column>\n")
    sys.stderr.write("Extracts line from full.txt if item in column number (base 0) is in subset.txt.\n")
    sys.exit(1)

keep = []
for line in open(sys.argv[2],'r'):
    keep.append(line.rstrip())

for line in open(sys.argv[1],'r'):
    if line.rstrip().split()[int(sys.argv[3])] in keep:
        sys.stdout.write(line)
