#!/usr/bin/env python

import os, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <remove.scaffolds> <in.gff>\n")
    sys.stderr.write("This script removes scaffolds with annotations from the input GFF.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

file = sys.stdin
if sys.argv[1] != '-':
    file = open(sys.argv[1], 'r')

remove = {}
for line in file:
    remove[line.rstrip()] = True

for line in open(sys.argv[2],'r'):
    if line.startswith('#'):
        sys.stdout.write(line)
    elif line.rstrip().split('\t')[0] not in remove:
        sys.stdout.write(line)
