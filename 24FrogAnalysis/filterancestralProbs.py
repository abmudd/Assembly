#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.ancestralProbs> <node> <min_p>\n")
    sys.stderr.write("This script extracts bases from input RAxML ancestralProbs of a particular node with a minimum probability and an unambiguous base.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

if sys.argv[1] == '-':
    in_file = sys.stdin
else:
    in_file = open(sys.argv[1],'r')
min_p = float(sys.argv[3])
for line in in_file:
    line = line.rstrip().split('\t')
    if line[0] == sys.argv[2]:
        if line[2] in ['A','C','G','T','a','c','g','t'] and (float(line[3]) >= min_p or float(line[4]) >= min_p or float(line[5]) >= min_p or float(line[6]) >= min_p):
            sys.stdout.write(line[2]+'\n')
        else:
            sys.stdout.write('-\n')
