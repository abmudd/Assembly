#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4 or sys.argv[1] == "--help":
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.bed> <in.fai> <chr_to_invert>\n")
    sys.stderr.write("This script inverts the bed locations for a particular chromosome.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

fai = {}
for line in open(sys.argv[2]):
    line = line.rstrip().split()
    fai[line[0]] = int(line[1])

for line in open(sys.argv[1]):
    line_split = line.rstrip().split()
    if line_split[0] == sys.argv[3]:
        sys.stdout.write(line_split[0]+'\t'+str(fai[line_split[0]]-int(line_split[2]))+'\t'+str(fai[line_split[0]]-int(line_split[1]))+'\t'+line_split[3]+'\t0\t+\n')
    else:
        sys.stdout.write(line)
