#!/usr/bin/env python

import os, sys

if len(sys.argv) != 2 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.genePred>\n")
    sys.stderr.write("This script extracts the largest transcript for each gene from the input genePred.\n")
    sys.exit(1)

seen_size = {}
seen_output = {}
for line in open(sys.argv[1],'r'):
    line_split = line.rstrip().split('\t')
    if line_split[11] in seen_size:
        if int(line_split[4])-int(line_split[3]) > seen_size[line_split[11]]:
            seen_size[line_split[11]] = int(line_split[4])-int(line_split[3])
            seen_output[line_split[11]] = line
    else:
        seen_size[line_split[11]] = int(line_split[4])-int(line_split[3])
        seen_output[line_split[11]] = line

for item in seen_output:
    sys.stdout.write(seen_output[item])
