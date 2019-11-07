#!/usr/bin/env python

import os, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: %s <in.config> <out.config>\n" % os.path.basename(sys.argv[0]))
    sys.stderr.write("This script converts the input config to an SSPACE config file.\n")
    sys.exit(1)

file = sys.stdin
if sys.argv[1] != '-':
    file = open(sys.argv[1], 'r')
output = open(sys.argv[2], 'w')

orient = {}
insave = {}
insstd = {}
read1  = {}
read2  = {}

for line in file:
    line = line.rstrip().split('\t')
    read1[line[0]]  = line[1]
    read2[line[0]]  = line[2]
    orient[line[0]] = line[3]
    insave[line[0]] = int(line[4])
    insstd[line[0]] = int(line[5])

count = 1
old = False

for key, value in sorted(insave.iteritems(), key=lambda (k,v): (v,k)):
    if not old:
        old = key
    else:
        if abs(insave[old] - value) > insstd[old]:
            count += 1
            old = key
    output.write('Lib'+str(count)+' bwa '+read1[key]+' '+read2[key]+' '+str(value)+' '+str(round((float(insstd[key])/value), 2))+' '+orient[key]+'\n')

if sys.argv[1] != '-':
    file.close()
output.close()
