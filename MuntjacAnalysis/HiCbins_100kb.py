#!/usr/bin/env python

import math, os, sys

if len(sys.argv) != 2 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.mnd>\n")
    sys.stderr.write("This script extracts number of contacts per 100 kb bin for the M. muntjak assembly from Juicer merged no dups output.\n")
    sys.exit(1)

bins = []
for n in range(0,11078):
    bins.append('Chr1.'+str(n+1))
for n in range(0,6825):
    bins.append('Chr2.'+str(n+1))
for n in range(0,6577):
    bins.append('Chr3.'+str(n+1))

bin_matrix = {}
for a in bins:
    for b in bins:
        bin_matrix[a+'-'+b] = 0

readline = sys.stdin
if sys.argv[1] != '-':
    readline = open(sys.argv[1],'r')

for line in readline:
    line = line.rstrip().split(' ')
    bin_matrix[line[1]+'.'+str(int(math.ceil(int(line[2])/100000.0)))+'-'+line[5]+'.'+str(int(math.ceil(int(line[6])/100000.0)))] += 1
    if line[1]+'.'+str(int(math.ceil(int(line[2])/100000.0))) != line[5]+'.'+str(int(math.ceil(int(line[6])/100000.0))):
        bin_matrix[line[5]+'.'+str(int(math.ceil(int(line[6])/100000.0)))+'-'+line[1]+'.'+str(int(math.ceil(int(line[2])/100000.0)))] += 1

sys.stdout.write('\t'+'\t'.join(bins)+'\n')
for a in bins:
    temp = []
    for b in bins:
        temp.append(str(bin_matrix[a+'-'+b]))
    sys.stdout.write(a+'\t'+'\t'.join(temp)+'\n')
