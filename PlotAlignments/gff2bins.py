#!/usr/bin/env python

import argparse, math, os, sys

parser = argparse.ArgumentParser(description='This script bins repeats based on the window.')
parser.add_argument("-b", "--binsize", metavar='INT', help="bin size [100000]", type=int, default=100000)
parser.add_argument("-o", "--output", metavar='STR', help="output file [stdout]", type=str, default='stdout')
required = parser.add_argument_group('required arguments')
required.add_argument("fai", help="input concat fai", type=str)
required.add_argument("gff", help="input concat repeat gff", type=str)
args = parser.parse_args()

if args.output == 'stdout':
    output = sys.stdout
else:
    output = open(args.output,'w')

bins = {}
for line in open(args.fai,'r'):
    line = line.rstrip().split('\t')
    for n in range(0,int(math.ceil(float(line[1])/args.binsize))):
        bins[line[0]+':'+str(n)] = 0

for line in open(args.gff,'r'):
    if not line.startswith('#'):
        line = line.rstrip().split('\t')
        bins[line[0]+':'+str(int(math.floor(float(line[3])/args.binsize)))] += 1

for item in bins:
    chrom = item.split(':')[0]
    size = str(int(item.split(':')[1])*args.binsize)
    output.write(chrom+'\t'+size+'\t'+str(bins[item])+'\n')
