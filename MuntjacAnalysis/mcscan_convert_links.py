#!/usr/bin/env python

import os, sys

if len(sys.argv) < 3 or len(sys.argv) > 6 or sys.argv[1] == "--help":
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.circos_links> <out.prefix> <count> <first> <second>\n")
    sys.stderr.write("This script converts the Circos links into JCVI files. Count is 1, first is A, and second is B if not provided.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

outlayout = open(sys.argv[2]+'.layout','w')
outlayout.write('# y, xstart, xend, rotation, color, label, va,  bed\n')
outlayout.write(' .6,     .1,    .8,       0,      ,      , top, '+sys.argv[2]+'_1.bed\n')
outlayout.write(' .4,     .1,    .8,       0,      ,      , top, '+sys.argv[2]+'_2.bed\n')
outlayout.write('# edges\n')
outlayout.write('e, 0, 1, '+sys.argv[2]+'.anchors\n')
outlayout.close()

outbed1 = open(sys.argv[2]+'_1.bed','w')
outbed2 = open(sys.argv[2]+'_2.bed','w')
outanchors = open(sys.argv[2]+'.anchors','w')
count = 1
first = 'A'
second = 'B'
if len(sys.argv) >= 4:
    count = int(sys.argv[3])
if len(sys.argv) >= 5:
    first = sys.argv[4]
if len(sys.argv) == 6:
    second = sys.argv[5]
chr1 = []
chr2 = []

for line in open(sys.argv[1],'r'):
    line = line.rstrip().split(' ')
    chrA = line[0].split('.')[1]
    chrB = line[3].split('.')[1]
    outbed1.write(chrA+'\t'+line[1]+'\t'+str(int(line[1])+3)+'\t'+first+str(count)+'\t0\t+\n')
    outbed1.write(chrA+'\t'+str(int(line[2])-3)+'\t'+line[2]+'\t'+first+str(count+1)+'\t0\t+\n')
    outbed2.write(chrB+'\t'+line[4]+'\t'+str(int(line[4])+3)+'\t'+second+str(count)+'\t0\t+\n')
    outbed2.write(chrB+'\t'+str(int(line[5])-3)+'\t'+line[5]+'\t'+second+str(count+1)+'\t0\t+\n')
    outanchors.write(first+str(count)+'\t'+first+str(count+1)+'\t'+second+str(count)+'\t'+second+str(count+1)+'\t500\t+\n')
    if chrA not in chr1:
        chr1.append(chrA)
    if chrB not in chr2:
        chr2.append(chrB)
    count += 2
outbed1.close()
outbed2.close()
outanchors.close()

outseqid = open(sys.argv[2]+'.seqid','w')
outseqid.write(','.join(chr1)+'\n')
outseqid.write(','.join(chr2)+'\n')
outseqid.close()
