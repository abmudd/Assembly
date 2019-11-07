#!/usr/bin/env python

import os, pysam, sys
from collections import defaultdict

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: %s <halfdepth.jellyfish.blast.out> <in.fa>\n" % os.path.basename(sys.argv[0]))
    sys.stderr.write("This script filters half-depth contigs and outputs the contigs to remove.\n")
    sys.exit(1)

file = sys.stdin
if sys.argv[1] != '-':
    file = open(sys.argv[1], 'r')
blast_count = {}
query_old = False
blast_col = []
count = 0

for blast in file:
    split_blast = blast.rstrip().split('\t')
    if not query_old or split_blast[0] != query_old:
        if count == 2:
            ref0 = blast_col[0].split('\t')[1]
            ref1 = blast_col[1].split('\t')[1]
            if ref0+':'+ref1 in blast_count:
                blast_count[ref0+':'+ref1] += 1
            elif ref1+':'+ref0 in blast_count:
                blast_count[ref1+':'+ref0] += 1
            else:
                blast_count[ref0+':'+ref1] = 1
        count = 0
        blast_col = []
    blast_col.append(blast.rstrip())
    query_old = split_blast[0]
    count += 1

reduce_blast_count = defaultdict(list)
count = []

for i in blast_count:
    reduce_blast_count[blast_count[i]].append(i)
    if blast_count[i] not in count:
        count.append(blast_count[i])

blast_count = {}
seen = {}
fasta = pysam.FastaFile(sys.argv[2])
length = dict(zip(fasta.references, fasta.lengths))

for value in sorted(count, reverse=True):
    for key in reduce_blast_count[value]:
        key = key.split(':')
        if key[0] not in seen and key[1] not in seen:
            length0 = length[key[0]]
            length1 = length[key[1]]
            if length0 < length1:
                sys.stdout.write(key[0]+'\n')
            else:
                sys.stdout.write(key[1]+'\n')
            seen[key[0]] = True
            seen[key[1]] = True
