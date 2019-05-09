#!/usr/bin/env python

import os, pysam, sys

if len(sys.argv) != 4 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <blast.out> <fasta> <approved.seqid>\n")
    sys.stderr.write("This script parses nt blast output to identify contigs without hits from the approved sequence id list.\n")
    sys.exit(1)

# Add approved sequence ids to approval_list
approval_list = {}
for field in open(sys.argv[3], 'r'):
    approval_list[field.rstrip()] = True

# Parse blast output and write out scaffolds not in approval_list
seen_scaffolds = {}
sys.stdout.write('# List of scaffolds with best hit not on the approved list\n')
for line in open(sys.argv[1], 'r'):
    line_split = line.rstrip().split('\t')
    if ' ' in line_split[1]:
        line_split[1] = line_split[1].split(' ')[0]
    if line_split[1] not in approval_list and line_split[0] not in seen_scaffolds:
        sys.stdout.write(line_split[0]+'\n')
        sys.stdout.write('# '+line)
        seen_scaffolds[line_split[0]] = True

# Parse fasta and check for any scaffolds without a hit
#fasta = pysam.FastaFile(sys.argv[2])
#for scaffold in zip(fasta.references):
#    query_scaffold = str(scaffold).split('\'')[1]
#    if query_scaffold not in approved_scaffolds and query_scaffold not in not_approved_scaffolds:
#        sys.stdout.write(query_scaffold+'\n')
