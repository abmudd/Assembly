#!/usr/bin/env python

import argparse, os, pysam, re, sys

parser = argparse.ArgumentParser(description='This script removes N at start and end of each scaffold and removes scaffolds only containing N.')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("in_fasta", help="input fasta file name", type=str)
args = parser.parse_args()

for record in pysam.FastxFile(args.in_fasta):
    start = False
    startN = False
    end = False
    endN = False
    if record.sequence[0] == 'N' or record.sequence[0] == 'n':
        startN = True
    if record.sequence[len(record.sequence)-1] == 'N' or record.sequence[len(record.sequence)-1] == 'n':
        endN = True
    if startN and endN:
        if bool(re.search('^[Nn]+$',record.sequence)):
            continue
        else:
            for n in range(0,len(record.sequence)):
                if not record.sequence[n] in ['N','n']:
                    start = n
                    break
            for n in range(len(record.sequence)-1,-1,-1):
                if not record.sequence[n] in ['N','n']:
                    end = n+1
                    break
            sys.stdout.write('>'+record.name+'\n')
            sys.stdout.write(record.sequence[start:end]+'\n')
    elif startN:
        for n in range(0,len(record.sequence)):
            if not record.sequence[n] in ['N','n']:
                start = n
                break
        sys.stdout.write('>'+record.name+'\n')
        sys.stdout.write(record.sequence[start:]+'\n')
    elif endN:
        for n in range(len(record.sequence)-1,-1,-1):
            if not record.sequence[n] in ['N','n']:
                end = n+1
                break
        sys.stdout.write('>'+record.name+'\n')
        sys.stdout.write(record.sequence[:end]+'\n')
    else:
        sys.stdout.write('>'+record.name+'\n')
        sys.stdout.write(record.sequence+'\n')
