#!/usr/bin/env python

# Import packages
import argparse, os, re, sys
from Bio import SeqIO

# Parse arguments
parser = argparse.ArgumentParser(description='This script outputs a bed file of the gaps in an assembly.')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("in_fasta", help="input fasta file name", type=str)
args = parser.parse_args()

# Open FASTA, search for masked regions, print in BED3 format
with open(args.in_fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('[Nn]+', str(record.seq)):
            sys.stdout.write(str(record.id)+'\t'+str(match.start())+'\t'+str(match.end())+'\n')
