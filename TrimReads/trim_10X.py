#!/usr/bin/env python

import argparse, os, pysam, sys
from itertools import izip

# Parse arguments
parser = argparse.ArgumentParser(description='This script adapter trims and outputs fastq files in the expected format of programs scaff10x and break10x.')
parser.add_argument("-r1", "--out_read1", metavar='STR', help="output fastq read 1 file name [read-BC_1.fastq]", type=str, default='read-BC_1.fastq')
parser.add_argument("-r2", "--out_read2", metavar='STR', help="output fastq read 2 file name [read-BC_2.fastq]", type=str, default='read-BC_2.fastq')
parser.add_argument("-n", "--out_name", metavar='STR', help="output name file name [read-BC_1.name]", type=str, default='read-BC_1.name')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("in_read1", help="input fastq read 1 file name", type=str)
required.add_argument("in_read2", help="input fastq read 2 file name", type=str)
args = parser.parse_args()

out_fq1 = open(args.out_read1,'w')
out_name = open(args.out_name,'w')
out_fq2 = open(args.out_read2,'w')
count = 0

for record1,record2 in izip(pysam.FastxFile(args.in_read1),pysam.FastxFile(args.in_read2)):
    out_fq1.write('@'+record1.name+'_'+record1.sequence[0:16]+'\n'+record1.sequence[23:]+'\n+\n'+record1.quality[23:]+'\n')
    out_fq2.write('@'+record2.name+'_'+record1.sequence[0:16]+'\n'+record2.sequence+'\n+\n'+record2.quality+'\n')
    out_name.write(str(count)+' '+record1.name+'_'+record1.sequence[0:16]+'\n')
    count += 1

out_fq1.close()
out_name.close()
out_fq2.close()
