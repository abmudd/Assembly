#!/usr/bin/env python

import argparse, gzip, os, pysam, sys

parser = argparse.ArgumentParser(description='This script converts fastq comments to either /1 or /2.')
parser.add_argument("-g", "--gzip", action='store_true', help="gzip output file")
parser.add_argument("-o", "--out", metavar='STR', help="output file name [input.rehead or input.rehead.gz]", type=str)
parser.add_argument("-p", "--pair", metavar='1/2', help="set pair end flag for entire file with either 1 or 2", type=str)
required = parser.add_argument_group('required arguments')
required.add_argument("input", help="input fastq file name", type=str)
args = parser.parse_args()

in_fq = pysam.FastxFile(args.input)
if args.out and args.gzip:
    out_fq = gzip.open(args.out, 'w')
elif args.out:
    out_fq = open(args.out, 'w')
elif args.gzip:
    out_fq = gzip.open(args.input+'.rehead.gz', 'w')
else:
    out_fq = open(args.input+'.rehead', 'w')

for record in in_fq:
    if args.pair:
        out_fq.write('@'+str(record.name)+'/'+args.pair+'\n'+str(record.sequence)+'\n+\n'+str(record.quality)+'\n')
    elif record.comment:
        record_pair = record.comment.split(':')[0]
        if record_pair == '1' or record_pair == '2':
            out_fq.write('@'+str(record.name)+'/'+record_pair+'\n'+str(record.sequence)+'\n+\n'+str(record.quality)+'\n')
        else:
            sys.stderr.write('Unable to determine pair from fastq header comment.\n')
            sys.exit()
    else:
        sys.stderr.write('Unable to find fastq header comment and pair end flag is not set.\n')
        sys.exit()

in_fq.close()
out_fq.close()
