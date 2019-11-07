#!/usr/bin/env python

import os, pysam, subprocess, sys, tempfile
from distutils.spawn import find_executable

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <blast.output> <in.fa>\n")
    sys.stderr.write("This script parses blast output and identifies contigs for removal.\n")
    sys.exit(1)

if find_executable('bedtools') == None:
    sys.stderr.write("Error: bedtools could not be located in the environment\n")
    sys.exit(1)

if os.stat(sys.argv[1]).st_size == 0:
    sys.stderr.write("Error: blast output "+sys.argv[1]+" is empty\n")
    sys.exit(1)

fasta = pysam.FastaFile(sys.argv[2])
contig_size = {}

for scaffold, length in zip(fasta.references,fasta.lengths):
    contig_size[scaffold] = length

blast_sort = tempfile.NamedTemporaryFile(delete=False)
subprocess.check_call('sort -k1,1 -k2,2 -k7,7n '+sys.argv[1]+' >'+blast_sort.name, shell=True)

bedtools_in = tempfile.NamedTemporaryFile(delete=False)
bedtools_out = tempfile.NamedTemporaryFile(delete=False)

for line in open(blast_sort.name, 'r'):
    line = line.rstrip().split('\t')
    if ':' in line[0]:
        temp = line[0].split(':')
        if len(temp) > 2:
            sys.stderr.write('Warning: Multiple ":" in name '+line[0]+'\n')
            line[0] = ':'.join(temp[:-1])
        else:
            line[0] = temp[0]
    if '+' in line[0] or '+' in line[1]:
        sys.stderr.write('Error: Multiple "+" in name '+line[0]+' or '+line[1]+'\n')
        sys.exit(1)
    bedtools_in.write(line[0]+'+'+line[1]+'\t'+line[6]+'\t'+line[7]+'\n')

bedtools_in.flush()
subprocess.check_call('bedtools merge -d 5 -i ' + bedtools_in.name, stdout=bedtools_out, shell=True)
bedtools_out.flush()

first_line = True
old_query = False
bases = 0
range = []
check_per_list = []
check_per_dict = {}

for line in open(bedtools_out.name, 'r'):
    line = line.rstrip().split('\t')
    query = line[0].split('+')
    if first_line or query == old_query:
        bases += (int(line[2]) - int(line[1])) + 1
        range.append(line[1]+'-'+line[2])
        first_line = False
    else:
        percentage = (bases/float(contig_size[old_query[0]]))*100
        if bases > 200 or percentage > 50:
            sys.stdout.write(old_query[0]+'\t'+old_query[1]+'\t'+str(bases)+'\t'+str(percentage)+'\t'+','.join(range)+'\t'+str(contig_size[old_query[0]])+'\n')
            if old_query[0] not in check_per_list or percentage > check_per_dict[old_query[0]]:
                check_per_dict[old_query[0]] = percentage
                check_per_list.append(old_query[0])
        bases = (int(line[2]) - int(line[1])) + 1
        range = [line[1]+'-'+line[2]]
    old_query = query

percentage = (bases/float(contig_size[old_query[0]]))*100
if bases > 200 or percentage > 50:
    sys.stdout.write(old_query[0]+'\t'+old_query[1]+'\t'+str(bases)+'\t'+str(percentage)+'\t'+','.join(range)+'\t'+str(contig_size[old_query[0]])+'\n')
    if old_query[0] not in check_per_list or percentage > check_per_dict[old_query[0]]:
        check_per_dict[old_query[0]] = percentage
        check_per_list.append(old_query[0])
        
item_remove = []
item_check = []

for item in check_per_list:
    if check_per_dict[item] >= 98 and item not in item_remove:
        item_remove.append(item)
        if item in item_check:
            item_check.remove(item)
    elif item not in item_remove and item not in item_check:
        item_check.append(item)

if item_remove or item_check:
    sys.stdout.write('\n')

if item_remove:
    sys.stdout.write('\n'+'Scaffolds to remove: '+','.join(item_remove))

if item_check:
    sys.stdout.write('\n'+'Scaffolds to check: '+','.join(item_check))

if item_remove or item_check:
    sys.stdout.write('\n')

subprocess.check_call(['rm', bedtools_in.name])
subprocess.check_call(['rm', bedtools_out.name])
subprocess.check_call(['rm', blast_sort.name])
