#!/usr/bin/env python

import os, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf> <BED/WIG>\n")
    sys.stderr.write("This script extracts the number of aligned species in each block of the input MAF file and outputs in either bed or wig format.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

header = False
species_ref = False
count = 0
total = 0

if sys.argv[2].upper() != 'BED' and sys.argv[2].upper() != 'WIG':
    sys.stderr.write('Unable to determine output format - please use either BED or WIG.\n')
    sys.exit(2)

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

for line in in_maf:
    if line.startswith('#'):
        header = True
    elif (header and not line.split()) or line.startswith('a'):
        header = False
        species_ref = False
        full_block = line
    elif line.startswith('s'):
        if not species_ref:
            line = line.rstrip().split()
            if sys.argv[2].upper() == 'BED':
                sys.stdout.write(line[1].split('.')[1]+'\t'+line[2]+'\t'+str(int(line[2])+int(line[3]))+'\t')
            elif sys.argv[2].upper() == 'WIG':
                sys.stdout.write('fixedStep chrom='+line[1].split('.')[1]+' start='+str(int(line[2])+1)+' step=1\n')
            total = int(line[3])
            species_ref = True
        count += 1
    elif not line.split():
        if sys.argv[2].upper() == 'BED':
            sys.stdout.write(str(count)+'\n')
        elif sys.argv[2].upper() == 'WIG':
            for n in range(0,total):
                sys.stdout.write(str(count)+'\n')
        count = 0
