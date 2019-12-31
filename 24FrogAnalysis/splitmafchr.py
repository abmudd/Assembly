#!/usr/bin/env python

import os, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf> <out.prefix>\n")
    sys.stderr.write("This script splits the input MAF file into unique MAF files for each reference scaffold and chromosome.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

header = False
species_ref = False
species_qry = False
full_block = False
seen_list = []

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

for line in in_maf:
    if line.startswith('#'):
        header = True
    elif (header and not line.split()) or line.startswith('a'):
        header = False
        species_qry = False
        full_block = line
    elif line.startswith('s'):
        if not species_qry:
            if not species_ref:
                species_ref = line.rstrip().split('.')[0]
            elif species_ref != line.rstrip().split('.')[0]:
                sys.err.write('Ref species is not consistent ... please see line starting with '+line.rstrip().split('.')[0]+'\n')
                sys.exit(2)
            species_qry = line.rstrip().split()[1].split('.')[1]
        full_block += line
    elif not line.split():
        if species_qry not in seen_list:
            seen_list.append(species_qry)
            output = open(sys.argv[2]+'.'+species_qry+'.maf','w')
            output.write('##maf version=1\n\n')
        else:
            output = open(sys.argv[2]+'.'+species_qry+'.maf','a')
        output.write(full_block+'\n')
        output.close()
