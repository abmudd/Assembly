#!/usr/bin/env python

import os, sys

if len(sys.argv) != 2 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf>\n")
    sys.stderr.write("This script extracts indels relative to reference species in maf alignment.\n")
    sys.exit(1)

species1 = False
species2 = False
header = False

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

sys.stdout.write('#Indel\tSize\n')

for line in in_maf:
    if line.startswith('#'):
        header = True
    elif header and not line.split():
        header = False
    elif line.startswith('a'):
        header = False
    elif line.startswith('s'):
        line_split = line.rstrip().split()
        if not species1:
            species1 = line_split
        else:
            species2 = line_split
    else:
        insertion = 0
        deletion = 0
        for i in range(0,len(species1[6])):
            out = False
            if species1[6][i] == '-' and species2[6][i] == '-':
                continue
            elif species1[6][i] == '-':
                deletion += 1
            elif species2[6][i] == '-':
                insertion += 1
            else:
                if deletion > 0:
                    sys.stdout.write('DEL\t'+str(deletion)+'\n')
                    deletion = 0
                    out = True
                elif insertion > 0:
                    sys.stdout.write('INS\t'+str(insertion)+'\n')
                    insertion = 0
                    out = True
            if i == len(species1[6])-1 and not out:
                if deletion > 0:
                    sys.stdout.write('DEL\t'+str(deletion)+'\n')
                elif insertion > 0:
                    sys.stdout.write('INS\t'+str(insertion)+'\n')
        species1 = False
        species2 = False
