#!/usr/bin/env python

import os, sys

if len(sys.argv) != 6 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf> <species1> <species2> <UCElength> <UCE%id>\n")
    sys.stderr.write("This script calculates the locations of UCEs and outputs the location as a bed file relative to species1. The script requires 1-to-1 orthologous alignment of the two species.\n")
    sys.stderr.write("Version: 1.0\n")
    sys.exit(1)

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

def similar(w1, w2):
    w1 = w1 + ' ' * (len(w2) - len(w1))
    w2 = w2 + ' ' * (len(w1) - len(w2))
    return sum(1 if i == j else 0 for i, j in zip(w1, w2))

def clear():
    species1 = False
    species1_info = False
    species1_count = 0
    species2 = False
    species2_info = False
    species2_count = 0
    return species1,species1_info,species1_count,species2,species2_info,species2_count

header = False
species1,species1_info,species1_count,species2,species2_info,species2_count = clear()

for line in in_maf:
    if line.startswith('#'):
        header = True
    elif (header and not line.split()) or line.startswith('a'):
        header = False
    elif line.startswith('s'):
        line_split = line.rstrip().split()
        if line_split[1].startswith(sys.argv[2]):
            species1 = line_split[6]
            species1_count += 1
            if line_split[4] == '+':
                species1_info = [line_split[1],line_split[4],line_split[2],str(int(line_split[2])+int(line_split[3]))]
            else:
                species1_info = [line_split[1],line_split[4],str(int(line_split[5])-int(line_split[2])-int(line_split[3])),str(int(line_split[5])-int(line_split[2]))]
        elif line_split[1].startswith(sys.argv[3]):
            species2 = line_split[6]
            species2_count += 1
            if line_split[4] ==	'+':
                species2_info = [line_split[1],line_split[4],line_split[2],str(int(line_split[2])+int(line_split[3]))]
            else:
                species2_info = [line_split[1],line_split[4],str(int(line_split[5])-int(line_split[2])-int(line_split[3])),str(int(line_split[5])-int(line_split[2]))]
    else:
        if species1_count != 1 or species2_count != 1:
            sys.stderr.write('Error: A block did not have both species as 1-to-1 orthologs.\n')
            sys.stderr.write(str(species1_info)+str(species2_info)+'\n')
            sys.exit(1)
        if len(species1) >= int(sys.argv[4]) and len(species2) >= int(sys.argv[4]):
            start = int(species1_info[2])
            for n in range(0,len(species1)-int(sys.argv[4])+1):
                if int(sys.argv[5]) == 100:
                    if species1[n:(n+int(sys.argv[4]))].upper() == species2[n:(n+int(sys.argv[4]))].upper().replace('-',':').replace('N',':'):
                        exact = len(species1[n:(n+int(sys.argv[4]))].upper())
                    else:
                        exact = 0
                else:
                    exact = similar(species1[n:(n+int(sys.argv[4]))].upper(),species2[n:(n+int(sys.argv[4]))].upper().replace('-',':').replace('N',':'))
                if exact >= (int(sys.argv[4]) * int(sys.argv[5]) / 100):
                    sys.stdout.write(species1_info[0]+'\t'+str(start)+'\t'+str(start+int(sys.argv[4]))+'\n')
                if species1[n] != '-':
                    start += 1
        species1,species1_info,species1_count,species2,species2_info,species2_count = clear()
