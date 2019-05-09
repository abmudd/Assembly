#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.maf> <species1> <species2>\n")
    sys.stderr.write("This script calculates the number of bases aligned and the percent ID for each maf block. The script requires 1-to-1 orthologous alignment of the two species.\n")
    sys.exit(1)

in_maf = sys.stdin
if sys.argv[1] != "-":
    in_maf = open(sys.argv[1],'r')

sys.stdout.write('#S1chr\tS1strand\tS1start\tS1end\tS2chr\tS2strand\tS2start\tS2end\tAlnsize\tAln100%ID\n')

def base_replace(string):
    return string.replace('A','X').replace('C','X').replace('G','X').replace('T','X')
    
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
    elif header and not line.split():
        header = False
    elif line.startswith('a'):
        header = False
    elif line.startswith('s'):
        line_split = line.rstrip().split()
        if line_split[1].startswith(sys.argv[2]):
            species1 = line_split[6]
            species1_count += 1
            if line_split[4] == '+':
                species1_info = line_split[1]+'\t'+line_split[4]+'\t'+line_split[2]+'\t'+str(int(line_split[2])+int(line_split[3]))+'\t'
            else:
                species1_info = line_split[1]+'\t'+line_split[4]+'\t'+str(int(line_split[5])-int(line_split[2])-int(line_split[3]))+'\t'+str(int(line_split[5])-int(line_split[2]))+'\t'
        elif line_split[1].startswith(sys.argv[3]):
            species2 = line_split[6]
            species2_count += 1
            if line_split[4] ==	'+':
                species2_info = line_split[1]+'\t'+line_split[4]+'\t'+line_split[2]+'\t'+str(int(line_split[2])+int(line_split[3]))+'\t'
            else:
                species2_info = line_split[1]+'\t'+line_split[4]+'\t'+str(int(line_split[5])-int(line_split[2])-int(line_split[3]))+'\t'+str(int(line_split[5])-int(line_split[2]))+'\t'
    else:
        if species1_count != 1 or species2_count != 1:
            sys.stderr.write('Error: A block did not have both species as 1-to-1 orthologs.\n')
            sys.stderr.write(str(species1_info)+str(species2_info)+'\n')
            sys.exit(1)
        exact = similar(species1.upper(),species2.upper().replace('-',':').replace('N',':'))
        total = similar(base_replace(species1.upper()),base_replace(species2.upper()).replace('-',':').replace('N',':'))
        sys.stdout.write(species1_info+species2_info+str(total)+'\t'+str(exact)+'\n')
        species1,species1_info,species1_count,species2,species2_info,species2_count = clear()
