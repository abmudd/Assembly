#!/usr/bin/env python

import argparse, sys, os

# Parse arguments
parser = argparse.ArgumentParser(description='Converts pseudo-testcross SNPs in AD format to JoinMap loc format.')
parser.add_argument("--output", metavar='STR', help="output name [stdout]", type=str, default='stdout')
parser.add_argument("--population", metavar='STR', help="population name [Blue]", type=str, default='Blue')
parser.add_argument("--remove", metavar='STR', help="remove pattern from locus ids [caffold]", type=str, default='caffold')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("parent1", help="sample ID for parent 1", type=str)
required.add_argument("parent2", help="sample ID for parent 2", type=str)
required.add_argument("input", help="input AD", type=str)
args = parser.parse_args()

# Open input file and set output file
in_AD = open(args.input, 'r')
out_loc = sys.stdout
if args.output != 'stdout':
    out_loc = open(args.output, 'w')

n_loc = 0
indv_name_list = []
converted_loc_list = []

# Convert vcf GT into loc format
for line in in_AD:

    if line.startswith('CHROM'): # select header line with individual names
        line = line.rstrip().split('\t')
        len_line = len(line)
        n_ind = len_line - 2 # find the number of individuals
        for n in range(2, len_line):
            indv_name_list.append(line[n]) # append all individual names to list
        p1idx = line.index(args.parent1) # find the index of parent1
        p2idx = line.index(args.parent2) # find the index of parent2

    else:
        line = line.rstrip().split('\t')

        # set flags based on parents' ad
        lmxll = False
        nnxnp = False
        p1_ad = line[p1idx].split(',')
        p2_ad = line[p2idx].split(',')
        min_pad = min(p1_ad[0],p1_ad[1],p2_ad[0],p2_ad[1])
        if min_pad > 3:
            error_on_line = True
        if p1_ad[1] == min_pad or p2_ad[1] == min_pad: # hom parent is 0/0, so nnxnp conversion
            nnxnp = True
        else: # hom parent is 1/1, so lmxll conversion
            lmxll = True

        # if parents are 0/1 and 1/1, convert ad to lm and ll
        if lmxll:
            header = '\t<lmxll>\t\t(ll,lm)'
            for n in range(2, len_line):
                if line[n] == '.': # missing call
                    line[n] = '--'
                else:
                    ad_info = line[n].split(',')
                    if ad_info[0] == '0' or (ad_info[0] == '1' and int(ad_info[1]) >= 4): # hom call
                        line[n] = 'll'
                    else: # het call
                        line[n] = 'lm'

        # if parents are 0/1 and 0/0, convert ad to np and nn
        else:
            header = '\t<nnxnp>\t\t(nn,np)'
            for n in range(2, len_line):
                if line[n] == '.': # missing call
                    line[n] = '--'
                else:
                    ad_info = line[n].split(',')
                    if ad_info[1] == '0' or (ad_info[1] == '1' and int(ad_info[0]) >= 4): # hom call
                        line[n] = 'nn'
                    else: # het call 
                        line[n] = 'np'

        if args.remove and args.remove in line[0]: # replace pattern in locus ids if set
            line[0] = line[0].replace('caffold','')
        len_seq = len(line[1]) # Get length of characers in base location
        if len_seq != 9: # Add leading zeros to set number of characters in base location to 9
            for n in range((len_seq), 9):
                line[1] = '0'+line[1]
        converted_loc_list.append(line[0]+':'+line[1]+header+'\t'+'\t'.join(map(str, line[2:len_line])))
        n_loc += 1 # sum the number of loci

# Write output
out_loc.write('name = '+args.population+'\npopt = CP\nnloc = '+str(n_loc)+'\nnind = '+str(n_ind)+'\n\n') # file header
for item in converted_loc_list:
    out_loc.write(item+'\n') # loc data
out_loc.write('\nindividual names:\n\n') # indv header
for name in indv_name_list:
    out_loc.write(name+'\n') # indv names

# Close input and output files 
in_AD.close()
if args.output != 'stdout':
    out_loc.close()
