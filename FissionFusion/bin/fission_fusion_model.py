#!/usr/bin/env python

import argparse, multiprocessing, random, sys, os
from argparse import RawDescriptionHelpFormatter

# Parse arguments
parser = argparse.ArgumentParser(description='Model chromosome fissions and fusions to count the number of fusions that\nperfectly reverse a prior fission.', epilog='This program has several assumptions:\n  1. Only one fission is allowed per chromosome.\n  2. All fissions occur first, followed by all fusions.\n  3. Chromosome selection for each fission is random.\n  4. Chromosome selection and orientation for each fusion is random.', formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--iterations", metavar='INT', help="number of iterations to run [1]", type=int, default=1)
parser.add_argument("-o", "--output", metavar='STR', help="output prefix [out]", type=str, default='out')
parser.add_argument("-p", "--processors", metavar='INT', help="number of processors [1]", type=int, default=1)
required = parser.add_argument_group('required arguments')
required.add_argument("karyotype", help="starting karyotype (n)", type=int)
required.add_argument("fission", help="number of fissions", type=int)
required.add_argument("fusion", help="number of fusions", type=int)
args = parser.parse_args()

# Error if fissions, fusions, or karyotype are < 1
if args.fission < 1 or args.fusion < 1 or args.karyotype < 1:
    sys.stderr.write('ERROR: fissions ('+str(args.fission)+'), fusions ('+str(args.fusion)+'), and karyotype ('+str(args.karyotype)+') must each be greater than zero.\n')
    sys.exit(1)

# Error if fusions are >= (karyotype + fissions)
elif args.fusion >= (args.karyotype + args.fission):
    sys.stderr.write('WARNING: fusions ('+str(args.fusion)+') must be less than sum of fissions and karyotype ('+str(args.fission+args.karyotype)+') ... Wrote negative output to file '+args.output+'.summary.\n')
    out_summary = open(args.output+'.summary','w')
    out_summary.write('0\t-1\n1\t-1\n')
    out_summary.close()
    sys.exit(0)

# Define function
def run_function(num):

    # Set starting karyotype
    starting_kar = range(1,args.karyotype+1)

    # Iterate through fissions
    fission_kar = []
    fission_count = []
    for n in range(0,args.fission):
        fission = random.choice(starting_kar)
        starting_kar.remove(fission)
        fission_count.append(str(fission))
        fission_kar.append(str(fission)+'A_'+str(fission)+'B')
        fission_kar.append(str(fission)+'C_'+str(fission)+'D')
    if len(starting_kar) > 0:
        for n in starting_kar:
            fission_kar.append(str(n)+'A_'+str(n)+'B_'+str(n)+'C_'+str(n)+'D')

    # Iterate through fusions
    fusion_kar = fission_kar
    positive_negative = [0,1]
    for n in range(0,args.fusion):
        fusion = random.sample(fusion_kar,2)
        fusion_kar.remove(fusion[0])
        fusion_kar.remove(fusion[1])
        A_order = random.choice(positive_negative)
        B_order = random.choice(positive_negative)
        if A_order == 1:
            temp_fusion = fusion[0].split('_')
            new = False
            for n in range(len(temp_fusion)-1,-1,-1):
                if not new:
                    new = temp_fusion[n]
                else:
                    new += '_'+ temp_fusion[n]
        else:
            new = fusion[0]
        if B_order == 1:
            temp_fusion = fusion[1].split('_')
            for n in range(len(temp_fusion)-1,-1,-1):
                new += '_'+ temp_fusion[n]
        else:
            new = new+'_'+fusion[1]
        fusion_kar.append(new)

    # Find fusions that perfectly reverse a prior fission
    count = 0
    for n in fission_count:
        if any(n+'A_'+n+'B_'+n+'C_'+n+'D' in s for s in fusion_kar):
            count += 1
    return count

# Run parallelized for set iterations
pool = multiprocessing.Pool(args.processors)
final_count = pool.map(run_function, [num for num in range(0,args.iterations)])
pool.close()

# Write raw counts and summary counts to output files
out_raw = open(args.output+'.raw','w')
out_summary = open(args.output+'.summary','w')
out_raw.write('\n'.join(map(str,final_count))+'\n')
for n in range(0,max(final_count)+1):
    out_summary.write(str(n)+'\t'+str(sum(i >= n for i in final_count))+'\n')
out_raw.close()
out_summary.close()
