#!/usr/bin/env python

import argparse, sys, os
from fisher import pvalue ## install from https://pypi.python.org/pypi/fisher/

# Exit after fatal error message
def fatal_error_exit(status=1):
    in_vcf.close()
    if args.output != 'stdout':
        out_AD.close()
    sys.exit(status)

# Extracts input vcf file
def extract_input(in_vcf, out_AD):
    contig = {}
    header_a_seen = False
    header_b_seen = False

    # header lines part a
    for line in in_vcf:
        if line.startswith('##fileformat'): # write file format header
            header_a_seen = True
        elif line.startswith('##contig'): # extract contig header and create dictionary
            line = line.rstrip().split('<', 1)[1].split('>')[0].replace('=', ',').split(',')

            # find the index of the contig ID term 
            try:
                header_ID = line.index('ID')
            except ValueError:
                out_log.write("["+os.path.basename(sys.argv[0])+"]: Error: Invalid vcf header - missing contig ID\n")
                fatal_error_exit();
            
            # find the index of the contig length term 
            try:
                header_length = line.index('length')
            except ValueError:
                out_log.write("["+os.path.basename(sys.argv[0])+"]: Error: Invalid vcf header - missing contig length\n")
                fatal_error_exit();

            # create dictionary with the contig ID and length
            contig[line[header_ID+1]] = line[header_length+1]
        elif line.startswith('##'): # extract remaining header
            continue
        else:
            break

    if not args.quiet:
        out_log.write(os.path.basename(sys.argv[0])+'Command='+' '.join(sys.argv)+'\n')

    in_vcf.seek(0)
    # header lines part b
    for line in in_vcf:
        if line.startswith('#CHROM'):
            header_b_seen = True # header has now been seen
            line = line.rstrip().split('\t') # split header for processing
            out_AD.write('CHROM\tPOS\t'+'\t'.join(line[9:])+'\n')
            len_header = len(line) # length of the header

            # find the index of parent1
            if args.parent1:
                try:
                    p1idx = line.index(args.parent1)
                except ValueError:
                    out_log.write("Invalid parent ID: %s\n" % args.parent1)
                    fatal_error_exit();

            # find the index of parent2
            if args.parent2:
                try:
                    p2idx = line.index(args.parent2)
                except ValueError:
                    out_log.write("Invalid parent ID: %s\n" % args.parent2)
                    fatal_error_exit();

            break

    if not header_a_seen or not header_b_seen: # if there was no fileformat header or CHROM header seen in the file
        out_log.write("["+os.path.basename(sys.argv[0])+"]: Error: Invalid vcf file - header not found\n")
        fatal_error_exit();

    in_vcf.seek(0)
    for line in in_vcf: # for the remaining lines in in_vcf, which contain locus genotype information
        if not line.startswith('#'):
            orig_line = line.strip()
            line = orig_line.split('\t')
            len_line = len(line)
            if len_line != len_header: # ensure that there are the right number of progeny
                if not args.quiet:
                    out_log.write("["+os.path.basename(sys.argv[0])+"]: Warning: Site " + line[0] + ":" + line[1] + " was skipped - incorrect number of progeny detected\n")
                else:
                    continue
            elif 'AD' not in line[8]: # ensure that the format cell specified the AD field
                if not args.quiet:
                    out_log.write("["+os.path.basename(sys.argv[0])+"]: Warning: Site " + line[0] + ":" + line[1] + " was skipped - no AD field detected\n")
                else:
                    continue
            elif int(line[1]) > int(contig[line[0]]): # ensure that site position is not greater than contig length
                if not args.quiet:
                    out_log.write("["+os.path.basename(sys.argv[0])+"]: Warning: Site " + line[0] + ":" + line[1] + " was skipped - site position is greater than contig length\n")
                else:
                    continue
            elif ',' in line[4]: # skip sites that aren't biallelic
                continue
            elif line[3] not in ('A','C','G','T') or line[4] not in ('A','C','G','T'): # skip indels
                continue
            elif (1 - (orig_line.count('./.'))/float(len_line-9)) < args.sites: # skip if the total fraction missing is less than specified
                continue
            else:
                ad_idx = 0 # if AD is only field
                if ':' in line[8]: # if multiple fields in format
                    line[8] = line[8].split(':')
                    ad_idx = line[8].index('AD')
                    line[8] = line[8][ad_idx] # sets format field to AD

                count_missing = 0
                count_valid = 0
                ref_p1 = 0
                alt_p1 = 0
                ref_p2 = 0
                alt_p2 = 0
                ref_off = 0
                alt_off = 0

                for n in range(9, len_line): # for each individual
                    if line[n] == './.':
                        line[n] = '.'
                        count_missing += 1
                    else:
                        data = line[n].split(':')[ad_idx].split(',')
                        if data[0] == '0' and data[1] == '0': # convert AD to . if ref and alt AD are both 0
                            line[n] = '.'
                            count_missing += 1
                        elif data[0] == '.' or data[1] == '.': # convert AD to . if ref or alt AD are .
                            line[n] = '.'
                            count_missing += 1
                        else: # extract AD for parent 1/2 and offspring
                            if args.parent1 and n == p1idx:
                                ref_p1 += float(data[0])
                                alt_p1 += float(data[1])
                            elif args.parent2 and n == p2idx:
                                ref_p2 += float(data[0])
                                alt_p2 += float(data[1])
                            else:
                                ref_off += float(data[0])
                                alt_off += float(data[1])
                            line[n] = (',').join(data[0:2])
                            count_valid += 1

                # Skip if sum AD for parent 1/2 is less than 3
                if (ref_p1 + alt_p1) < 3 or (ref_p2 + alt_p2) < 3:
                    continue

                # Skip if more than one parental AD is zero:
                count_zero = 0
                if ref_p1 == 0:
                    count_zero += 1
                if alt_p1 == 0:
                    count_zero += 1
                if ref_p2 == 0:
                    count_zero += 1
                if alt_p2 == 0:
                    count_zero += 1
                if count_zero > 1:
                    continue

                Minor = False
                # If offspring ref AD > offspring alt AD, alt is minor allele
                if ref_off > alt_off:
                    Minor = '1'
                    # Calculate Fisher's exact test for offspring
                    p_off = pvalue(ref_off, alt_off, int((ref_off+alt_off)*0.75), int((ref_off+alt_off)*0.25))

                    # If parent 1 alt AD > parent 2 alt AD, parent 1 is het and parent 2 is hom
                    if alt_p1 > alt_p2:
                        # Calculate Fisher's exact test for parents
                        p_p1 = pvalue(ref_p1, alt_p1, int((ref_p1+alt_p1)*0.50), int((ref_p1+alt_p1)*0.50))
                        p_p2 = pvalue(ref_p2, alt_p2, (ref_p2+alt_p2), 0)

                    # If parent 1 alt AD < parent 2 alt AD, parent 1 is hom and parent 2 is het 
                    elif alt_p1 < alt_p2:
                        # Calculate Fisher's exact test for parents
                        p_p1 = pvalue(ref_p1, alt_p1, (ref_p1+alt_p1), 0)
                        p_p2 = pvalue(ref_p2, alt_p2, int((ref_p2+alt_p2)*0.50), int((ref_p2+alt_p2)*0.50))

                    # Skip if parent 1 alt AD = parent 2 alt AD
                    else:
                        continue

                # If offspring ref AD < offspring alt AD, ref is minor allele
                elif ref_off < alt_off:
                    Minor = '0'
                    # Calculate Fisher's exact test for offspring
                    p_off = pvalue(ref_off, alt_off, int((ref_off+alt_off)*0.25), int((ref_off+alt_off)*0.75))

                    # If parent 1 ref AD > parent 2 ref AD, parent 1 is het and parent 2 is hom
                    if ref_p1 > ref_p2:
                        # Calculate Fisher's exact test for parents
                        p_p1 = pvalue(ref_p1, alt_p1, int((ref_p1+alt_p1)*0.50), int((ref_p1+alt_p1)*0.50))
                        p_p2 = pvalue(ref_p2, alt_p2, 0, (ref_p2+alt_p2))

                    # If parent 1 ref AD < parent 2 ref AD, parent 1 is hom and parent 2 is het
                    elif ref_p1 < ref_p2:
                        # Calculate Fisher's exact test for parents
                        p_p1 = pvalue(ref_p1, alt_p1, 0, (ref_p1+alt_p1))
                        p_p2 = pvalue(ref_p2, alt_p2, int((ref_p2+alt_p2)*0.50), int((ref_p2+alt_p2)*0.50))

                    # Skip if parent 1 ref AD = parent 2 ref AD
                    else:
                        continue

                # Skip if offspring ref AD = offspring alt AD
                else:
                    continue

                # Write out locus to AD file
                if (count_valid/float(count_missing+count_valid)) >= args.sites and float(str(p_off.two_tail)) >= args.p_value and float(str(p_p1.two_tail)) >= args.p_value and float(str(p_p2.two_tail)) >= args.p_value:
                    out_AD.write(line[0]+'\t'+line[1]+'\t'+'\t'.join(map(str, line[9:len_line]))+'\n')

# Parse arguments
parser = argparse.ArgumentParser(description='This script extracts biallelic pseudo-testcross SNPs from a vcf file and outputs the sites in AD format.')
parser.add_argument("-p1", "--parent1", metavar='STR', help="sample ID for parent 1", type=str)
parser.add_argument("-p2", "--parent2", metavar='STR', help="sample ID for parent 2", type=str)
parser.add_argument("-p", "--p_value", metavar='FLOAT', help="p-value cutoff for pseudo-testcross determination [0.10]", type=float, default='0.10')
parser.add_argument("-q", "--quiet", help="run without stderr output", action='store_true')
parser.add_argument("-s", "--sites", metavar='FLOAT',help="minimum fraction of sites with data [0.50]", type=float, default='0.50')
parser.add_argument("-l", "--log", metavar='STR', help="log name [stderr]", type=str, default='stderr')
parser.add_argument("-o", "--output", metavar='STR', help="output name [stdout]", type=str, default='stdout')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("input", help="input vcf", type=str)
args = parser.parse_args()

# Open input file and set output file
in_vcf = open(args.input, 'r')
out_AD = sys.stdout
if args.output != 'stdout':
    out_AD = open(args.output, 'w')
out_log = sys.stderr
if args.log != 'stderr':
    out_log = open(args.log, 'w')

extract_input(in_vcf, out_AD)

# Close input and output files 
in_vcf.close()
if args.output != 'stdout':
    out_AD.close()
if args.log != 'stderr':
    out_log.close()
