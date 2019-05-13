#!/usr/bin/env python

###############################################################################
# Setup
###############################################################################

# Import modules
# ============================================================

import argparse, datetime, os, pysam, subprocess, sys


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script reformats and rewrites the GeMoMa gff output and then extracts the CDS and protein sequences using Kent tools.')
parser.add_argument('-a', '--annotation', metavar='STR', help='annotation version number', type=str)
parser.add_argument('-g', '--genome', metavar='STR', help='genome version number', type=str)
parser.add_argument("-k", "--kent", metavar='STR', help="full path to kent directory [~/tools/kent]", type=str, default='~/tools/kent')
parser.add_argument('-o', '--output', metavar='STR', help='output file prefix [output]', type=str, default='output')
parser.add_argument('-p', '--replace', metavar='STR', help='replace any evidence old names with new names in this tab delimited file', type=str, default=False)
parser.add_argument('-r', '--remove', metavar='STR', help='remove any GeMoMa genes using any gene names in this file as evidence', type=str, default=False)
parser.add_argument('-s', '--species', metavar='STR', help='species name', type=str)
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument('fasta', metavar='fasta', help='Input genome fasta file', type=str)
required.add_argument('gff', metavar='gff', help='Input GeMoMa gff file', type=str)
args = parser.parse_args()


# Function to check executables
# ============================================================

def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return os.path.abspath(program)
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def check_exe_return(function):
    if not which(function):
        sys.stderr.write('ERROR: '+function+' not found in specified location or not executable\n')
        sys.exit(127)
    else:
        return which(function)


# Check executables
# ============================================================

faToTwoBit = check_exe_return(os.path.expanduser(args.kent)+'/src/faToTwoBit')
gff3ToGenePred = check_exe_return(os.path.expanduser(args.kent)+'/src/gff3ToGenePred')
genePredToProt = check_exe_return(os.path.expanduser(args.kent)+'/src/genePredToProt')


# Set variables
# ============================================================

cds_count = 1
gene_count = 0
fai = {}
with pysam.FastxFile(args.fasta) as fg:
    for entry in fg:
        fai[entry.name] = entry.sequence
new_gene = False
remove_evidence = []
if args.remove:
    for line in open(args.remove):
        remove_evidence.append(line.strip())
replace_evidence = {}
if args.replace:
    for	line in	open(args.replace):
        line = line.strip().split()
        replace_evidence[line[0]] = line[1]
remove_found = False


# Open output
# ============================================================

out_gff = open(args.output+'.gff','w')
out_sh = open(args.output+'.sh','w')


###############################################################################
# Run
###############################################################################

# Set header
# ============================================================

now = datetime.datetime.now()
out_gff.write('##gff-version 3\n')
out_gff.write('#species '+args.species+'\n')
out_gff.write('#genome version '+args.genome+'\n')
out_gff.write('#annotation version '+args.annotation+'\n')
out_gff.write('#date last updated '+now.strftime('%Y-%m-%d')+'\n')


# Parse gff
# ============================================================

for line in open(args.gff):
    if not line.startswith('#'):
        line = line.rstrip().split('\t')
        if line[2] == 'prediction': # Reformat prediction line to gene line
            new_gene = True
            remove_found = False
            for item in remove_evidence: # Exclude any genes with evidence in remove arg
                if item in line[8]:
                    remove_found = True
            if not remove_found:
                for item in replace_evidence:
                    if item in line[8]:
                        line[8] = line[8].replace(item,replace_evidence[item])
                gene_count += 1
                line[2] = 'gene'
                evidence = line[8].split(';')[1].split('=')[1]
                if 'alternative=' in line[8]:
                    evidence = evidence+','+line[8].split(';')[-1].split('"')[1].replace('gene:','').replace('TRANSCRIPT:','')
                line[8] = 'ID=gene'+str(gene_count)+';evidence='+evidence
                out_gff.write('\t'.join(line)+'\n')
        elif line[2] == 'CDS' and not remove_found:
            if new_gene:
                cds_count = 1
                new_gene = False
            else:
                cds_count += 1
            line[8] = 'ID=gene'+str(gene_count)+'_cds'+str(cds_count)+';Parent=gene'+str(gene_count)
            out_gff.write('\t'.join(line)+'\n')


# Write bash script
# ============================================================

out_sh.write('#!/bin/bash\n')
out_sh.write(faToTwoBit+' '+args.fasta+' '+args.output+'.2bit\n')
out_sh.write('sed \'s/GeMoMa\tgene/GeMoMa\tmRNA/\' '+args.output+'.gff | '+gff3ToGenePred
             +' stdin stdout | sort -k2,2 -k4,4n >'+args.output+'.genePred\n')
out_sh.write(genePredToProt+' -cdsFa='+args.output+'.cds.fasta '+args.output+'.genePred '+args.output
             +'.2bit '+args.output+'.pro.aa\n')


# Close output files and run bash script
# ============================================================

out_gff.close()
out_sh.close()
subprocess.call('sh '+args.output+'.sh',shell=True)
exit()
