#!/usr/bin/env python

import argparse, os, pysam, sys
from Bio.Seq import Seq

# Parse arguments
parser = argparse.ArgumentParser(description='Extract bases with four fold degeneracy from input fasta and genePred annotation of reference species, informed by conserved amino acids in other species using input MAF alignment with the reference species on top.')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("fasta", help="input fasta file", type=str)
required.add_argument("genepred", help="input genepred file", type=str)
required.add_argument("maf", help="input maf file", type=str)
required.add_argument("species", help="comma-separated list of maf species with reference species first", type=str)
args = parser.parse_args()

# Get list of species
species_list = [str(item) for item in args.species.split(',')]

# Define fourfold degenerate site codons
fourfold_deg = ["CTT", "CTA", "CTG", "CTC", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "GGT", "GGC", "GGA", "GGG"]

# Import reference fasta file
fai = {}
with pysam.FastxFile(args.fasta) as fg:
    for entry in fg:
        fai[entry.name] = entry.sequence

# Import genepred file and extract 4D codon locations
fourDsite = {}
for entry in open(args.genepred,'r'):
    entry = entry.rstrip().split('\t')
    strand = entry[2]
    start = entry[8].split(',')
    end = entry[9].split(',')
    # Extract CDS bases and locations
    for i in range(0,int(entry[7])):
        CDS = fai[entry[1]][int(start[i]):int(end[i])] # Get CDS bases
        sites = []
        for j in range(int(start[i]),int(end[i])): # Get CDS locations
            sites.append(j)
        # Parse through each codon in the CDS beginning with exon frame offset
        for k in range(int(entry[14].split(',')[0]),len(CDS),3):
            if k+3 >= len(CDS): # Skip ending sequence that doesn't form a 3 base codon
                continue
            elif sites[k+3] - sites[k] != 3: # Skip codons spanning introns
                continue
            elif strand == '+' and CDS[k:k+3].upper() in fourfold_deg:
                fourDsite[entry[1]+':'+str(sites[k])] = '+'
            elif strand == '-' and str(Seq(CDS[k:k+3]).reverse_complement()).upper() in fourfold_deg:
                fourDsite[entry[1]+':'+str(sites[k])] = '-'

# Set blank dictionary for final set of 4D bases
final = {}
for species in species_list:
    final[species] = ''

# Set empty variables
def blank():
    src = False
    start = False
    all_seq = {}
    return src,start,all_seq
    
# Import maf file and identify conserved 4D codons
header = False
src,start,all_seq = blank()
for line in open(args.maf,'r'):
    if line.startswith('#'): # Skip header lines
        header = True
    elif header and not line.split():
        header = False
    elif line.startswith('a'):
        header = False
    elif line.startswith('s'): # Parse through each block and record species names/sequences as well as top (ref) species src/start site
        line = line.rstrip().split()
        if line[1].split('.')[0] == species_list[0]: # Top (ref) species
            src = line[1]
            start = int(line[2])
            all_seq[species_list[0]] = line[6]
        else: # All other species  
            all_seq[line[1].split('.')[0]] = line[6]
#    elif len(all_seq) == 1: # Skip blocks with only one (presumably ref) species
#        sys.stderr.write('ERROR: Skipping block with sequence '+'\t'.join(all_seq)+' due to only having one species.\n')
#        src,start,all_seq = blank()
    elif species_list[0] not in all_seq: # Skip block without ref species
        sys.stderr.write('ERROR: Skipping block with sequence '+'\t'.join(all_seq)+' due to missing reference.\n')
        src,start,all_seq = blank()
    else:
        for i in range(0,len(all_seq[species_list[0]])): # For base in ref sequence
            if all_seq[species_list[0]][i] == '-': # Exclude sites where ref is a gap
                start -= 1
            elif src+':'+str(start+i) in fourDsite: # Select 4D codons
                sites = []
                for j in range(i,len(all_seq[src.split('.')[0]])): # Remove gaps in next ref sequence
                    if len(sites) == 3: # If sites for three bases are identified
                        break
                    elif all_seq[species_list[0]][j] != '-': # Skip sites with gaps in ref
                        sites.append(j)
                conserved = True
                codon = {}
                for species in species_list:
                    if species in all_seq: # Only look at species with sequences
                        codon[species] = ''
                        for k in sites: # Get codon after removing ref seq gaps
                            codon[species] += all_seq[species][k:k+1].upper()
                        if fourDsite[src+':'+str(start+i)] == '-': # Reverse complement if negative strand
                            codon[species] = str(Seq(codon[species]).reverse_complement())
                        if codon[species] not in fourfold_deg or Seq(codon[species]).translate() != Seq(codon[species_list[0]]).translate(): # Check for codon being 4D with conserved amino acid
                            conserved = False
                if conserved: # Only keep 4D codons with conserved amino acid
                    for species in species_list:
                        if species in codon: # Save 4D base for species with data
                            final[species] += codon[species][-1]
                        else: # Save blank for species without data
                            final[species] += '*'
        src,start,all_seq = blank()

# Write out final aligned fasta output
for species in species_list:
    sys.stdout.write('>'+species+'\n'+''.join(final[species])+'\n')
