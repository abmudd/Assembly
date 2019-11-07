# Homology scripts

Version 1.0

Citing: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS. Muntjac chromosome evolution and architecture. bioRxiv: 772343. doi: [10.1101/772343](https://doi.org/10.1101/772343)

## 1. 4Dextract.py

Python script to extract bases with four fold degeneracy from input fasta and genePred annotation of reference species, informed by conserved amino acids in other species using input MAF alignment with the reference species on top.

### Prerequisite Python modules:

```
argparse
Bio.Seq
os
pysam
sys
```

### Help message:

```
usage: 4Dextract.py [-h] fasta genepred maf species

Python script to extract bases with four fold degeneracy from input fasta and genePred annotation of reference species, informed by conserved amino acids in other species using input MAF alignment with the reference species on top.

optional arguments:
  -h, --help  show this help message and exit
  -v, --version  show program's version number and exit

required arguments:
  fasta       input fasta file
  genepred    input genepred file
  maf         input maf file
  species     comma-separated list of maf species with reference species first
```

## 2. extractOrthoVenn.py

Python script to extract amino acid sequences from OrthoVenn 1-to-1 orthologs.

### Prerequisite Python modules:

```
os
pysam
sys
```

### Help message:

```
Usage: extractOrthoVenn.py <in.orthologs> <species1> <species1.pro.aa> <species2> <species2.pro.aa>
This script extracts amino acid sequences from OrthoVenn 1-to-1 orthologs.
```
