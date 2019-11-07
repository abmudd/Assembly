# 4Dextract.py

Python script to extract bases with four fold degeneracy from input fasta and genePred annotation of reference species, informed by conserved amino acids in other species using input MAF alignment with the reference species on top.

Version 1.0

Citing: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS. Muntjac chromosome evolution and architecture. bioRxiv: 772343. doi: [10.1101/772343](https://doi.org/10.1101/772343)

## Prerequisite Python modules:

```
argparse
Bio.Seq
os
pysam
sys
```

## Usage:

```
usage: 4Dextract.py [-h] fasta genepred maf species
```

## Required arguments:

```
  fasta       input fasta file
  genepred    input genepred file
  maf         input maf file
  species     comma-separated list of maf species with reference species first
```

## Optional arguments:

```
  -h, --help  show this help message and exit
  -v, --version  show program's version number and exit
```
