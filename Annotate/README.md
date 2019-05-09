# postGeMoMa.py

Python pipeline to reformat and rewrite the GeMoMa gff output and then extract the CDS and protein sequences using Kent tools.

Version 1.0

Publication status: unpublished

## Prerequisites:

```
kent (tested on binaries downloaded 2019-03-05)
```

Python modules:
```
argparse
datetime
os
pysam
subprocess
sys
```

## Usage: 

```
usage: postGeMoMa.py [-h] [-a STR] [-g STR] [-k STR] [-o STR] [-p STR]
                     [-r STR] [-s STR] [-v]
                     fasta gff
```

## Required arguments:

```
fasta                 Input genome fasta file
gff                   Input GeMoMa gff file
```

## Optional arguments:

```
-h, --help            show this help message and exit
-a STR, --annotation STR
                      annotation version number
-g STR, --genome STR  genome version number
-k STR, --kent STR    full path to kent directory [~/tools/kent]
-o STR, --output STR  output file prefix [output]
-p STR, --replace STR
                      replace any evidence old names with new names in this
                      tab delimited file
-r STR, --remove STR  remove any GeMoMa genes using any gene names in this
                      file as evidence
-s STR, --species STR
                      species name
-v, --version         show program's version number and exit
```
