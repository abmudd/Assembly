# Annotation scripts

Version 1.0

Citing:
* For postGeMoMa.py: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS. Muntjac chromosome evolution and architecture. bioRxiv: 772343. doi: [10.1101/772343](https://doi.org/10.1101/772343)
* For largestgenePred.py and filter_trinity.py: unpublished


## 1. postGeMoMa.py

Python pipeline to reformat and rewrite the GeMoMa gff output and then extract the CDS and protein sequences using Kent tools.

Version 1.0

### Prerequisites:

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

### Usage: 

```
usage: postGeMoMa.py [-h] [-a STR] [-g STR] [-k STR] [-o STR] [-p STR]
                     [-r STR] [-s STR] [-v]
                     fasta gff
```

### Required arguments:

```
fasta                 Input genome fasta file
gff                   Input GeMoMa gff file
```

### Optional arguments:

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

## 2. filter_trinity.py

Python script that splits up and filters single- and multi-exon Trinity transcripts from the input bam. Single-exon transcripts should be later filtered with TransDecoder. Multi-exon transcripts are tossed if the first and/or last exon is less than 60 bp, if an intron is less than 60 bp or greater than 300,000 bp, or if the aligned transcript size is less than 250 bp.

### Prerequisite Python modules:

```
argparse, os, pysam, sys
```

### Usage:

```
usage: filter_trinity.py [-h] [-o STR] [-v] bam
```

### Required arguments:

```
bam                   input bam of Trinity transcripts aligned to genome
```

### Optional arguments:

```
  -h, --help            show this help message and exit
  -o STR, --output STR  output file prefix [output]
  -v, --version         show program's version number and exit
```

## 3. largestgenePred.py

Python script that extracts the largest transcript for each gene from the input genePred.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: largestgenePred.py <in.genePred>
This script extracts the largest transcript for each gene from the input genePred.
```
