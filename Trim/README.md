# Read trimming scripts

## 1. nxtrim_pipeline.sh

Bash pipeline to trim Nextera mate pair sequencing and separate mate pair reads and short insert pair end reads in the output files prefix.mp.final_1.fastq.gz, prefix.mp.final_2.fastq.gz, prefix.pe.final_1.fastq.gz, and prefix.pe.final_2.fastq.gz.

Version 1.0

Citing: unpublished

### Prerequisites:

```
nxtrim (tested on v0.4.2)
perl (tested on v5.18.2)
python (tested on v2.7-anaconda)
```

Python modules:
```
argparse
gzip
os
pysam
sys
```

Perl modules:
```
Getopt::Std
IO::Compress::Gzip
IO::Uncompress::Gunzip
```

### Help message:

```
Usage: nxtrim_pipeline.sh <prefix path for MP fastq.gz files>
This script assumes that the files are named prefixpath_1.fastq.gz and prefixpath_2.fastq.gz
```

## 2. trim_10X.py

Python script to adapter trim 10X reads and output fastq files in the expected format of programs scaff10x and break10x.

Version 1.0

Citing: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS. Muntjac chromosome evolution and architecture. bioRxiv: 772343. doi: [10.1101/772343](https://doi.org/10.1101/772343)

### Prerequisite Python modules:

```
argparse
itertools
os
pysam
sys
```

### Help message:

```
usage: trim_10X.py [-h] [-r1 STR] [-r2 STR] [-n STR] in_read1 in_read2

This script adapter trims and outputs fastq files in the expected format of
programs scaff10x and break10x that use trimmed 10X data.

optional arguments:
  -h, --help            show this help message and exit
  -r1 STR, --out_read1 STR
                        output fastq read 1 file name [read-BC_1.fastq]
  -r2 STR, --out_read2 STR
                        output fastq read 2 file name [read-BC_2.fastq]
  -n STR, --out_name STR
                        output name file name [read-BC_1.name]
  -v, --version         show program's version number and exit

required arguments:
  in_read1              input fastq read 1 file name
  in_read2              input fastq read 2 file name
```
