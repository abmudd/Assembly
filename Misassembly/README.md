# LoReM (Long Read Misassemblies)

Python script to bin and identify regions of the genome assembly with low spanning coverage and high number of read terminals (starts and ends) using 10X linked reads, PacBio reads, or Nanopore reads. For 10X linked reads, the input bam is the sorted and indexed bam with a BX optional field, such as the output from longranger align. For PacBio or Nanopore reads, the input bam is the sorted and indexed bam from your choice of aligner.

Version 1.0

Citing: unpublished

## Prerequisite Python modules:

```
argparse, collections, gzip, numpy, os, pysam, scipy, sys
```

## Usage

```
usage: LoReM.py [-h] [-b INT] [-l INT] [-o STR] [-s FLOAT] [-t FLOAT] [-v]
                [-d INT] [-r INT] -f STR -i STR -y STR
```

## Required arguments

```
-f STR, --fasta STR        genome fasta
-i STR, --in_bam STR       sorted and indexed bam
-y STR, --type STR         type of input: 10X, PacBio, or Nanopore
```

## 10X optional arguments

```
-d INT, --distance INT     minimum distance between reads with same 10X barcode being called new 10X molecules [50000]
-r INT, --reads INT        minimum reads required per 10X molecule [4]
```

## Optional arguments

```
-h, --help                 show help message and exit
-b INT, --binsize INT      bin size for analysis of sequence [1000]
-l INT, --length INT       minimum 10X molecule or PacBio/Nanopore read length [1000]
-o STR, --output STR       output file prefix [output]
-s FLOAT, --spanning FLOAT spanning coverage threshold: calculates from p-value if <1; otherwise sets cutoff value [0.05]
-t FLOAT, --terminal FLOAT minimum terminal (starts/ends) threshold: calculates from p-value if <1; otherwise sets cutoff value [0.05]
-v, --version              show version info and exit
```

## Output

LoReM writes out four files along with information to stderr:

```
1. output.full.bed.gz      (if -y 10X) bgzip bed file of 10X linked read from bam; notes position, start position, stop position, barcode name, read coverage, and density of starts on read
2. output.full.bed.gz.tbi  (if -y 10X) tabix file of output.full.bed.gz
3. output.tsv.gz           gzip tsv file with chromosome, start position, stop position, long read spanning coverage, number of long read starts, and number of long read ends for all bins in the assembly
4. output.bed              bed file with chromosome, start position, stop position, and name (coverage:##;starts:##;ends:##) for bins with potential misassemblies
5. stderr                  log of number of scaffolds in fasta, reads in bam, identified reads with spanning coverage / starts / ends, and set thresholds
```

## Suggestion for rerunning 10X data

To rerun 10X data without the time-consuming recalculation of linked read locations, convert the output.full.bed.gz to a bam file with bedtools bedToBam and then rerun LoReM using the bam file with flag -y set to either PacBio or Nanopore.
