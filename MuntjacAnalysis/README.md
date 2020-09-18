# Muntjac analysis scripts

Version 1.0

Citing: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS (2020). Analysis of muntjac deer genome and chromatin architecture reveals rapid karyotype evolution. Communications Biology 3, 480. doi: [10.1038/s42003-020-1096-9](https://doi.org/10.1038/s42003-020-1096-9).

## 1. run_fission_fusion.sh

Bash script to model all fission and fusion possibilities from an initial karyotype to one million iterations in order to count the number of fusions that perfectly reverse a prior fission. This script submits jobs to a SLURM or SGE cluster.

### Prerequisites:

```
python (tested on v2.7-anaconda)
qbatch (tested on v0.4.0-mudd; https://bitbucket.org/rokhsar-lab/gbs-analysis/src)
R (tested on v3.5.3)
```

Python modules:
```
argparse
multiprocessing
os
random
sys
```

R libraries:
```
gplots
RColorBrewer
```

### Usage: 

```
Usage: run_fission_fusion.sh <1n karyotype>
```

### Notes:

```
This program has several assumptions:
  1. Only one fission is allowed per chromosome.
  2. All fissions occur first, followed by all fusions.
  3. Chromosome selection for each fission is random.
  4. Chromosome selection and orientation for each fusion is random.

This script assumes a standard UNIX/Linux install and access to a SGE or SLURM job-scheduling system via qbatch (see https://bitbucket.org/rokhsar-lab/gbs-analysis/src).
```

## 2. HiCbins_1Mb.py

Python script that extracts number of contacts per 1 Mb bin for the M. muntjak assembly from Juicer merged no dups output.

### Prerequisite Python modules:

```
math
os
sys
```

### Help message:

```
Usage: HiCbins_1Mb.py <in.mnd>
This script extracts number of contacts per 1 Mb bin for the M. muntjak assembly from Juicer merged no dups output.
```

## 3. HiCbins_100kb.py

Python script that extracts number of contacts per 100 kb bin for the M. muntjak assembly from Juicer merged no dups output.

### Prerequisite Python modules:

```
math
os
sys
```

### Help message:

```
Usage: HiCbins_1Mb.py <in.mnd>
This script extracts number of contacts per 100 kb bin for the M. muntjak assembly from Juicer merged no dups output.
```

## 4. extractOrthoVenn.py

Python script that extracts amino acid sequences of 1-to-1 orthologs from OrthoVenn output and outputs an interleaved fasta file.

### Prerequisite Python modules:

```
os
pysam
sys
```

### Help message:

```
Usage: extractOrthoVenn.py <in.orthologs> <species1> <species1.pro.aa> <species2> <species2.pro.aa>
This script extracts amino acid sequences of 1-to-1 orthologs from OrthoVenn output and outputs an interleaved fasta file.
```

## 5. extract2speciesmaf.py

Python script that filters maf alignments where each block contains exactly one species 1 alignment and one species 2 alignment and requires the alignment to have at least 1 base of overlap. The output is either a new pairwise MAF file or a BEDPE file.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: extract2speciesmaf.py <in.maf> <species1> <species2> <MAF/BED>
This script filters maf alignments where each block contains exactly one species 1 alignment and one species 2 alignment and requires the alignment to have at least 1 base of overlap. The output is either a new pairwise MAF file or a BEDPE file.
```

## 6. mcscan_convert_links.py

Python script that converts the Circos links into JCVI files. Count is 1, first is A, and second is B if not provided.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: mcscan_convert_links.py <in.circos_links> <out.prefix> <count> <first> <second>
This script converts the Circos links into JCVI files. Count is 1, first is A, and second is B if not provided.
```

## 7. mcscan_invert_chr.py

Python script that inverts the bed locations for a particular chromosome.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: mcscan_invert_chr.py <in.bed> <in.fai> <chr_to_invert>
This script inverts the bed locations for a particular chromosome.
```
