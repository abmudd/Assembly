# Muntjac analysis scripts

Version 1.0

Citing: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS. Muntjac chromosome evolution and architecture. bioRxiv: 772343. doi: [10.1101/772343](https://doi.org/10.1101/772343)

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

Python script that extracts number of contacts per 1 Mb bin for the M. muntjak assembly from mnd.

### Prerequisite Python modules:

```
math
os
sys
```

### Help message:

```
Usage: HiCbins_1Mb.py <in.mnd>
Python script that extracts number of contacts per 1 Mb bin for the M. muntjak assembly from Juicer merged no dups output.
```

## 3. HiCbins_100kb.py

Python script that extracts number of contacts per 100 kb bin for the M. muntjak assembly from mnd.

### Prerequisite Python modules:

```
math
os
sys
```

### Help message:

```
Usage: HiCbins_1Mb.py <in.mnd>
Python script that extracts number of contacts per 100 kb bin for the M. muntjak assembly from Juicer merged no dups output.
```
