# cactus_filter.py

Python pipeline to run cactus with two species and filter the output. Cactus jobs are submitted to a SLURM cluster with sbatch using the defined memory, processor, and time flags, whereas filter jobs are run locally. All cactus variables must either be in the appropriate paths before running or have the -l flag set to load cactus as a module.

Version 1.0

Citing: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS (2020). Analysis of muntjac deer genome and chromatin architecture reveals rapid karyotype evolution. Communications Biology 3, 480. doi: [10.1038/s42003-020-1096-9](https://doi.org/10.1038/s42003-020-1096-9).

## Prerequisites:

```
cactus (tested on commit e4d0859)
chainCleaner (tested on commit aacca59)
HAL (tested on commit f7287c8)
kent (tested on binaries downloaded 2019-03-05)
```

Python modules:
```
argparse
datetime
os
subprocess
sys
```

## Usage: 

```
usage: cactus_filter.py [-h] [-a STR] [-c STR] [-f STR] [-g STR] [-k STR] [-l]
                        [-m INT] [-n STR] [-o STR] [-p INT] [-r INT] [-s STR]
                        [-t INT] [-v]
                        qry_name qry_fasta ref_name ref_fasta tree
```

## Required arguments:

```
qry_name              name of query species in tree
qry_fasta             query fasta file
ref_name              name of reference species in tree
ref_fasta             reference fasta file
tree                  phylogenetic tree [string or file]
```

## Optional arguments:

```
-h, --help            show this help message and exit
-a STR, --hal STR     full path to hal directory [~/tools/hal]
-c STR, --chain STR   chain file to liftover reference assembly [none]
-f STR, --filter STR  netFilter synteny level: chimp for human/chimp, mouse
                      for human/mouse, or none for no synteny netFilter
                      [none]
-g STR, --gat STR     full path to GenomeAlignmentTools directory
                      [~/tools/GenomeAlignmentTools]
-k STR, --kent STR    full path to kent directory [~/tools/kent]
-l, --load            load cactus and samtools modules
-m INT, --memory INT  RAM memory for job submission in GB [110]
-n STR, --nodedir STR
                      SLURM node temp directory [$TMPDIR]
-o STR, --output STR  output prefix [out]
-p INT, --processor INT
                      processors for job submission [64]
-r INT, --colinrun INT
                      cutoff size for runs of colinearity [1000]
-s STR, --samtools STR
                      path to samtools [samtools]
-t INT, --time INT    hours for job submission [168]
-v, --version         show program's version number and exit
```
