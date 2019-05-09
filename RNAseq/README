# STARalign.sh

Bash pipeline to align RNA-seq reads to a genome with splice junctions using STAR.

Version 1.0

Publication status: unpublished

## Prerequisites:

```
qbatch (tested on v0.3.1-mudd; https://bitbucket.org/rokhsar-lab/gbs-analysis/src)
samtools (tested on v1.6)
STAR (tested on v2.5.3a)
```

## Usage: 

```
STARalign.sh --workdir STR --RNA-seq-dir STR --ref-genome STR
             [--cthreads INT] [--index-mem INT] [--lthreads INT] [--module STR] 
             [--runtime HH:MM:SS] [--sj-support INT] [--help] [-h]
```

## Required arguments:

```
-w, --workdir STR       path for working directory
-d, --RNA-seq-dir STR   path for directory with RNA-seq files
-r, --ref-genome STR    path for genome fasta file
```

## Optional arguments:

```
-a, --align-mem INT     memory for STAR align step in GB [110]
-c, --cthreads INT      run pipeline on cluster with INT number of processors [32]
-i, --index-mem INT     memory for STAR index step in GB [240]
-l --lthreads INT       run pipeline locally with INT number of processors [0]
-m, --module STR        module load STR (recommend star/2.5.3a and samtools/1.6)
-t, --runtime HH:MM:SS  runtime for job submissions [12:0:0]
-s, --sj-support INT    minimum support for splice junctions [20]
-h, --help              show this help message and exit
```

## Notes:
1. This script assumes a standard UNIX/Linux install and access to a SGE or SLURM job-scheduling system via qbatch (see https://bitbucket.org/rokhsar-lab/gbs-analysis/src).
2. Running locally/lthreads takes priority over running on the cluster/cthreads. If running locally, caution should be taken about the lthreads flag as the STAR align step is memory intensive. Too many lthreads will cause a segmentation fault.
3. This script assumes that fastq files end in the format [0,1,2][fastq format], are the only fastq files in the RNA-seq-dir, and are all only one fastq format. Accepted fastq formats are .fastq, .fq, .fastq.gz, and .fq.gz. These files should already be adapter trimmed.
4. It is highly recommended that absolute paths be passed via the required flags, as those paths may be used in submission to the cluster.
