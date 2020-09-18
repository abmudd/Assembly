# Decontamination pipelines

Bash pipelines to flag contaminant scaffolds using (1) the archaea, bacteria, and virus RefSeq databases as well as the UniVec vector database in general_decon.sh, (2) the mitochrondria assembly in mt_decon.sh, and (3) the nt database in nt_decon.sh. The output of all three pipelines are *.filter files, which contain flagged contaminant scaffolds.

Version 1.0

Citing:
* For general_decon.sh and mt_decon.sh: Mudd AB, Bredeson JV, Baum R, Hockemeyer D, and Rokhsar DS. Muntjac chromosome evolution and architecture. bioRxiv: 772343. doi: [10.1101/772343](https://doi.org/10.1101/772343)
* For nt_decon.sh: unpublished

## Prerequisites:

```
bedtools (tested on v2.27.0)
blast+ (tested on v2.6.0)
python (tested on v2.7-anaconda)
qbatch (tested on v0.3.1-mudd; https://bitbucket.org/rokhsar-lab/gbs-analysis/src)
samtools (tested on v1.6)
```

Python modules:
```
distutils
os
pysam
subprocess
sys
tempfile
```

## 1. general_decon.sh

```
general_decon.sh v1.0 

Usage: general_decon.sh --fasta STR [--archaeadb STR] [--bacteriadb STR] [--module STR]
       [--vectordb STR] [--virusdb STR] [--workdir STR] [--help]

Pipeline to flag contaminant scaffolds using archaea, bacteria, vector, and virus databases.

Required arguments:
       -f, --fasta STR         path to assembly

Optional arguments:
       -a, --archaeadb STR     path to archaea blastdb [scriptsdir/Archaea_RefSeq]
       -b, --bacteriadb STR    path to bacteria blastdb [scriptsdir/Bacteria_RefSeq]
       -m, --module STR        module load STR (such as bedtools/2.27.0)
       -v, --vectordb STR      path to virus blastdb [scriptsdir/UniVec]
       -i, --virusdb STR       path to virus blastdb [scriptsdir/Viral_RefSeq]
       -w, --workdir STR       path to working directory [pwd]
       -h, --help              show this help message and exit

Notes:
   1. This script assumes a standard UNIX/Linux install and access to a SGE or SLURM job-
      scheduling system via qbatch (see https://bitbucket.org/rokhsar-lab/gbs-analysis/src).
   2. This script requires blastn, bedtools, and qbatch. These can be loaded using the
      module flag, manually hard set at the top of the script, or included in the PATH.
   3. The output has the following columns: (1) Scaffold name, (2) Blast hit name, (3) Size
      of blast hit, (4) Percent of scaffold bases in blast hit, (5) Range of scaffold bases
      in blast hit, and (6) Size of scaffold. The output then lists scaffolds to remove and
      scaffolds to check.
```

## 2. mt_decon.sh

```
Usage: mt_decon.sh <assembly.fa> <mitodb>
Pipeline to flag contaminant scaffolds using the mitochrondrial assembly blast database.
Requires blastn, which can be hard set or included in the PATH.
```

## 3. nt_decon.sh

```
nt_decon.sh v1.0 

Usage: nt_decon.sh --approveddb STR --fasta STR --approvedseqid STR [--module STR]
       [--ntdb STR] [--workdir STR] [--help]

Pipeline to flag contaminant scaffolds using the nt database.

Required arguments:
       -a, --approveddb STR    path for approved species db
       -f, --fasta STR         path for assembly
       -s, --approvedseqid STR path for approved sequence ids

Optional arguments:
       -m, --module STR        module load STR (such as bedtools/2.27.0)
       -n, --ntdb STR          path for nt database [scriptsdir/nt]
       -w, --workdir STR       path for working directory [pwd]
       -h, --help              show this help message and exit

Notes:
   1. This script assumes a standard UNIX/Linux install and access to a SGE or SLURM job-
      scheduling system via qbatch (see https://bitbucket.org/rokhsar-lab/gbs-analysis/src).
   2. This script requires blastn, bedtools, qbatch, and samtools. These can be loaded
      using the module flag, manually hard set at the top of the script, or included in
      the PATH.
   3. The approved sequence ids should be the output from blastdbcmd -outfmt "0" of
      those sequences in the ntdb that are verified as close to the sequenced species.
   4. The approved species db must be a txt file with each line containing a db name and
      full path to the db location, split by a tab. These should be assemblies or fasta,
      such as EST or GSS sequences, from closely related species that will be used as a
      second pass for flagging contaminant sequences.
```
