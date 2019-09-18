# align_pipeline.sh

Bash pipeline to check single haplotype, identify misjoins, and rescaffold. Intended for use with adapter trimmed Illumina or 10X data. Inclusion of PB data is in development and not recommended for current use.

Version 1.0

Citing: unpublished

## Prerequisites:

```
bedtools (tested on v2.27.0)
blast+ (tested on v2.6.0)
bwa (tested on v0.7.17-r1188)
jellyfish (tested on v2.2.7)
python (tested on v2.7-anaconda)
qbatch (tested on v0.3.1-mudd; https://bitbucket.org/rokhsar-lab/gbs-analysis/src)
R (tested on v3.3.2)
samtools (tested on v1.6)
SSPACE (tested on v3.0)
```

Python modules:
```
argparse
collections
itertools
numpy
os
pysam
sys
```

## Usage: 

```
Usage: align_pipeline.sh --config STR --fasta STR --mode STR [--halfdepth INT]
       [--kminusone INT] [--meraculous] [--module STR] [--workdir STR] [--help]
```

## Required arguments:

```
-c, --config STR        path for config file
-f, --fasta STR         path for assembly
-o, --mode STR          designate mode: prep, repropair, haplotype, scaffold,
                        or finish
```

## Optional arguments:

```
-a, --halfdepth INT     max cutoff for median coverage of half depth [0]
-k, --kminusone INT     kmer size used in assembly minus one [60]
-e, --meraculous        config file is from meraculous [False]
-m, --module STR        module load STR (such as bedtools/2.17.0)
-w, --workdir STR       path for working directory [pwd]
-h, --help              show this help message and exit
```

## Notes:
1. This script requires bedtools, blast, bwa, qbatch, samtools, SSPACE, and R. These can be loaded using the module flag, manually hard set at the top of the script, or included in the PATH.
2. If not a meraculous config file, the config file must be a txt file with each line containing a library name; read 1 full path; read 2 full path; [FR,RF,10X,PB]; insert size; and insert standard deviation; all split by tabs. 10X reads are assumed to be in FR orientation.
3. If you are unsure of the insert size and standard deviation, please run the script estimate_insert_size.pl included with SSPACE.

## Instructions:
1. Start by running prep mode.
2. If you want to change the insert size and standard deviation and rerun proper pairing, alter the input config file and then run repropair mode. Proper pairing is used in misassembly mode.
3. Run haplotype mode if sufficient half-depth scaffolds are present in the *.hist.pdf output.
4. Run scaffold mode to scaffold using SSPACE. Will only run if previous haplotype mode was run.
5. Run finish mode to clean up the temporary files.
