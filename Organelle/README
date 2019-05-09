# organelle_pipeline.py

Bash pipeline to download organelle sequences from species close to the sequenced species using NCBI accessions, create the respective NOVOPlasty config files, and run NOVOPlasty.

Version 1.0

Publication status: unpublished

## Prerequisites:

```
efetch (tested on v7.10) 
NOVOPlasty (tested on v2.6.3)
perl (tested on v5.18.2)
seqkit (tested on v0.7.2-dev)
```

Python modules:
```
argparse
os
subprocess
sys
```

## Usage: 

```
usage: organelle_pipeline.py [-h] [-e STR] [-n STR] [-s STR] [-v] [-x] -a STR
                             -f STR -i INT -l INT -r STR
```

## Required arguments:

```
-a STR, --accessions STR
                      list of NCBI accessions
-f STR, --freads STR  full path to forward reads
-i INT, --isize INT   insert size
-l INT, --rlen INT    read length
-r STR, --rreads STR  full path to reverse reads
```

## Optional arguments:

```
-h, --help            show this help message and exit
-e STR, --efetch STR  full path to efetch [efetch]
-n STR, --novoplasty STR
                      full path to NOVOPlasty.pl [NOVOPlasty.pl]
-s STR, --seqkit STR  full path to seqkit [seqkit]
-v, --version         show program's version number and exit
-x, --tenx            flag if input data is 10X [False]
```
