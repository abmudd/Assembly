# Genetic map scripts

Python scripts F1call.py and F2call.py extract biallelic pseudo-testcross SNPs in vcf format from a F1 cross or F2 cross, respectively. The output AD (allele depth) format files are converted into JoinMap loc format with Python script AD2loc.py.

Version 1.0

Publication status: unpublished

## Prerequisite Python modules:

```
argparse
fisher
os
sys
```

## 1. F1call.py

```
usage: F1call.py [-h] [-p1 STR] [-p2 STR] [-p FLOAT] [-q] [-s FLOAT] [-l STR]
                 [-o STR] [-v]
                 input

This script extracts biallelic pseudo-testcross SNPs from a vcf file and
outputs the sites in AD format.

optional arguments:
  -h, --help            show this help message and exit
  -p1 STR, --parent1 STR
                        sample ID for parent 1
  -p2 STR, --parent2 STR
                        sample ID for parent 2
  -p FLOAT, --p_value FLOAT
                        p-value cutoff for pseudo-testcross determination
                        [0.10]
  -q, --quiet           run without stderr output
  -s FLOAT, --sites FLOAT
                        minimum fraction of sites with data [0.50]
  -l STR, --log STR     log name [stderr]
  -o STR, --output STR  output name [stdout]
  -v, --version         show program's version number and exit

required arguments:
  input                 input vcf
```

## 2. F2call.py

```
usage: F2call.py [-h] [-p1 STR] [-p2 STR] [-g1 STR] [-g2 STR] [-p FLOAT] [-q]
                 [-s FLOAT] [-l STR] [-o STR] [-v]
                 input mode

This script extracts biallelic pseudo-testcross SNPs from a vcf file and
outputs the sites in AD format.

optional arguments:
  -h, --help            show this help message and exit
  -p1 STR, --parent1 STR
                        sample ID for parent 1
  -p2 STR, --parent2 STR
                        sample ID for parent 2
  -g1 STR, --grandparent1 STR
                        sample ID for grandparent 1
  -g2 STR, --grandparent2 STR
                        sample ID for grandparent 2
  -p FLOAT, --p_value FLOAT
                        p-value cutoff for pseudo-testcross determination
                        [0.10]
  -q, --quiet           run without stderr output
  -s FLOAT, --sites FLOAT
                        minimum fraction of sites with data [0.50]
  -l STR, --log STR     log name [stderr]
  -o STR, --output STR  output name [stdout]
  -v, --version         show program's version number and exit

required arguments:
  input                 input vcf
  mode                  analysis mode: b (grand/parents), g (grand)
```

## 3. AD2loc.py

```
usage: AD2loc.py [-h] [--output STR] [--population STR] [--remove STR]
                  parent1 parent2 input

Converts pseudo-testcross SNPs in AD format to JoinMap loc format.

optional arguments:
  -h, --help        show this help message and exit
  --output STR      output name [stdout]
  --population STR  population name [Blue]
  --remove STR      remove pattern from locus ids [caffold]
  -v, --version     show program's version number and exit

required arguments:
  parent1           sample ID for parent 1
  parent2           sample ID for parent 2
  input             input AD
```
