# alignment_plots.py

Python pipeline to writes out alignments of three species chromosomes with centromere locations and repeat histograms.

Version 1.0

Citing: unpublished

## Prerequisite Python modules:

```
argparse
matplotlib
numpy
os
pandas
sys
```

## Usage: 

```
usage: alignment_plots.py [-h] [-o STR] [-x] [-v] -c STR -l STR -f STR -r STR
                          -s STR
```

## Required arguments:

```
-c STR, --centromere STR  concat centromere positions
-f STR, --fai STR         concat fai
-l STR, --links STR       links bedpe
-r STR, --repeatbins STR  concat repeat bins output from gff2bins.py
-s STR, --seqconfig STR   config with all included sequence
```

## Optional arguments:

```
-h, --help                show this help message and exit
-o STR, --output STR      output directory [pwd]
-v, --version             show program's version number and exit
-x, --xtropicalis         flag if mid is Xtropicalis
```

## Example standard input files:

### Centromere

```
Xtropicalis.Chr1	88644696
XlaevisL.Chr1	86160887
XlaevisS.Chr1	75445275
```

### Fai

```
Xtropicalis.Chr1	217471166	114	100	101
XlaevisL.Chr1	219879705	7	60	61
XlaevisS.Chr1	180018555	223544381	60	61
```

### Links

```
Xtropicalis.Chr1	9916668	10108055	XlaevisL.Chr1	666827	863361
Xtropicalis.Chr1	10108483	10137569	XlaevisL.Chr1	864695	895504
Xtropicalis.Chr1	10137589	10191968	XlaevisL.Chr1	898531	943596
Xtropicalis.Chr1	10194355	10197744	XlaevisL.Chr1	944869	953000
XlaevisL.Chr1	666827	863361	Xtropicalis.Chr1	9916668	10108055
XlaevisL.Chr1	864695	895504	Xtropicalis.Chr1	10108483	10137569
XlaevisL.Chr1	898531	943596	Xtropicalis.Chr1	10137589	10191968
XlaevisL.Chr1	944869	953000	Xtropicalis.Chr1	10194355	10197744
Xtropicalis.Chr1	130008	132245	XlaevisS.Chr1	22820	24720
Xtropicalis.Chr1	246768	341080	XlaevisS.Chr1	81736	116290
Xtropicalis.Chr1	341097	364955	XlaevisS.Chr1	116317	149755
Xtropicalis.Chr1	400575	403009	XlaevisS.Chr1	155610	157053
XlaevisS.Chr1	155610	157053	Xtropicalis.Chr1	400575	403009
XlaevisS.Chr1	164183	225233	Xtropicalis.Chr1	423329	580186
XlaevisS.Chr1	431373	488056	Xtropicalis.Chr1	618317	738324
XlaevisS.Chr1	495670	507125	Xtropicalis.Chr1	924794	944654
```

## Sequence config:

The config file has four columns: (1) The chromosomes to visualize in src notation (i.e. species.chromosome), (2) True or False for reversing the chromosome, (3) whether to add additional Mb before the particular chromosome on that line, and (4) whether to place the chromosome on the top (high), middle (mid), or low (bottom) lines. Only one chromosome can be visualized on the middle line. If multiple chromosomes are visualized on the top or bottom lines, the order from left to right will match the order in the config file from top to bottom. Each line is scaled independently.

### Example file

```
#Src	Rev	ExtraMb	HighMidLow
XlaevisL.Chr1	False	0	High
Xtropicalis.Chr1	False	0	Mid
XlaevisS.Chr1	False	0	Low
```

## Repeat bins input:

For the repeat bins input, (1) extract the lines containing the desired repeat class from the repeat GFF file for each species, (2) convert the standard chromosome naming into src notation (i.e. species.chromosome) for these extracted lines, (3) concat these lines together, and (4) run gff2bins.py to bin the repeats.

### Example repeat GFF of L1 repeats used as input to gff2bins.py

```
Xtropicalis.Chr1	RepeatMasker	similarity	72760	73981	12.1	-	.	Target	"Motif:L1-49_XT"	4493	5716;class	"L1"
XlaevisL.Chr1	RepeatMasker	similarity	290962	291097	27.1	-	.	Target	"Motif:L1-26_XT"	5361	5508
XlaevisS.Chr1	RepeatMasker	similarity	16135	16460	25.1	+	.	Target	"Motif:L1-56_XT"	4616	4958
```
