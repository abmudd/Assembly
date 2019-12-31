# 24 frog analysis scripts

Version 1.0

Citing: unpublished

## 1. splitmafchr.py

Python script that splits the input MAF file into unique MAF files for each reference scaffold and chromosome.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: splitmafchr.py <in.maf> <out.prefix>
This script splits the input MAF file into unique MAF files for each reference scaffold and chromosome.
```

## 2. countalnspecies.py

Python script that extracts the number of aligned species in each block of the input MAF file and outputs in either bed or wig format.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: countalnspecies.py <in.maf> <BED/WIG>
This script extracts the number of aligned species in each block of the input MAF file and outputs in either bed or wig format.
```

## 3. filterancestralProbs.py

Python script that extracts bases from input RAxML ancestralProbs of a particular node with a minimum probability and an unambiguous base.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: filterancestralProbs.py <in.ancestralProbs> <node> <min_p>
This script extracts bases from input RAxML ancestralProbs of a particular node with a minimum probability and an unambiguous base.
```

## 4. maf2uce.py

Python script that calculates the locations of UCEs and outputs the location as a bed file relative to species1. The script requires 1-to-1 orthologous alignment of the two species.

### Prerequisite Python modules:

```
os
sys
```

### Help message:

```
Usage: maf2uce.py <in.maf> <species1> <species2> <UCElength> <UCE%id>
This script calculates the locations of UCEs and outputs the location as a bed file relative to species1. The script requires 1-to-1 orthologous alignment of the two species.
```
