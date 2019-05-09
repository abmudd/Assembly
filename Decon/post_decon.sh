#!/bin/bash
module load blast+/2.6.0
module load samtools/1.6
less [A,B,V,m]*filter
emacs remove.scaffolds
# cut -f1 mt_blast.out | sort | uniq >>remove.scaffolds
grep -v '^ ' approved_blast.unknown.scaffolds.filter | tr ':' '\t' | tr '-' '\t' | awk '{if ($3 < 1000 || ($5/$3) > 0.5) {print $1}}' >>remove.scaffolds
sort remove.scaffolds | uniq >temp && mv temp remove.scaffolds
~/Assembly/Decon/remove_scaffolds.py remove.scaffolds final.scaffolds.unfiltered.single-haplotype.fa >assembly.decon.fa
samtools faidx assembly.decon.fa
wc -l final.scaffolds.unfiltered.single-haplotype.fa.fai
wc -l remove.scaffolds
wc -l assembly.decon.fa.fai
