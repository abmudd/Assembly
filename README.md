# Assembly

Pipelines and scripts related to genome assembly, annotation, and analysis. Brief summaries of each directory are listed below, and each directory has its own README with further information.

Version 1.0

## 24FrogAnalysis

Scripts related to the 24 frog analysis:
* countalnspecies.py
* filterancestralProbs.py
* maf2uce.py
* splitmafchr.py

## AlignRNAseq

STARalign.sh: Bash pipeline to align RNA-seq reads to a genome with splice junctions using STAR.

## AssembleOrganelle

organelle_pipeline.py: Bash pipeline to download organelle sequences from species close to the sequenced species using NCBI accessions, create the respective NOVOPlasty config files, and run NOVOPlasty.

## CactusWrapper

cactus_filter.py: Python pipeline to run cactus with two species and filter the output. Cactus jobs are submitted to a SLURM cluster with sbatch using the defined memory, processor, and time flags, whereas filter jobs are run locally. All cactus variables must either be in the appropriate paths before running or have the -l flag set to load cactus as a module.

## Decontamination

Scripts related to genome assembly decontamination:
* general_decon.sh
* mt_decon.sh
* nt_decon.sh

## Extract4DSites

4Dextract.py: Python script to extract bases with four fold degeneracy from input fasta and genePred annotation of reference species, informed by conserved amino acids in other species using input MAF alignment with the reference species on top.

## FilterAnnotation

Scripts related to genome annotation filtering:
* filter_trinity.py
* largestgenePred.py
* postGeMoMa.py

## FilterVCF

Scripts related to filtering VCF files for a genetic map:
* AD2loc.py
* F1call.py
* F2call.py

## IdentifyMisassemblies

LoReM.py: Python script to bin and identify regions of the genome assembly with low spanning coverage and high number of read terminals (starts and ends) using 10X linked reads, PacBio reads, or Nanopore reads. For 10X linked reads, the input bam is the sorted and indexed bam with a BX optional field, such as the output from longranger align. For PacBio or Nanopore reads, the input bam is the sorted and indexed bam from your choice of aligner.

## Miscellaneous

Miscellaneous Python scripts:
* extractmatch.py
* getgapbed.py
* removefromgff.py
* replacegffname.py
* replacescafname.py
* startendNs.py

## MuntjacAnalysis

Scripts related to the muntjac analysis:
* extract2speciesmaf.py
* extractOrthoVenn.py
* HiCbins_1Mb.py
* HiCbins_100kb.py
* mcscan_convert_links.py
* mcscan_invert_chr.py
* run_fission_fusion.sh

## PlotAlignments

alignment_plots.py: Python pipeline to writes out alignments of three species chromosomes with centromere locations and repeat histograms.

## RemoveDuplicateHaplotype

align_pipeline.sh: Bash pipeline to remove duplicate haplotypes and then rescaffold. Intended for use with adapter trimmed Illumina or 10X data.

## TrimReads

Scripts related to read trimming:
* nxtrim_pipeline.sh
* trim_10X.py

## Copyright

Copyright (c)2019. The Regents of the University of California (Regents). All Rights Reserved. Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact the Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
