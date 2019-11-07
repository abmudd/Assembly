#!/bin/bash

# Copyright (c)2017. The Regents of the University of California (Regents).
# All Rights Reserved. Permission to use, copy, modify, and distribute this
# software and its documentation for educational, research, and
# not-for-profit purposes, without fee and without a signed licensing
# agreement, is hereby granted, provided that the above copyright notice,
# this paragraph and the following two paragraphs appear in all copies,
# modifications, and distributions. Contact the Office of Technology
# Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA
# 94720-1620, (510) 643-7201, for commercial licensing opportunities.

# IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
# ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
# REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

# REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
# HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


###############################################################################
# Setup
###############################################################################

# Set hard variables
# ============================================================

balanceFasta=`dirname $0`/bin/balanceFasta
bam2bed=`dirname $0`/bin/bam2bed.py
blast_filter=`dirname $0`/bin/blastfilter.py
combine_2less_1less=`dirname $0`/bin/combine_2less_1less.py
keep_scaffolds=`dirname $0`/bin/keep_scaffolds.py
jellyfish_filter=`dirname $0`/bin/jellyfishfilter.py
medavedepth=`dirname $0`/bin/medavedepth.py
plotr=`dirname $0`/bin/plot.R
remove_scaffolds=`dirname $0`/bin/remove_scaffolds.py
repairbam=`dirname $0`/bin/repairbam.py
sspace_config=`dirname $0`/bin/sspace_config.py
sspace_tab=`dirname $0`/bin/sspace_tab.py
BEDTOOLS=;
BLASTN=;
BWA=;
JELLYFISH=;
HALFDEPTH=0;
KMER=60;
MAKEBLASTDB=;
#MISASSEMBLYSIZE=1000;
QBATCH=;
RSCRIPT=;
SAMTOOLS=;
SSPACE=;
WORKDIR=`pwd`;
VERSION='1.0'

printf "[%s] " `basename $0` >&2;
echo $@ >&2;


# Log function
# ============================================================

function log () {
    printf "[%s] " `basename $0` >&2;
    echo "$1." >&2;
}


# Error function
# ============================================================

function error () {
    printf "[%s] ERROR: " `basename $0` >&2;
    echo "$2." >&2;
    exit $1;
}


# Help usage function
# ============================================================

function usage () {
    printf "\n" >&2;
    printf "%s v%s \n" `basename $0` $VERSION >&2;
    printf "\n" >&2;
    printf "Usage: %s --config STR --fasta STR --mode STR [--halfdepth INT]\n" `basename $0` >&2;
    printf "       [--kminusone INT] [--meraculous] [--module STR] [--workdir STR] [--help]\n" >&2;
    printf "\n" >&2;
    printf "Pipeline to remove duplicate haplotypes and then rescaffold.\n" >&2;
    printf "Intended for use with adapter trimmed Illumina or 10X data.\n" >&2;
#    printf "Inclusion of PB data is in development and not recommended for current use.\n" >&2;
    printf "\n" >&2;
    printf "Required arguments:\n" >&2;
    printf "       -c, --config STR        path for config file\n" >&2;
    printf "       -f, --fasta STR         path for assembly\n" >&2;
    printf "       -o, --mode STR          designate mode: prep, repropair, haplotype, scaffold,\n" >&2;
    printf "                               or finish\n" >&2;
    printf "\n" >&2;
    printf "Optional arguments:\n" >&2;
    printf "       -a, --halfdepth INT     max cutoff for median coverage of half depth [0]\n" >&2;
    printf "       -k, --kminusone INT     kmer size used in assembly minus one [60]\n" >&2;
    printf "       -e, --meraculous        config file is from meraculous [False]\n" >&2;
#    printf "       -i, --misasize INT      min PB read or 10X molecule length for misassembly [1000]\n" >&2;
    printf "       -m, --module STR        module load STR (such as bedtools/2.17.0)\n" >&2;
    printf "       -w, --workdir STR       path for working directory [pwd]\n" >&2;
    printf "       -h, --help              show this help message and exit\n" >&2;
    printf "\n" >&2;
    printf "Notes:\n" >&2;
    printf "   1. This script requires bedtools, blast, bwa, qbatch, samtools, SSPACE, and R. These\n" >&2;
    printf "      can be loaded using the module flag, manually hard set at the top of the script, or\n" >&2;
    printf "      included in the PATH.\n" >&2
    printf "   2. If not a meraculous config file, the config file must be a txt file with each\n" >&2
    printf "      line containing a library name; read 1 full path; read 2 full path; [FR,RF,10X,PB];\n" >&2
    printf "      insert size; and insert standard deviation; all split by tabs. 10X reads are\n" >&2;
    printf "      assumed to be in FR orientation.\n" >&2
#    printf "      For PB reads, leave blank tabs in place of the\n" >&2
#    printf "      read 2 full path; insert size; and insert standard deviation.\n" >&2
#    printf "   3. The user must manually align the PB data to the assembly and put the bam file in\n" >&2;
#    printf "      the properbam folder with the name library_name.bam. We recommend aligning using\n" >&2;
#    printf "      the map4cns pipeline (https://bitbucket.org/rokhsar-lab/map4cns).\n" >&2;
    printf "   3. If you are unsure of the insert size and standard deviation, please run the script\n" >&2
    printf "      estimate_insert_size.pl included with SSPACE.\n" >&2
    printf "\n" >&2;
    printf "Instructions:\n" >&2;
    printf "   1. Start by running prep mode.\n" >&2;
    printf "   2. If you want to change the insert size and standard deviation and rerun proper \n" >&2
    printf "      pairing, alter the input config file and then run repropair mode. Proper pairing\n" >&2
    printf "      is used in misassembly mode.\n" >&2;
    printf "   3. Run haplotype mode if sufficient half-depth scaffolds are present in the *.hist.pdf\n" >&2
    printf "      output.\n" >&2;
#    printf "   4. Run misassembly mode to identify and mask potential misassemblies.\n" >&2;
    printf "   4. Run scaffold mode to scaffold using SSPACE. Will only run if previous haplotype\n" >&2
    printf "      mode was run.\n" >&2;
    printf "   5. Run finish mode to clean up the temporary files.\n" >&2;
    printf "\n" >&2;
}


## Retrieve options on the command line and check for errors
# ============================================================

CONFIG=;
FASTA=;
HELP_MESSAGE=;
MERACULOUS=;
MODE=;

while [[ -n $@ ]]; do
    case "$1" in
        '--config') shift; CONFIG=$1;;
        '--fasta') shift; FASTA=$1;;
        '--halfdepth') shift; HALFDEPTH=$1;;
        '--kmer') shift; KMER=$1;;
#        '--misasize') shift; MISASSEMBLYSIZE=$1;;
        '--mode') shift; MODE=$1;;
        '--module') shift; module load $1;;
        '--meraculous') MERACULOUS=1;;
        '--workdir') shift; WORKDIR=$1;;
        '--help') HELP_MESSAGE=1;;
        '-c') shift; CONFIG=$1;;
        '-f') shift; FASTA=$1;;
        '-a') shift; HALFDEPTH=$1;;
        '-k') shift; KMER=$1;;
#        '-i') shift; MISASSEMBLYSIZE=$1;;
        '-o') shift; MODE=$1;;
        '-m') shift; module load $1;;
        '-e') MERACULOUS=1;;
        '-w') shift; WORKDIR=$1;;
        '-h') HELP_MESSAGE=1;;
        -*) usage; error 2 "Invalid option: ${1}";;
        *) break;;
    esac;
    shift;
done

if [[ -n "${HELP_MESSAGE}" ]]; then
    usage;
    exit 1;

elif [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
    usage;
    error 1 "CONFIG not defined or does not exist";

elif [[ -z "${FASTA}" || ! -f "${FASTA}" ]]; then
    usage;
    error 1 "FASTA not defined or does not exist";

elif [[ -z "${WORKDIR}" || ! -d "${WORKDIR}" ]]; then
    usage;
    error 1 "WORKDIR not defined or does not exist";

elif [[ -z "${HALFDEPTH}" ]]; then
    usage;
    error 1 "HALFDEPTH not defined";

elif [[ -z "${KMER}" ]]; then
    usage;
    error 1 "KMER not defined";

#elif [[ -z "${MISASSEMBLYSIZE}" ]]; then
#    usage;
#    error 1 "MISASSEMBLYSIZE not defined";

elif [[ -z "${MODE}" ]]; then
    usage;
    error 1 "MODE not defined";

elif [[ -z "${balanceFasta}" || ! -f "${balanceFasta}" ]]; then
    error 1 "balanceFasta not defined or does not exist";

elif [[ -z "${bam2bed}" || ! -f "${bam2bed}" ]]; then
    error 1 "bam2bed not defined or does not exist";

elif [[ -z "${blast_filter}" || ! -f "${blast_filter}" ]]; then
    error 1 "blast_filter not defined or does not exist";

elif [[ -z "${combine_2less_1less}" || ! -f "${combine_2less_1less}" ]]; then
    error 1 "combine_2less_1less not defined or does not exist";

elif [[ -z "${keep_scaffolds}" || ! -f "${keep_scaffolds}" ]]; then
    error 1 "keep_scaffolds not defined or does not exist";

elif [[ -z "${jellyfish_filter}" || ! -f "${jellyfish_filter}" ]]; then
    error 1 "jellyfish_filter not defined or does not exist";

elif [[ -z "${medavedepth}" || ! -f "${medavedepth}" ]]; then
    error 1 "medavedepth not defined or does not exist";

elif [[ -z "${plotr}" || ! -f "${plotr}" ]]; then
    error 1 "plotr not defined or does not exist";

elif [[ -z "${remove_scaffolds}" || ! -f "${remove_scaffolds}" ]]; then
    error 1 "remove_scaffolds not defined or does not exist";

elif [[ -z "${repairbam}" || ! -f "${repairbam}" ]]; then
    error 1 "repairbam not defined or does not exist";

elif [[ -z "${sspace_config}" || ! -f "${sspace_config}" ]]; then
    error 1 "sspace_config not defined or does not exist";

elif [[ -z "${sspace_tab}" || ! -f "${sspace_tab}" ]]; then
    error 1 "sspace_tab not defined or does not exist";

fi

MODE=$(tr '[A-Z]' '[a-z]' <<< $MODE)
ORIGCONFIG=$CONFIG


# Check that external tools are accessible
# ============================================================

if [[ -z "${BEDTOOLS}" || ! -x "${BEDTOOLS}" ]]; then
    BEDTOOLS=`which bedtools`;
    if [[ -z "${BEDTOOLS}" || ! -x "${BEDTOOLS}" ]]; then
        error 127 "BEDTOOLS not found or not executable";
    fi
fi

if [[ -z "${BWA}" || ! -x "${BWA}" ]]; then
    BWA=`which bwa`;
    if [[ -z "${BWA}" || ! -x "${BWA}" ]]; then
        error 127 "BWA not found or not executable";
    fi
fi

if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
    BLASTN=`which blastn`;
    if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
        error 127 "BLASTN not found or not executable";
    fi
fi

if [[ -z "${JELLYFISH}" || ! -x "${JELLYFISH}" ]]; then
    JELLYFISH=`which jellyfish`;
    if [[ -z "${JELLYFISH}" || ! -x "${JELLYFISH}" ]]; then
        error 127 "JELLYFISH not found or not executable";
    fi
fi

if [[ -z "${MAKEBLASTDB}" || ! -x "${MAKEBLASTDB}" ]]; then
    MAKEBLASTDB=`which makeblastdb`;
    if [[ -z "${MAKEBLASTDB}" || ! -x "${MAKEBLASTDB}" ]]; then
        error 127 "MAKEBLASTDB not found or not executable";
    fi
fi

if [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
    QBATCH=`which qbatch`;
    if [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
        error 127 "QBATCH not found or not executable";
    fi
fi

if [[ -z "${SAMTOOLS}" || ! -x "${SAMTOOLS}" ]]; then
    SAMTOOLS=`which samtools`;
    if [[ -z "${SAMTOOLS}" || ! -x "${SAMTOOLS}" ]]; then
        error 127 "SAMTOOLS not found or not executable";
    fi
fi

if [[ -z "${SSPACE}" || ! -x "${SSPACE}" ]]; then
    SSPACE=`which SSPACE_Standard_v3.0.pl`;
    if [[ -z "${SSPACE}" || ! -x "${SSPACE}" ]]; then
        error 127 "SSPACE_Standard_v3.0.pl not found or not executable";
    fi
fi

if [[ -z "${RSCRIPT}" || ! -x "${RSCRIPT}" ]]; then
    RSCRIPT=`which Rscript`;
    if [[ -z "${RSCRIPT}" || ! -x "${RSCRIPT}" ]]; then
        error 127 "RSCRIPT not found or not executable";
    fi
fi


###############################################################################
# Run:
###############################################################################

# Write out initial command info and remove done files
# ============================================================ 

log "Running mode $MODE .."

# For rerun of proper pair
if [[ $MODE = "repropair" ]]; then
    if [[ ! -f "${WORKDIR}/log/prep.5.repairbam.done" || ! -f "${WORKDIR}/log/prep.6.samtools_merge.done" ]]; then
	error 1 "Must run prep before repropair"
    else
	rm ${WORKDIR}/log/prep.5.repairbam.done
	rm ${WORKDIR}/log/prep.6.samtools_merge.done
    fi
fi


# Make bwa index
# ============================================================

if [[ ! -f "${FASTA}.sa" ]]; then
    $BWA index $FASTA
fi


# Convert config file
# ============================================================

TENX=;
WGS=;
PB=;

# Convert meraculous config
if [[ -n "${MERACULOUS}" ]]; then
    grep '^lib_seq' $CONFIG | grep ',' | sed 's/\,/\ /' | awk 'OFS="\t" { if ( $9 == "0" ) { print $4,$2,$3,"FR",$5,$6 } else { print $4,$2,$3,"RF",$5,$6 } }' | sort >${WORKDIR}/convert.config
    CONFIG=${WORKDIR}/convert.config
    WGS=1
    log "Converted Meraculous config file"

# Set variables if WGS and/or 10X libraries
elif grep -Pq '\t10X\t' $CONFIG || grep -Pq '\tPB\t' $CONFIG; then
    if [[ $(grep -Pc '\t10X\t' $CONFIG) == $(wc -l <$CONFIG) ]]; then
	TENX=1
	log "Identified only 10X libraries"
    elif [[ $(grep -Pc '\tPB\t' $CONFIG) == $(wc -l <$CONFIG) ]]; then
	PB=1
	log "Identified only PB libraries"
    elif [[ $(awk '{ if ( $4 == "10X" || $4 == "PB" ) {count++}} END { print count }' $CONFIG) == $(wc -l <$CONFIG) ]]; then
	PB=1
	TENX=1
	log "Identified PB and 10X libraries"
    elif [[ $(grep -Pc '\tPB\t' $CONFIG) == 0 ]]; then
	TENX=1
	WGS=1
	log "Identified WGS and 10X libraries"
    elif [[ $(grep -Pc '\t10X\t' $CONFIG) == 0 ]]; then
	PB=1
	WGS=1
	log "Identified WGS and PB libraries"
    else
	PB=1
	TENX=1
	WGS=1
	log "Identified WGS, PB, and 10X libraries"
    fi
    grep -Pv '\tPB\t' $CONFIG | sed 's/\t10X/\tFR/' >${WORKDIR}/convert.config
    CONFIG=${WORKDIR}/convert.config
else
    WGS=1
fi


# Run bwa aln/samse, samtools merge/sort/fixmate, repairbam, and medavedepth
# ============================================================

if [[ $MODE = "prep" ]]; then

    # Make directories if needed
    mkdir -p ${WORKDIR}/batch
    mkdir -p ${WORKDIR}/properbam
    mkdir -p ${WORKDIR}/log
    mkdir -p ${WORKDIR}/temp

    # Align each end with bwa aln
    if [[ -f ${WORKDIR}/log/prep.1.bwa_aln.done ]]; then
	log "BWA aln previously completed"
    else
	cat $CONFIG | while read name read1 read2 orient insave insstd; do echo "$BWA aln -t 32 $FASTA $read1 > ${WORKDIR}/temp/${name}.1.sai" && echo "$BWA aln -t 32 $FASTA $read2 > ${WORKDIR}/temp/${name}.2.sai"; done >${WORKDIR}/batch/bwa_aln.batch
	$QBATCH submit -W -p 16 -t 120:0:0 ${WORKDIR}/batch/bwa_aln.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/bwa_aln.batch.e* | wc -l) != $(wc -l <${WORKDIR}/batch/bwa_aln.batch) ]]; then
	    error 127 "BWA aln was not successful"
	else
	    touch ${WORKDIR}/log/prep.1.bwa_aln.done
	    log "BWA aln was successful"
	fi
    fi

    # Align each end with bwa samse and filter for mapQ >=30 and not flag 3844
    if [[ -f ${WORKDIR}/log/prep.2.bwa_samse.done ]]; then
	log "BWA samse previously completed"
    else
	cat $CONFIG | while read name read1 read2 orient insave insstd; do echo "$BWA samse ${FASTA} ${WORKDIR}/temp/${name}.1.sai $read1 | $SAMTOOLS view -@ 8 -bS -q 30 -F 3844 - >${WORKDIR}/temp/${name}.1.bam" && echo "$BWA samse ${FASTA} ${WORKDIR}/temp/${name}.2.sai $read2 | $SAMTOOLS view -@ 8 -bS -q 30 -F 3844 - >${WORKDIR}/temp/${name}.2.bam"; done >${WORKDIR}/batch/bwa_samse.batch
	$QBATCH submit -W -p 8 -t 12:0:0 ${WORKDIR}/batch/bwa_samse.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/bwa_samse.batch.e* | wc -l) != $(wc -l <${WORKDIR}/batch/bwa_samse.batch) ]]; then
	    error 127 "BWA samse was not successful"
	else
	    touch ${WORKDIR}/log/prep.2.bwa_samse.done
	    rm ${WORKDIR}/temp/*sai
	    log "BWA samse was successful"
	fi
    fi

    # Merge and name sort bam files for each library
    if [[ -f ${WORKDIR}/log/prep.3.bwa_mrgsrt.done ]]; then
	log "SAMTOOLS merge and sort previously completed"
    else
	cat $CONFIG | while read name read1 read2 orient insave insstd; do echo "$SAMTOOLS merge -@ 8 - ${WORKDIR}/temp/${name}.1.bam ${WORKDIR}/temp/${name}.2.bam | $SAMTOOLS sort -@ 8 -n - >${WORKDIR}/temp/${name}.mrgsrt.bam"; done >${WORKDIR}/batch/samtools_mrgsrt.batch
	$QBATCH submit -W -p 8 -t 12:0:0 ${WORKDIR}/batch/samtools_mrgsrt.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/samtools_mrgsrt.batch.e* | wc -l) != $(wc -l <${WORKDIR}/batch/samtools_mrgsrt.batch) ]]; then
	    error 127 "SAMTOOLS merge and sort were not successful"
	else
	    touch ${WORKDIR}/log/prep.3.bwa_mrgsrt.done
	    log "SAMTOOLS merge and sort were successful"
	fi
    fi

    # Fix mate for each library
    if [[ -f ${WORKDIR}/log/prep.4.samtools_fix.done ]]; then
	log "SAMTOOLS fixmate previously completed"
    else
	cat $CONFIG | while read name read1 read2 orient insave insstd; do echo "$SAMTOOLS fixmate -rp ${WORKDIR}/temp/${name}.mrgsrt.bam ${WORKDIR}/temp/${name}.mrgsrt.fix.bam"; done >${WORKDIR}/batch/samtools_fix.batch
	$QBATCH submit -T 20 ${WORKDIR}/batch/samtools_fix.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/samtools_fix.batch.e* | wc -l) != $(wc -l <${WORKDIR}/batch/samtools_fix.batch) ]]; then
	    error 127 "SAMTOOLS fixmate was not successful"
	else
	    touch ${WORKDIR}/log/prep.4.samtools_fix.done
	    log "SAMTOOLS fixmate was successful"
	fi
    fi

    # Determine proper pair for each library and output insert size and standard deviation to new.config
    if [[ -f ${WORKDIR}/log/prep.5.repairbam.done ]]; then
	log "REPAIRBAM previously completed"
    else
	cat $CONFIG | while read name read1 read2 orient insave insstd; do printf "${name}\t" && $SAMTOOLS view ${WORKDIR}/temp/${name}.mrgsrt.fix.bam | head -n1 | wc -l; done | grep -v '1$' >${WORKDIR}/temp/count
	if [[ $(wc -l <${WORKDIR}/temp/count) != 0 ]]; then
	    error 127 "Empty bam file found - remove library from config and resume"
	fi
	cat $CONFIG | while read name read1 read2 orient insave insstd; do echo "$repairbam -i $insave -s $insstd -p $orient ${WORKDIR}/temp/${name}.mrgsrt.fix.bam $CONFIG ${WORKDIR}/properbam/${name} && $RSCRIPT $plotr hist ${WORKDIR}/properbam/${name}.insert 1 50 ${WORKDIR}/properbam/${name}.insert"; done >${WORKDIR}/batch/repairbam.batch
	$QBATCH submit -T 20 ${WORKDIR}/batch/repairbam.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/repairbam.batch.e* | wc -l) != $(wc -l <${WORKDIR}/batch/repairbam.batch) ]]; then
	    error 127 "REPAIRBAM was not successful"
	else
	    cat ${WORKDIR}/properbam/*config >${WORKDIR}/new.config
	    CONFIG=${WORKDIR}/new.config
            touch ${WORKDIR}/log/prep.5.repairbam.done
	    log "REPAIRBAM was successful"
	fi
    fi

    # Merge, sort, and index bam files
    if [[ -f ${WORKDIR}/log/prep.6.samtools_merge.done ]]; then
	log "SAMTOOLS merge, sort, and index previously completed"
    else
	if [[ -n $PB ]]; then
	    log "PB bam file must be in ${WORKDIR}/properbam/ or you must rerun samtools merge and medave"
	fi
	$SAMTOOLS merge -@ 20 ${WORKDIR}/properbam/all.bam ${WORKDIR}/properbam/*.bam
	$SAMTOOLS sort -@ 20 -m 1G -o ${WORKDIR}/properbam/all.sort.bam ${WORKDIR}/properbam/all.bam
	$SAMTOOLS index ${WORKDIR}/properbam/all.sort.bam
	if [[ -f ${WORKDIR}/properbam/all.sort.bam.bai ]]; then
	    touch ${WORKDIR}/log/prep.6.samtools_merge.done
	    log "SAMTOOLS merge, sort, and index were successful"
	else
	    error 127 "SAMTOOLS merge, sort, and index were not successful"
	fi
    fi

    # Split fasta files and run medavedepth
    if [[ -f ${WORKDIR}/log/prep.7.medave.done ]]; then
	log "MEDAVE previously completed"
    else
	$balanceFasta -p 100 $FASTA ${WORKDIR}/temp/sequence
	for n in {0001..0100}; do echo "$medavedepth -b ${WORKDIR}/temp/sequence.0${n}.bed ${WORKDIR}/properbam/all.sort.bam $FASTA >${WORKDIR}/temp/sequence.0${n}.medave"; done >${WORKDIR}/batch/medave.batch
	$QBATCH submit -W -R 20 -t 12:0:0 ${WORKDIR}/batch/medave.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/medave.batch.e* | wc -l) != $(wc -l <${WORKDIR}/batch/medave.batch) ]]; then
	    error 127 "MEDAVE was not successful"
	else
	    # Plot histograms of median and average output as well as cut <= 100 output
	    cat ${WORKDIR}/temp/sequence.*.medave >${WORKDIR}/medave
	    awk '{if ($3 <= 100 && $4 <= 100 && $3 > 0 && $4 > 0) {print}}' ${WORKDIR}/medave >${WORKDIR}/medave.100
	    perl -ane 'BEGIN{%hist=()}{$hist{$F[2]}+=$F[1]} END{map{print("$_\t$hist{$_}\n")} sort {$a<=>$b} keys(%hist)}' ${WORKDIR}/medave.100 >${WORKDIR}/median.100.weighted
	    perl -ane 'BEGIN{%hist=()}{$hist{$F[3]}+=$F[1]} END{map{print("$_\t$hist{$_}\n")} sort {$a<=>$b} keys(%hist)}' ${WORKDIR}/medave.100 >${WORKDIR}/average.100.weighted
	    $RSCRIPT $plotr hist ${WORKDIR}/medave 3 100 ${WORKDIR}/median
	    $RSCRIPT $plotr hist ${WORKDIR}/medave 4 100 ${WORKDIR}/average
	    $RSCRIPT $plotr hist ${WORKDIR}/medave.100 3 100 ${WORKDIR}/median.100
	    $RSCRIPT $plotr hist ${WORKDIR}/medave.100 4 100 ${WORKDIR}/average.100
	    $RSCRIPT $plotr line ${WORKDIR}/median.100.weighted 1 2
	    $RSCRIPT $plotr line ${WORKDIR}/average.100.weighted 1 2
	    touch ${WORKDIR}/log/prep.7.medave.done
	    log "MEDAVE was successful"
	fi
    fi
fi


# Run single haplotype filter
# ============================================================

if [[ $MODE = "haplotype" ]]; then
    # Confirm that median and average depths were run successfully
    if [[ ! -f ${WORKDIR}/log/prep.7.medave.done ]]; then
	error 1 "Must run prep before haplotype"
    elif [[ -f ${WORKDIR}/log/hap.1.done ]]; then
	error 1 "Haplotype previously completed"
    fi

    # Extract scaffolds less than HALFDEPTH and convert to contigs
    mkdir -p ${WORKDIR}/haplotype
    awk -v var=$HALFDEPTH '{if ($3<=var) {print $1}}' ${WORKDIR}/medave | tee ${WORKDIR}/haplotype/halfdepth.list | $keep_scaffolds - $FASTA | perl -e 'my ($name,$seq)=("","");while(<>){chomp;if ((/\>/)||(eof)){if (($seq ne "")||(eof)){$seq.=$_ if (eof);my @Seq=split /N+/,$seq;if (@Seq>1){for my $i(1..@Seq){print "${name}_contig$i\n";print $Seq[$i-1],"\n";}}else{print "$name\n$seq\n";}}$name=$_;$seq="";}else{chomp;$seq.=$_;}}' >${WORKDIR}/haplotype/halfdepth.fa

    # Make blastdb and samtools faidx on contigs
    $MAKEBLASTDB -dbtype nucl -in ${WORKDIR}/haplotype/halfdepth.fa
    $SAMTOOLS faidx ${WORKDIR}/haplotype/halfdepth.fa

    # Blast all vs all contigs and remove those with a larger blast match
    $BLASTN -query ${WORKDIR}/haplotype/halfdepth.fa -db ${WORKDIR}/haplotype/halfdepth.fa -task blastn -outfmt 6 -num_threads 20 -word_size $KMER | awk '{if ($1!=$2) {print}}' | $blast_filter - ${WORKDIR}/haplotype/halfdepth.fa >${WORKDIR}/haplotype/blast.rm.list
    $JELLYFISH count -C -m 31 -s 10000000 -L 2 -U 100 -t 20 -o ${WORKDIR}/haplotype/halfdepth.mer ${WORKDIR}/haplotype/halfdepth.fa
    if [[ -f ${WORKDIR}/haplotype/halfdepth.mer_00 ]]; then
	$JELLYFISH merge -o ${WORKDIR}/haplotype/halfdepth.mer.jf ${WORKDIR}/haplotype/halfdepth.mer_*
	$JELLYFISH dump -o ${WORKDIR}/haplotype/halfdepth.mer.dump ${WORKDIR}/haplotype/halfdepth.mer.jf
    else
	$JELLYFISH dump -o ${WORKDIR}/haplotype/halfdepth.mer.dump ${WORKDIR}/haplotype/halfdepth.mer
    fi
    grep -v '^>' ${WORKDIR}/haplotype/halfdepth.mer.dump | awk '{print ">Contig"(++i)"\n" $0}' | $BLASTN -db ${WORKDIR}/haplotype/halfdepth.fa -task blastn -outfmt 6 -num_threads 20 -word_size 31 -perc_identity 100 | $jellyfish_filter - ${WORKDIR}/haplotype/halfdepth.fa >${WORKDIR}/haplotype/jellyfish.rm.list
    cat ${WORKDIR}/haplotype/jellyfish.rm.list ${WORKDIR}/haplotype/blast.rm.list | sort | uniq -d | $remove_scaffolds - ${WORKDIR}/haplotype/halfdepth.fa >${WORKDIR}/haplotype/singlehap.fa
    touch ${WORKDIR}/log/hap.1.done
fi


# Run misassembly filter
# ============================================================

#if [[ $MODE = "misassembly" ]]; then
    # Confirm that median and average depths were run successfully and that data contains 10X or PB
#    if [[ ! -f ${WORKDIR}/log/prep.7.medave.done ]]; then
#	error 1 "Must run prep before misassembly"
#    elif [[ ! -n $PB ]]; then
#	error 1 "Can only run misassembly with PB data"
#    elif [[ -f ${WORKDIR}/log/mis.1.done ]]; then
#	error 1 "Misassembly previously completed"
#    fi

    # Filter libraries to 10X or PB
#    mkdir -p ${WORKDIR}/misassembly
#    awk '{if ($4=="PB" || $4=="10X") {print $1 "\t" $4}}' $ORIGCONFIG >${WORKDIR}/misassembly/libraries
#    cat ${WORKDIR}/misassembly/libraries | while read library orient; do $SAMTOOLS sort -@ 20 -m 1G -o ${WORKDIR}/properbam/${library}.sort.bam ${WORKDIR}/properbam/${library}.bam && $SAMTOOLS index ${WORKDIR}/properbam/${library}.sort.bam; done

    # If single haplotype filter was run, remove scaffolds from that analysis
#    if [[ -f ${WORKDIR}/haplotype/singlehap.fa ]]; then
#	$remove_scaffolds ${WORKDIR}/haplotype/halfdepth.list $FASTA > ${WORKDIR}/misassembly/exclhalfdepth.fa
#	$SAMTOOLS faidx ${WORKDIR}/misassembly/exclhalfdepth.fa
#	awk '{print $1 "\t0\t" $2}' ${WORKDIR}/misassembly/exclhalfdepth.fa.fai >${WORKDIR}/misassembly/filter.bed
#	FASTA=${WORKDIR}/misassembly/exclhalfdepth.fa

    # If single haplotype filter wasn't run, use full fasta file
#    else
#	$SAMTOOLS faidx $FASTA
#	awk '{print $1 "\t0\t" $2}' ${FASTA}.fai >${WORKDIR}/misassembly/filter.bed
#    fi

    # Convert bam to bed
    # If PB
#    if [[ -n $PB ]]; then
#	grep -P '\tPB' ${WORKDIR}/misassembly/libraries | while read library orient; do $bam2bed -o ${WORKDIR}/misassembly/${library} -i ${WORKDIR}/properbam/${library}.sort.bam -t PacBio -l $MISASSEMBLYSIZE; done
#    fi
    # If 10X
#    if [[ -n $TENX ]]; then
#	grep -P '\t10X' ${WORKDIR}/misassembly/libraries | while read library orient; do $bam2bed -o ${WORKDIR}/misassembly/${library} -i ${WORKDIR}/properbam/${library}.sort.bam -t 10X -l $MISASSEMBLYSIZE; done
#    fi

    # Merge bed files and run tabix
#    cat ${WORKDIR}/misassembly/libraries | while read library orient; do zcat ${WORKDIR}/misassembly/${library}.bed.gz; done | $BEDTOOLS intersect -wa -a stdin -b ${WORKDIR}/misassembly/filter.bed >${WORKDIR}/misassembly/all.bed
#    python -c "import pysam; pysam.tabix_index('${WORKDIR}/misassembly/all.bed', force=True, preset='bed')"
#    log "BAM2BED was successful"

    # Determine regions with low coverage
#    $BEDTOOLS genomecov -d -i ${WORKDIR}/misassembly/all.bed.gz -g ${FASTA}.fai | awk '{if ($3<=2) {print $1 "\t" $2-1 "\t" $2}}' | $BEDTOOLS merge -i stdin >${WORKDIR}/misassembly/2orLess_all.depth.bed
#    $BEDTOOLS genomecov -d -i ${WORKDIR}/misassembly/all.bed.gz -g ${FASTA}.fai | awk '{if ($3<=1) {print $1 "\t" $2-1 "\t" $2}}' | $BEDTOOLS merge -i stdin >${WORKDIR}/misassembly/1orLess_all.depth.bed
#    $combine_2less_1less ${WORKDIR}/misassembly/2orLess_all.depth.bed ${WORKDIR}/misassembly/1orLess_all.depth.bed $FASTA ${WORKDIR}/misassembly/out
#    $BEDTOOLS maskfasta -fi $FASTA -bed ${WORKDIR}/misassembly/out.mask.bed -fo ${WORKDIR}/misassembly/final.masked.fa
#    $SAMTOOLS faidx ${WORKDIR}/misassembly/final.masked.fa
#    $BEDTOOLS sort -faidx ${WORKDIR}/misassembly/final.masked.fa.fai -i ${WORKDIR}/misassembly/out.trim.bed | $BEDTOOLS complement -i stdin -g ${WORKDIR}/misassembly/final.masked.fa.fai > ${WORKDIR}/misassembly/out.trim.complement.bed
#    $BEDTOOLS getfasta -fi ${WORKDIR}/misassembly/final.masked.fa -bed ${WORKDIR}/misassembly/out.trim.complement.bed -fo ${WORKDIR}/misassembly/final.masked_trim.fa
#    touch ${WORKDIR}/log/mis.1.done
#fi


# Scaffold assembly
# ============================================================

if [[ $MODE = "scaffold" ]]; then
    mkdir -p ${WORKDIR}/scaffold

    # Set fasta file if previously run haplotype and misassembly modes
    if [[ ! -f ${WORKDIR}/scaffold/input.fa ]]; then
	if [[ -f ${WORKDIR}/haplotype/singlehap.fa && -f ${WORKDIR}/misassembly/final.masked_trim.fa ]]; then
	    cat ${WORKDIR}/haplotype/singlehap.fa ${WORKDIR}/misassembly/final.masked_trim.fa > ${WORKDIR}/scaffold/input.fa
	elif [[ -f ${WORKDIR}/haplotype/singlehap.fa && ! -f ${WORKDIR}/misassembly/final.masked_trim.fa ]]; then
	    $remove_scaffolds ${WORKDIR}/haplotype/halfdepth.list $FASTA > ${WORKDIR}/scaffold/input.fa
	    cat ${WORKDIR}/haplotype/singlehap.fa >> ${WORKDIR}/scaffold/input.fa
	elif [[ ! -f ${WORKDIR}/haplotype/singlehap.fa && -f ${WORKDIR}/misassembly/final.masked_trim.fa ]]; then
	    cp ${WORKDIR}/misassembly/final.masked_trim.fa ${WORKDIR}/scaffold/input.fa
	else
	    cp $FASTA ${WORKDIR}/scaffold/input.fa
	    log "Unable to find output from single haplotype or misassembly"
	fi
    fi

    FASTA=${WORKDIR}/scaffold/input.fa

    # Check if SSPACE has previously run
    if [[ -f ${WORKDIR}/log/scaf.SSPACE.done ]]; then
	log "SSPACE previously completed"

    # Create SSPACE config and run SSPACE
    elif [[ -n $WGS && ! -f ${WORKDIR}/log/scaf.SSPACE.done ]]; then
	grep -Pv '\t10X\t' $ORIGCONFIG | grep -Pv '\tPB\t' | $sspace_config - ${WORKDIR}/scaffold/SSPACE.config
	echo "cd ${WORKDIR}/scaffold/ && perl $SSPACE -l SSPACE.config -s $FASTA -x 0 -n $KMER -k 10 -T 32 -b SSPACE" > ${WORKDIR}/batch/SSPACE.batch
	$QBATCH submit -W -p 32 -R 15 -t 96:0:0 ${WORKDIR}/batch/SSPACE.batch ${WORKDIR}/log

	if [[ $(grep done ${WORKDIR}/log/SSPACE.batch.e* | wc -l) != 1 ]]; then
	    error 127 "SSPACE was not successful"
	else
	    cp ${WORKDIR}/scaffold/SSPACE/SSPACE.final.scaffolds.fasta ${WORKDIR}/scaffold/scaffold.fa
	    touch ${WORKDIR}/log/scaf.SSPACE.done
	    log "SSPACE was successful"
	fi
    fi
fi


# Finish align pipeline
# ============================================================

if [[ $MODE = "finish" ]]; then
    mkdir -p ${WORKDIR}/pdf/
    mv ${WORKDIR}/*.pdf ${WORKDIR}/pdf/
    mv ${WORKDIR}/properbam/*.pdf ${WORKDIR}/pdf/
    rm ${WORKDIR}/properbam/*
    rm ${WORKDIR}/temp/*
    rm ${WORKDIR}/log/*.batch*
fi

log "Finished mode $MODE"

exit 0
