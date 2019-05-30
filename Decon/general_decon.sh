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

filter_decon_blast=`dirname $0`/bin/filter.decon_blast.py
ARCHAEADB=`dirname $0`/Archaea_RefSeq;
BACTERIADB=`dirname $0`/Bacteria_RefSeq;
BEDTOOLS=;
BLASTN=;
QBATCH=;
VECTORDB=`dirname $0`/UniVec;
VIRUSDB=`dirname $0`/Viral_RefSeq;
WORKDIR=`pwd`;
VERSION='1.0';


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
    printf "Usage: %s --fasta STR [--archaeadb STR] [--bacteriadb STR] [--module STR]\n" `basename $0` >&2;
    printf "       [--vectordb STR] [--virusdb STR] [--workdir STR] [--help]\n" >&2;
    printf "\n" >&2;
    printf "Pipeline to flag contaminant scaffolds using archaea, bacteria, vector, and virus databases.\n" >&2;
    printf "\n" >&2;
    printf "Required arguments:\n" >&2;
    printf "       -f, --fasta STR         path to assembly\n" >&2;
    printf "\n" >&2;
    printf "Optional arguments:\n" >&2;
    printf "       -a, --archaeadb STR     path to archaea blastdb [scriptsdir/Archaea_RefSeq]\n" >&2;
    printf "       -b, --bacteriadb STR    path to bacteria blastdb [scriptsdir/Bacteria_RefSeq]\n" >&2;
    printf "       -m, --module STR        module load STR (such as bedtools/2.27.0)\n" >&2;
    printf "       -v, --vectordb STR      path to virus blastdb [scriptsdir/UniVec]\n" >&2;
    printf "       -i, --virusdb STR       path to virus blastdb [scriptsdir/Viral_RefSeq]\n" >&2;
    printf "       -w, --workdir STR       path to working directory [pwd]\n" >&2;
    printf "       -h, --help              show this help message and exit\n" >&2;
    printf "\n" >&2;
    printf "Notes:\n" >&2;
    printf "   1. This script assumes a standard UNIX/Linux install and access to a SGE or SLURM job-\n" >&2;
    printf "      scheduling system via qbatch (see https://bitbucket.org/rokhsar-lab/gbs-analysis/src).\n" >&2;
    printf "   2. This script requires blastn, bedtools, and qbatch. These can be loaded using the\n" >&2;
    printf "      module flag, manually hard set at the top of the script, or included in the PATH.\n" >&2;
    printf "\n" >&2;
}


## Retrieve options on the command line and check for errors
# ============================================================

FASTA=;
HELP_MESSAGE=;

while [[ -n $@ ]]; do
    case "$1" in
        '--archaeadb') shift; ARCHAEADB=$1;;
        '--bacteriadb') shift; BACTERIADB=$1;;
        '--fasta') shift; FASTA=$1;;
        '--module') shift; module load $1;;
        '--vectordb') shift; VECTORDB=$1;;
        '--virusdb') shift; VIRUSDB=$1;;
        '--workdir') shift; WORKDIR=$1;;
        '--help') HELP_MESSAGE=1;;
        '-a') shift; ARCHAEADB=$1;;
        '-b') shift; BACTERIADB=$1;;
        '-f') shift; FASTA=$1;;
        '-m') shift; module load $1;;
        '-v') shift; VECTORDB=$1;;
        '-i') shift; VIRUSDB=$1;;
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

elif [[ -z "${FASTA}" || ! -f "${FASTA}" ]]; then
    error 1 "FASTA not defined or does not exist";

elif [[ -z "${ARCHAEADB}" ]]; then
    error 1 "ARCHAEADB not defined";

elif [[ -z "${BACTERIADB}" ]]; then
    error 1 "BACTERIADB not defined";

elif [[ -z "${VECTORDB}" ]]; then
    error 1 "VECTORDB not defined";

elif [[ -z "${VIRUSDB}" ]]; then
    error 1 "VIRUSDB not defined";

elif [[ -z "${filter_decon_blast}" || ! -f "${filter_decon_blast}" ]]; then
    error 1 "filter_decon_blast not defined or does not exist";

fi


# Check that external tools are accessible
# ============================================================

if [[ -z "${BEDTOOLS}" || ! -x "${BEDTOOLS}" ]]; then
    BEDTOOLS=`which bedtools`;
    if [[ -z "${BEDTOOLS}" || ! -x "${BEDTOOLS}" ]]; then
        error 127 "BEDTOOLS not found or not executable";
    fi
fi

if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
    BLASTN=`which blastn`;
    if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
        error 127 "BLASTN not found or not executable";
    fi
fi

if [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
    QBATCH=`which qbatch`;
    if [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
        error 127 "QBATCH not found or not executable";
    fi
fi


###############################################################################
# Run:
###############################################################################

# Make and submit blast scripts
# ============================================================

if [[ ! -f ${WORKDIR}/decon_blast.done ]]; then
    # Archaea
    echo "$BLASTN -evalue 0.0000000001 -query $FASTA -db $ARCHAEADB -task blastn -outfmt 6 -perc_identity 95 -num_threads 15 1>${WORKDIR}/Archaea_blast.out 2>${WORKDIR}/Archaea_blast.log" >${WORKDIR}/decon_blast.batch;
    # Vector
    echo "$BLASTN -evalue 0.0000000001 -query $FASTA -db $VECTORDB -task blastn -outfmt 6 -perc_identity 95 -num_threads 15 1>${WORKDIR}/Vector_blast.out 2>${WORKDIR}/Vector_blast.log" >>${WORKDIR}/decon_blast.batch;
    # Virus
    echo "$BLASTN -evalue 0.0000000001 -query $FASTA -db $VIRUSDB -task blastn -outfmt 6 -perc_identity 95 -num_threads 15 1>${WORKDIR}/Virus_blast.out 2>${WORKDIR}/Virus_blast.log" >>${WORKDIR}/decon_blast.batch;
    #Bacteria
    echo "$BLASTN -evalue 0.0000000001 -query $FASTA -db $BACTERIADB -task blastn -outfmt 6 -perc_identity 95 -num_threads 15 1>${WORKDIR}/Bacteria_blast.out 2>${WORKDIR}/Bacteria.${n}_blast.log" >>${WORKDIR}/decon_blast.batch;
    $QBATCH submit -W -p 8 -t 168:0:0 ${WORKDIR}/decon_blast.batch;

    if [[ $(grep done ${WORKDIR}/decon_blast.batch.log/decon_blast.batch.e* | wc -l) != $(wc -l < ${WORKDIR}/decon_blast.batch) ]]; then
	error 127 "Decon blast was not successful"
    else
	touch ${WORKDIR}/decon_blast.done
    fi
fi


# Parse all blast output
# ============================================================

for n in Archaea Vector Virus Bacteria; do $filter_decon_blast ${WORKDIR}/${n}_blast.out $FASTA >${WORKDIR}/${n}_blast.out.filter; done;

exit 0;
