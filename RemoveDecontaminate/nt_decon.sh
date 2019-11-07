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

APPROVEDDB=;
APPROVEDSEQID=;
NTDB=`dirname $0`/nt
balanceFasta=`dirname $0`/bin/balanceFasta
filter_decon_blast=`dirname $0`/bin/filter.decon_blast.py
filter_nt_blast=`dirname $0`/bin/filter.nt_blast.py
keep_scaffolds=`dirname $0`/bin/keep_scaffolds.py
BEDTOOLS=;
BLASTDBCMD=;
BLASTN=;
QBATCH=;
SAMTOOLS=;
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
    printf "Usage: %s --approveddb STR --fasta STR --approvedseqid STR [--module STR]\n" `basename $0` >&2;
    printf "       [--ntdb STR] [--workdir STR] [--help]\n" >&2;
    printf "\n" >&2;
    printf "Pipeline to flag contaminant scaffolds using the nt database.\n" >&2;
    printf "\n" >&2;
    printf "Required arguments:\n" >&2;
    printf "       -a, --approveddb STR    path for approved species db\n" >&2;
    printf "       -f, --fasta STR         path for assembly\n" >&2;
    printf "       -s, --approvedseqid STR path for approved sequence ids\n" >&2;
    printf "\n" >&2;
    printf "Optional arguments:\n" >&2;
    printf "       -m, --module STR        module load STR (such as bedtools/2.27.0)\n" >&2;
    printf "       -n, --ntdb STR          path for nt database [scriptsdir/nt]\n" >&2;
    printf "       -w, --workdir STR       path for working directory [pwd]\n" >&2;
    printf "       -h, --help              show this help message and exit\n" >&2;
    printf "\n" >&2;
    printf "Notes:\n" >&2;
    printf "   1. This script assumes a standard UNIX/Linux install and access to a SGE or SLURM job-\n" >&2;
    printf "      scheduling system via qbatch (see https://bitbucket.org/rokhsar-lab/gbs-analysis/src).\n" >&2;
    printf "   2. This script requires blastn, bedtools, qbatch, and samtools. These can be loaded\n" >&2;
    printf "      using the module flag, manually hard set at the top of the script, or included in\n" >&2;
    printf "      the PATH.\n" >&2;
    printf "   3. The approved sequence ids should be the output from blastdbcmd -outfmt \"%i\" of\n" >&2;
    printf "      those sequences in the ntdb that are verified as close to the sequenced species.\n" >&2;
    printf "   4. The approved species db must be a txt file with each line containing a db name and\n" >&2;
    printf "      full path to the db location, split by a tab. These should be assemblies or fasta,\n" >&2;
    printf "      such as EST or GSS sequences, from closely related species that will be used as a\n" >&2;
    printf "      second pass for flagging contaminant sequences.\n" >&2;
    printf "\n" >&2;
}


## Retrieve options on the command line and check for errors
# ============================================================

FASTA=;
HELP_MESSAGE=;

while [[ -n $@ ]]; do
    case "$1" in
        '--approveddb') shift; APPROVEDDB=$1;;
        '--approvedseqid') shift; APPROVEDSEQID=$1;;
        '--ntdb') shift; NTDB=$1;;
        '--fasta') shift; FASTA=$1;;
        '--module') shift; module load $1;;
        '--workdir') shift; WORKDIR=$1;;
        '--help') HELP_MESSAGE=1;;
        '-a') shift; APPROVEDDB=$1;;
        '-s') shift; APPROVEDSEQID=$1;;
        '-n') shift; NTDB=$1;;
        '-f') shift; FASTA=$1;;
        '-m') shift; module load $1;;
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

elif [[ -z "${APPROVEDDB}" || ! -f "${APPROVEDDB}" ]]; then
    error 1 "APPROVEDDB not defined or does not exist";

elif [[ -z "${APPROVEDSEQID}" || ! -f "${APPROVEDSEQID}" ]]; then
    error 1 "APPROVEDSEQID not defined or does not exist";

elif [[ -z "${FASTA}" || ! -f "${FASTA}" ]]; then
    error 1 "FASTA not defined or does not exist";

elif [[ -z "${NTDB}" || ! -f "${NTDB}" ]]; then
    error 1 "NTDB not defined or does not exist";

elif [[ -z "${balanceFasta}" || ! -f "${balanceFasta}" ]]; then
    error 1 "balanceFasta not defined or does not exist";

elif [[ -z "${filter_decon_blast}" || ! -f "${filter_decon_blast}" ]]; then
    error 1 "filter_decon_blast not defined or does not exist";

elif [[ -z "${filter_nt_blast}" || ! -f "${filter_nt_blast}" ]]; then
    error 1 "filter_nt_blast not defined or does not exist";

elif [[ -z "${keep_scaffolds}" || ! -f "${keep_scaffolds}" ]]; then
    error 1 "keep_scaffolds not defined or does not exist";

fi


# Check that external tools are accessible
# ============================================================

if [[ -z "${BEDTOOLS}" || ! -x "${BEDTOOLS}" ]]; then
    BEDTOOLS=`which bedtools`;
    if [[ -z "${BEDTOOLS}" || ! -x "${BEDTOOLS}" ]]; then
        error 127 "BEDTOOLS not found or not executable";
    fi
fi

if [[ -z "${BLASTDBCMD}" || ! -x "${BLASTDBCMD}" ]]; then
    BLASTDBCMD=`which blastdbcmd`;
    if [[ -z "${BLASTDBCMD}" || ! -x "${BLASTDBCMD}" ]]; then
        error 127 "BLASTDBCMD not found or not executable";
    fi
fi

if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
    BLASTN=`which blastn`;
    if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
        error 127 "BLASTN not found or not executable";
    fi
fi

if [[ -z "${SAMTOOLS}" || ! -x "${SAMTOOLS}" ]]; then
    SAMTOOLS=`which samtools`;
    if [[ -z "${SAMTOOLS}" || ! -x "${SAMTOOLS}" ]]; then
        error 127 "SAMTOOLS not found or not executable";
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

# Make samtools faidx
# ============================================================

if [[ ! -f "${FASTA}.fai" ]]; then
    $SAMTOOLS faidx $FASTA;
fi


# Split fasta into 1000
# ============================================================

if [[ ! -f ${WORKDIR}/split.done ]]; then
    mkdir -p ${WORKDIR}/split;
    $balanceFasta -p 1000 $FASTA ${WORKDIR}/split/sequence;
    touch ${WORKDIR}/split.done
fi


# Extract each sequence and blast against nt
# ============================================================

if [[ ! -f ${WORKDIR}/nt_blast.done ]]; then
    for n in {0001..1000}; do echo "$BEDTOOLS getfasta -fi $FASTA -bed ${WORKDIR}/split/sequence.0${n}.bed -fo ${WORKDIR}/split/sequence.0${n}.fa && $BLASTN -evalue 0.0000000001 -query ${WORKDIR}/split/sequence.0${n}.fa -db $NTDB -outfmt 6 -perc_identity 75 -num_threads 8 1>${WORKDIR}/split/nt.0${n}.blast.out"; done >${WORKDIR}/nt_blast.batch
    $QBATCH submit -W -p 8 -t 24:0:0 ${WORKDIR}/nt_blast.batch

    if [[ $(grep done ${WORKDIR}/nt_blast.batch.log/nt_blast.batch.e* | wc -l) != $(wc -l < ${WORKDIR}/nt_blast.batch) ]]; then
	error 127 "NT blast was not successful"
    else
	cat ${WORKDIR}/split/nt.0*.blast.out >${WORKDIR}/nt_blast.out
	rm ${WORKDIR}/split/*
	touch ${WORKDIR}/nt_blast.done
    fi
fi


# Filter output of nt blast and extract hits without approved best match
# ============================================================

if [[ ! -f ${WORKDIR}/nt_blast.out.nomatch.fa.fai ]]; then
    $filter_nt_blast ${WORKDIR}/nt_blast.out $FASTA $APPROVEDSEQID >${WORKDIR}/nt_blast.out.filter
    grep -v '^#' ${WORKDIR}/nt_blast.out.filter | cut -f1 -d':' | $keep_scaffolds - $FASTA >${WORKDIR}/nt_blast.out.nomatch.fa
    $SAMTOOLS faidx ${WORKDIR}/nt_blast.out.nomatch.fa
fi


# Blast hits against approved db
# ============================================================

if [[ ! -f ${WORKDIR}/approved_blast.done ]]; then
    mkdir -p ${WORKDIR}/nt_blast.out.nomatch
    cat $APPROVEDDB | while read name file; do echo "$BLASTN -evalue 0.0000000001 -query ${WORKDIR}/nt_blast.out.nomatch.fa -db $file -outfmt 6 -perc_identity 75 -num_threads 8 1>${WORKDIR}/nt_blast.out.nomatch/${name}_blast.out 2>${WORKDIR}/nt_blast.out.nomatch/${name}_blast.log"; done >${WORKDIR}/approved_blast.batch
    $QBATCH submit -W -p 8 -t 72:0:0 ${WORKDIR}/approved_blast.batch

    if [[ $(grep done ${WORKDIR}/approved_blast.batch.log/approved_blast.batch.e* | wc -l) != $(wc -l < ${WORKDIR}/approved_blast.batch) ]]; then
	error 127 "Approved blast was not successful"
    else
	touch ${WORKDIR}/approved_blast.done
    fi
fi

echo "# List of nt no-match scaffolds with sufficient hit to the approved db" >${WORKDIR}/approved_blast.validated.scaffolds
cat $APPROVEDDB | while read name file; do cat ${WORKDIR}/nt_blast.out.nomatch/${name}_blast.out; done | cut -f1 | cut -f1 -d':' | sort | uniq >>${WORKDIR}/approved_blast.validated.scaffolds
echo "# List of nt no-match scaffolds without sufficient hit to the approved db" >${WORKDIR}/approved_blast.unknown.scaffolds
cat ${WORKDIR}/nt_blast.out.filter ${WORKDIR}/approved_blast.validated.scaffolds | grep -v '^#' | cut -f1 -d':' | sort | uniq -u >>${WORKDIR}/approved_blast.unknown.scaffolds
grep -v '^#' ${WORKDIR}/approved_blast.unknown.scaffolds | while read f; do grep "# ${f}:" ${WORKDIR}/nt_blast.out.filter | cut -c3- | cut -f1,3-4,11-12 && printf "  " && grep "# ${f}:" ${WORKDIR}/nt_blast.out.filter | cut -f2 | while read n; do $BLASTDBCMD -db $NTDB -entry "$n" -outfmt "%t"; done; done >${WORKDIR}/approved_blast.unknown.scaffolds.filter
grep -v '^ ' ${WORKDIR}/approved_blast.unknown.scaffolds.filter | tr ':' '\t' | tr '-' '\t' | awk '{if ($3 < 1000 || ($5/$3) > 0.5) {print $1}}' >${WORKDIR}/approved_blast.unknown.scaffolds.1kb_50perc.filter

exit 0
