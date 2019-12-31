#!/bin/bash

BLASTN=;
WORKDIR=`pwd`;

if [[ ! -f $1 || ! -f $2 ]]; then
  echo "Usage: `basename $0` <assembly.fa> <mitodb>"
  echo "Pipeline to flag contaminant scaffolds using the mitochrondrial assembly blast database."
  echo "Requires blastn, which can be hard set or included in the PATH."
  exit 1
fi

if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
    BLASTN=`which blastn`;
    if [[ -z "${BLASTN}" || ! -x "${BLASTN}" ]]; then
        echo "BLASTN not found or not executable";
	exit 127;
    fi
fi

$BLASTN -evalue 0.0000000001 -query $1 -db $2 -task blastn -outfmt 6 -perc_identity 99 -num_threads 35 1>$WORKDIR/mt_blast.out 2>$WORKDIR/mt_blast.log
cut -f1 $WORKDIR/mt_blast.out | sort | uniq | while read f; do grep ${f} $WORKDIR/mt_blast.out | cut -f1,3-4,7-12 && printf "  " && grep -P "^${f}\t" $1.fai | cut -f1-2; done >$WORKDIR/mt_blast.out.filter
