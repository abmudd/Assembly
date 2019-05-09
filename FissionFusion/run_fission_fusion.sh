#!/bin/bash

###############################################################################
# Setup
###############################################################################

# Set hard variables
# ============================================================

Prefix=`dirname $0`
Karyotype=$1
Iterations=$2
Fusion="$(($Karyotype+$Karyotype-1))"
VERSION='1.0';
RSCRIPT=;
QBATCH=;


# Error function
# ============================================================

function error () {
    printf "[%s] ERROR: " `basename $0` >&2;
    echo "$2." >&2;
    exit $1;
}


# Help usage function
# ============================================================

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ -z "$1" ]; then
    echo "Usage: `basename $0` <1n karyotype>"
    echo "This script models all fission and fusion possibilities from an initial karyotype to one million iterations in order to count the number of fusions that perfectly reverse a prior fission. This script submits jobs to a SLURM or SGE cluster."
    exit 0
fi


## Check for file errors and that external tools are accessible
# ============================================================

if [[ -z "${RSCRIPT}" || ! -x "${RSCRIPT}" ]]; then
    RSCRIPT=`which Rscript`;
    if [[ -z "${RSCRIPT}" || ! -x "${RSCRIPT}" ]]; then
        error 127 "Rscript not found or not executable";
    fi
fi

if [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
    QBATCH=`which qbatch`;
    if [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
        error 127 "qbatch not found or not executable";
    fi
fi

###############################################################################
# Run:
###############################################################################

for i in $(seq 1 $Karyotype); do for u in $(seq 1 $Fusion); do echo "$Prefix/bin/fission_fusion_model.py -o ${i}_${u} -p 16 -i 1000000 $Karyotype ${i} ${u}"; done; done >script.sh
$QBATCH submit -W -t 1:0:0 -p 16 -R 1 script.sh
for i in $(seq 1 $Karyotype); do printf ",${i}"; done >Kar_$Karyotype.csv && echo "" >>Kar_$Karyotype.csv
for u in $(seq $Fusion -1 1); do printf "${u}" >>Kar_$Karyotype.csv && for i in $(seq 1 $Karyotype); do head -2 ${i}_${u}.summary | tail -n1 | awk '{if ($2>=0) {printf ","$2/1000000} else {printf ","$2/10}}'; done >>Kar_$Karyotype.csv && echo "" >>Kar_$Karyotype.csv; done
$RSCRIPT $Prefix/bin/heatmap.R Kar_$Karyotype.csv
