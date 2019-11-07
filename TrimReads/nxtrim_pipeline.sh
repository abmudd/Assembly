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

NXTRIM=;
reheader=`dirname $0`/bin/reheader.py;
single2pair=`dirname $0`/bin/single2pair.pl;
VERSION='1.0';

printf "[%s] " `basename $0` >&2;
echo $@ >&2;


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
    echo "Usage: `basename $0` <prefix path for MP fastq.gz files>"
    echo "This script assumes that the files are named prefixpath_1.fastq.gz and prefixpath_2.fastq.gz"
    exit 0
fi


## Check for file errors and that external tools are accessible
# ============================================================

if [[ ! -f "$1_1.fastq.gz" || ! -f "$1_2.fastq.gz" ]]; then
    error 1 "Could not find $1_1.fastq.gz or $1_2.fastq.gz."
elif [[ -z "${reheader}" || ! -f "${reheader}" ]]; then
    error 1 "reheader.py not defined or does not exist";
elif [[ -z "${single2pair}" || ! -f "${single2pair}" ]]; then
    error 1 "single2pair.py not defined or does not exist";
fi

if [[ -z "${NXTRIM}" || ! -x "${NXTRIM}" ]]; then
    NXTRIM=`which nxtrim`;
    if [[ -z "${NXTRIM}" || ! -x "${NXTRIM}" ]]; then
        error 127 "NXTRIM not found or not executable";
    fi
fi


###############################################################################
# Run:
###############################################################################

$NXTRIM --separate --similarity 0.9 --joinreads -1 $1_1.fastq.gz -2 $1_2.fastq.gz -O $1

$single2pair $1.se.fastq.gz $1.se &
$reheader -g -o $1_R1.mp.rh.fastq.gz $1_R1.mp.fastq.gz &
$reheader -g -o $1_R2.mp.rh.fastq.gz $1_R2.mp.fastq.gz &
$reheader -g -o $1_R1.unknown.rh.fastq.gz $1_R1.unknown.fastq.gz &
$reheader -g -o $1_R2.unknown.rh.fastq.gz $1_R2.unknown.fastq.gz &
$reheader -g -o $1_R1.pe.rh.fastq.gz $1_R1.pe.fastq.gz &
$reheader -g -o $1_R2.pe.rh.fastq.gz $1_R2.pe.fastq.gz &
wait

cat $1_R1.mp.rh.fastq.gz $1_R1.unknown.rh.fastq.gz >$1.mp.final_1.fastq.gz &
cat $1_R2.mp.rh.fastq.gz $1_R2.unknown.rh.fastq.gz >$1.mp.final_2.fastq.gz &
cat $1_R1.pe.rh.fastq.gz $1.se_1.fastq.gz >$1.pe.final_1.fastq.gz &
cat $1_R2.pe.rh.fastq.gz $1.se_2.fastq.gz >$1.pe.final_2.fastq.gz &
wait

exit 0
