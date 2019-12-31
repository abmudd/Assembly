#!/bin/bash

###############################################################################
# Setup
###############################################################################

# Set hard and optional variables
# ============================================================

RUNTIME='12:0:0';
SJSUPPORT='20';
VERSION='1.0';
DATE=`date +%F`;
COMMAND="$0 $*";
CTHREADS="32";
LTHREADS="0";
INDEXMEM="240";
ALIGNMEM="110";


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
    printf "Usage: %s [--workdir STR] [--RNA-seq-dir STR] [--ref-genome STR]\n" `basename $0` >&2;
    printf "       [--cthreads INT] [--index-mem INT] [--lthreads INT] [--module STR] \n" >&2;
    printf "       [--runtime HH:MM:SS] [--sj-support INT] [--help] [-h]\n" >&2;
    printf "\n" >&2;
    printf "Bash pipeline to align RNA-seq reads to a genome with splice junctions using STAR.\n" >&2;
    printf "\n" >&2;
    printf "Required arguments:\n" >&2;
    printf "       -w, --workdir STR       path for working directory\n" >&2;
    printf "       -d, --RNA-seq-dir STR   path for directory with RNA-seq files\n" >&2;
    printf "       -r, --ref-genome STR    path for genome fasta file\n" >&2;
    printf "\n" >&2;
    printf "Optional arguments:\n" >&2;
    printf "       -a, --align-mem INT     memory for STAR align step in GB [110]\n" >&2;
    printf "       -c, --cthreads INT      run pipeline on cluster with INT number of processors [32]\n" >&2;
    printf "       -i, --index-mem INT     memory for STAR index step in GB [240]\n" >&2;
    printf "       -l, --lthreads INT      run pipeline locally with INT number of processors [0]\n" >&2;
    printf "       -m, --module STR        module load STR (recommend star/2.5.3a and samtools/1.6)\n" >&2;
    printf "       -t, --runtime HH:MM:SS  runtime for job submissions [12:0:0]\n" >&2;
    printf "       -s, --sj-support INT    minimum support for splice junctions [20]\n" >&2;
    printf "       -h, --help              show this help message and exit\n" >&2;
    printf "\n" >&2;
    printf "Notes:\n" >&2;
    printf "   1. This script assumes a standard UNIX/Linux install and access to a SGE or SLURM\n" >&2;
    printf "      job-scheduling system via qbatch.\n" >&2;
    printf "   2. Running locally/lthreads takes priority over running on the cluster/cthreads.\n" >&2;
    printf "      If running locally, caution should be taken about the lthreads flag as the STAR\n" >&2;
    printf "      align step is memory intensive. Too many lthreads will cause a segmentation fault.\n" >&2;
    printf "   3. This script assumes that fastq files end in the format [0,1,2][fastq format], are\n" >&2;
    printf "      the only fastq files in the RNA-seq-dir, and are all only one fastq format.\n" >&2;
    printf "      Accepted fastq formats are .fastq, .fq, .fastq.gz, and .fq.gz. These files should\n" >&2;
    printf "      already be adapter trimmed.\n" >&2;
    printf "   4. It is highly recommended that absolute paths be passed via the required flags, as\n" >&2;
    printf "      those paths may be used in submission to the cluster.\n" >&2;
    printf "\n" >&2;
}


## Retrieve options on the command line and check for errors
# ============================================================

HELP_MESSAGE=;

while [[ -n $@ ]]; do
    case "$1" in
        '--workdir') shift; WORKDIR=$1;;
        '--RNA-seq-dir') shift; RNADIR=$1;;
        '--ref-genome') shift; GENOME=$1;;
        '--align-mem') shift; ALIGNMEM=$1;;
        '--cthreads') shift; CTHREADS=$1;;
        '--index-mem') shift; INDEXMEM=$1;;
        '--lthreads') shift; LTHREADS=$1;;
        '--module') shift; module load $1;;
        '--runtime') shift; RUNTIME=$1;;
        '--sj-support') shift; SJSUPPORT=$1;;
        '--help') HELP_MESSAGE=1;;
        '-w') shift; WORKDIR=$1;;
        '-d') shift; RNADIR=$1;;
        '-r') shift; GENOME=$1;;
        '-a') shift; ALIGNMEM=$1;;
        '-c') shift; CTHREADS=$1;;
        '-i') shift; INDEXMEM=$1;;
        '-l') shift; LTHREADS=$1;;
        '-m') shift; module load $1;;
        '-t') shift; RUNTIME=$1;;
        '-s') shift; SJSUPPORT=$1;;
        '-h') HELP_MESSAGE=1;;
        -*) usage; error 2 "Invalid option: ${1}";;
        *) break;;
    esac;
    shift;
done

if [[ -n "${HELP_MESSAGE}" ]]; then
    usage;
    exit 1;

elif [[ -z "${WORKDIR}" || ! -d "${WORKDIR}" ]]; then
    usage; error 1 "WORKDIR not defined or does not exist";

elif [[ -z "${RNADIR}" || ! -d "${RNADIR}" ]]; then
    usage; error 1 "RNADIR not defined or does not exist";

elif [[ -z "${GENOME}" || ! -f "${GENOME}" ]]; then
    usage; error 1 "GENOME not defined or does not exist";

else
    printf "[%s] Starting %s %s %s %s %s %s\n" `basename $0` `date`     >&2;
    printf "[%s] Command-line: $COMMAND\n" `basename $0` >&2;
    printf "[%s] Version: $VERSION\n" `basename $0` >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "WORKDIR"     $WORKDIR   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "RNADIR"      $RNADIR    >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "GENOME"      $GENOME    >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "RUNTIME"     $RUNTIME   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "SJSUPPORT"   $SJSUPPORT >&2;
fi


# Check that external tools are accessible
# ============================================================

CAT=`which cat`;
ZCAT=`which zcat`;
PERL=`which perl`;
QBATCH=`which qbatch`;
STAR=`which STAR`;
SAMTOOLS=`which samtools`;

if [[ -z "${PERL}" || ! -x "${PERL}" ]]; then
    error 127 "perl not in PATH env variable or not executable";

elif [[ -z "${CAT}" || ! -x "${CAT}" ]]; then
    error 127 "cat not in PATH env variable or not executable";

elif [[ -z "${ZCAT}" || ! -x "${ZCAT}" ]]; then
    error 127 "zcat not in PATH env variable or not executable";

elif [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
    error 127 "qbatch not in PATH env variable or not executable";

elif [[ -z "${STAR}" || ! -x "${STAR}" ]]; then
    error 127 "star not in PATH env variable or not executable";

elif [[ -z "${SAMTOOLS}" || ! -x "${SAMTOOLS}" ]]; then
    error 127 "samtools not in PATH env variable or not executable";

fi


# Determine fastq format
# ============================================================

GZ=;
UNGZ=;
CATPATH=;

if [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name '*.fastq.gz' -size +1c | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fastq.gz";
    CATPATH=$ZCAT;
    GZ=1;

elif [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name '*.fq.gz' -size +1c | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fq.gz";
    CATPATH=$ZCAT;
    GZ=1;
fi

if [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name '*.fastq' -size +1c | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fastq";
    CATPATH=$CAT;
    UNGZ=1;

elif [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name '*.fq' -size +1c | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fq";
    CATPATH=$CAT;
    UNGZ=1;
fi

if [[ -n $GZ && -n $UNGZ ]]; then
    error 1 "Mixed compressed and uncompressed fastq files detected. Please standardize your files to either all compressed or uncompressed"

elif [[ -z $GZ && -z $UNGZ ]]; then
    usage; error 1 "Fastq files could not be found in RNADIR";
fi


# Determine if fastq files are single, paired, or both
# ============================================================

if [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name "*0${FQFORMAT}" -size +1c | wc -l)" -ne 0 ]]; then
    READEND="single";
fi

if [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name "*1${FQFORMAT}" -size +1c | wc -l)" -ne 0 && "$(find ${RNADIR}/* -maxdepth 0 -type f -name "*2${FQFORMAT}" -size +1c | wc -l)" -ne 0 ]]; then
    if [[ $READEND == "single" ]]; then
	READEND="both";
    else
	READEND="paired";
    fi
fi

if [[ $READEND == "paired" || $READEND == "both" ]]; then
    if [[ "$(find ${RNADIR}/* -maxdepth 0 -type f -name "*1${FQFORMAT}" -size +1c | wc -l)" != "$(find ${RNADIR}/* -maxdepth 0 -type f -name "*2${FQFORMAT}" -size +1c | wc -l)" ]]; then
	error 1 "Unequal number of paired files detected";
    fi
fi


# Determine max read length - 1
# ============================================================

READLEN=`for f in $RNADIR/*[0-2]$FQFORMAT; do $CATPATH ${f} | head -n4000; done | awk '{if (NR % 4 == 2) { print length($1)-1}}' | sort -k1,1nr | head -1`;


# Determine threads
# ============================================================

if [[ $LTHREADS -ne 0 ]]; then
    STAR_THREADS=$(printf "%.0f\n" "$((($(grep -c processor /proc/cpuinfo) - 1) / $LTHREADS))");
    CTHREADS="0";
else
    STAR_THREADS=$CTHREADS;
fi

if [[ -z $CTHREADS ]]; then
    error 1 "Unable to determine CTHREADS";

elif [[ -z $LTHREADS ]]; then
    error 1 "Unable to determine LTHREADS";

elif [[ -z $READLEN ]]; then
    error 1 "Unable to determine READLEN";

elif [[ -z $READEND ]]; then
    error 1 "Unable to determine READEND";

elif [[ -z $FQFORMAT ]]; then
    error 1 "Unable to determine FQFORMAT";

else
    printf "[%s] PARAM: %s = %s\n" `basename $0` "CTHREADS"    $CTHREADS  >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "LTHREADS"    $LTHREADS  >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "STARTHREADS" $STAR_THREADS  >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "ALIGNMEM"    $ALIGNMEM  >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "INDEXMEM"    $INDEXMEM  >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "READLEN-1"   $READLEN   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "READEND"     $READEND   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "FQFORMAT"    $FQFORMAT  >&2;
fi


# Make directories
# ============================================================

mkdir -p $WORKDIR/Index/1st_Index;
mkdir -p $WORKDIR/Index/2nd_Index;
mkdir -p $WORKDIR/Split_Fastq;
mkdir -p $WORKDIR/RNA_Mapping/1st_Mapping;
mkdir -p $WORKDIR/RNA_Mapping/2nd_Mapping;
mkdir -p $WORKDIR/Final_Bam;


###############################################################################
# First Round:
###############################################################################

# Split transcriptome files
# ============================================================

printf "[%s] Splitting transcriptome files \n" `basename $0` >&2;

if [[ -f "${WORKDIR}/Split_Fastq.submit" && -f "${WORKDIR}/Split_Fastq.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    for f in $RNADIR/*$FQFORMAT; do 
	sample=`basename $f $FQFORMAT`; 
	echo "$CATPATH $f | split -dl 4000000 -a 4 - $WORKDIR/Split_Fastq/${sample}_";
    done >$WORKDIR/Split_Fastq.submit;
    
    $QBATCH submit -T $STAR_THREADS -W -n Split_Fastq.submit $WORKDIR/Split_Fastq.submit $WORKDIR/BATCH_Split_Fastq;

    if [[ "$(grep 'done' ${WORKDIR}/BATCH_Split_Fastq/Split_Fastq.submit.e* | wc -l)" == "$(ls ${WORKDIR}/BATCH_Split_Fastq/Split_Fastq.submit.e* | wc -l)" ]]; then
	touch $WORKDIR/Split_Fastq.done;
    else
	error 1 "Splitting transcriptome files failed; please identify error and restart";
    fi
fi


# Index genome
# ============================================================

printf "[%s] FIRST ROUND \n" `basename $0` >&2;
printf "[%s] Indexing genome \n" `basename $0` >&2;

if [[ -f "${WORKDIR}/1st_Index.submit" && -f "${WORKDIR}/1st_Index.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    echo "cd $WORKDIR && $STAR --runMode genomeGenerate --runThreadN 1 --genomeDir $WORKDIR/Index/1st_Index --genomeFastaFiles $GENOME --limitGenomeGenerateRAM=$(($INDEXMEM * 1000000000))" >$WORKDIR/1st_Index.submit;
    $QBATCH submit -T $LTHREADS -W -R $INDEXMEM -t $RUNTIME -n Index_1st.submit $WORKDIR/1st_Index.submit $WORKDIR/BATCH_1st_Index;

    if [[ "$(grep 'done' ${WORKDIR}/BATCH_1st_Index/1st_Index.submit.e* | wc -l)" == "$(wc -l <${WORKDIR}/1st_Index.submit)" ]]; then
	touch $WORKDIR/1st_Index.done;
    else
	error 1 "Indexing genome failed; please identify error and restart";
    fi
fi


# Align split transcriptome files
# ============================================================

printf "[%s] Aligning split transcriptome files \n" `basename $0` >&2;

if [[ -f "$WORKDIR/1st_Align.submit" && -f "$WORKDIR/1st_Align.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    printf '' >$WORKDIR/1st_Align.submit;

    # For paired reads
    # ============================================================

    if [[ $READEND == "paired" || $READEND == "both" ]]; then
	ls $WORKDIR/Split_Fastq/*1_[0-9][0-9][0-9][0-9] | rev | cut -c7- | cut -d'/' -f1 | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
	while read sample N; do 
	    for n in `seq 0 $N`; do 
		n=`printf %04d $n`; 
		echo "$STAR --genomeDir $WORKDIR/Index/1st_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}1_${n} $WORKDIR/Split_Fastq/${sample}2_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/1st_Mapping/${sample}_${n} --runThreadN $STAR_THREADS --chimOutType SeparateSAMold --chimSegmentMin 20 --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >$WORKDIR/1st_Align.submit;
    fi

    # For single reads
    # ============================================================

    if [[ $READEND == "single" || $READEND == "both" ]]; then
	ls $WORKDIR/Split_Fastq/*0_[0-9][0-9][0-9][0-9] | rev | cut -c7- | cut -d'/' -f1 | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
        while read sample N; do 
	    for n in `seq 0 $N`; 
	    do n=`printf %04d $n`;
		echo "$STAR --genomeDir $WORKDIR/Index/1st_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}0_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/1st_Mapping/${sample}_${n} --runThreadN $STAR_THREADS --chimOutType SeparateSAMold --chimSegmentMin 20 --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >>$WORKDIR/1st_Align.submit;
    fi

    $QBATCH submit -T $LTHREADS -W -t $RUNTIME -p $CTHREADS -R $ALIGNMEM -n Align_1.submit $WORKDIR/1st_Align.submit $WORKDIR/BATCH_1st_Align;

    if [[ "$(grep -o 'done' ${WORKDIR}/BATCH_1st_Align/1st_Align.submit.e* | wc -l)" == "$(wc -l <${WORKDIR}/1st_Align.submit)" ]]; then
	touch $WORKDIR/1st_Align.done;
    else
	error 1 "Aligning split transcriptome files failed; please identify error and restart";
    fi
fi


# Identify splice junctions
# ============================================================

printf "[%s] Identifying splice junctions \n" `basename $0` >&2;

if [[ -f "${WORKDIR}/all_sj.hist" && -f "${WORKDIR}/Find_SJ.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    cat $WORKDIR/RNA_Mapping/1st_Mapping/*SJ.out.tab | awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if ($7>0) print $1,$2,$3,strChar[$4],$5,$7}' | $PERL -ane '$d{"$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]"} += $F[5]; END{map{print("$_\t$d{$_}\n")}keys(%d)}' | sort -k1,1 -k2,2n >$WORKDIR/all_sj.out;
    cut -f6 $WORKDIR/all_sj.out | sort -n | uniq -c | awk '{print $2 "\t" $1}' >$WORKDIR/all_sj.hist
    if [[ "$(find ${WORKDIR}/* -maxdepth 0 -type f -name 'all_sj.hist' -size +1c | wc -l)" -ne 0 ]]; then
	touch $WORKDIR/Find_SJ.done;
    else
	error 1 "Identifying splice junctions failed; please identify error and restart";
    fi
fi


###############################################################################
# Second Round:
###############################################################################

# Filter splice junctions
# ============================================================

printf "[%s] Filtering splice junctions \n" `basename $0` >&2;

if [[ -f "${WORKDIR}/sjdb.out" && -f "${WORKDIR}/Filter_SJ.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    awk -v sjs="$SJSUPPORT" '{if($6>sjs) print $0}' $WORKDIR/all_sj.out >$WORKDIR/sjdb.out;
    if [[ "$(find ${WORKDIR}/* -maxdepth 0 -type f -name 'sjdb.out' -size +1c | wc -l)" -ne 0 ]]; then
	touch $WORKDIR/Filter_SJ.done;
    else
	error 1 "Filtering splice junctions failed; please identify error and restart";
    fi
fi


# Index genome
# ============================================================

printf "[%s] SECOND ROUND \n" `basename $0` >&2;
printf "[%s] Indexing genome \n" `basename $0` >&2;

if [[ -f "${WORKDIR}/2nd_Index.submit" && -f "${WORKDIR}/2nd_Index.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    echo "cd $WORKDIR && $STAR --runMode genomeGenerate --runThreadN 1 --genomeDir $WORKDIR/Index/2nd_Index --genomeFastaFiles $GENOME --sjdbFileChrStartEnd $WORKDIR/sjdb.out --limitGenomeGenerateRAM=$(($INDEXMEM * 1000000000)) --sjdbOverhang $READLEN" >$WORKDIR/2nd_Index.submit;

    $QBATCH submit -T $LTHREADS -W -R $INDEXMEM -t $RUNTIME -n Index_2nd.submit $WORKDIR/2nd_Index.submit $WORKDIR/BATCH_2nd_Index;

    if [[ "$(grep 'done' ${WORKDIR}/BATCH_2nd_Index/2nd_Index.submit.e* | wc -l)" != "$(wc -l <$WORKDIR/2nd_Index.submit)" ]]; then
	error 1 "Indexing genome failed; please identify error and restart";
    else
	touch $WORKDIR/2nd_Index.done;
    fi
fi


# Align split transcriptome files
# ============================================================

printf "[%s] Aligning split transcriptome files \n" `basename $0` >&2;

if [[ -f "$WORKDIR/2nd_Align.submit" && -f "$WORKDIR/2nd_Align.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    printf '' >$WORKDIR/2nd_Align.submit;

    # For paired reads
    # ============================================================

    if [[ $READEND == "paired" || $READEND == "both" ]]; then
	ls $WORKDIR/Split_Fastq/*1_[0-9][0-9][0-9][0-9] | rev | cut -c7- | cut -d'/' -f1 | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
        while read sample N; do 
	    for n in `seq 0 $N`; do
		n=`printf %04d $n`; 
		echo "$STAR --genomeDir $WORKDIR/Index/2nd_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}1_${n} $WORKDIR/Split_Fastq/${sample}2_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/2nd_Mapping/${sample}_${n} --runThreadN $STAR_THREADS --chimOutType SeparateSAMold --chimSegmentMin 20 --sjdbFileChrStartEnd $WORKDIR/sjdb.out --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >$WORKDIR/2nd_Align.submit;
    fi

    # For single reads
    # ============================================================

    if [[ $READEND == "single" || $READEND == "both" ]]; then
	ls $WORKDIR/Split_Fastq/*0_[0-9][0-9][0-9][0-9] | rev | cut -c7- | cut -d'/' -f1 | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
	while read sample N; do 
	    for n in `seq 0 $N`; do 
		n=`printf %04d $n`; 
		echo "$STAR --genomeDir $WORKDIR/Index/2nd_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}0_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/2nd_Mapping/${sample}_${n} --runThreadN $STAR_THREADS --chimOutType SeparateSAMold --chimSegmentMin 20 --sjdbFileChrStartEnd $WORKDIR/sjdb.out --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >>$WORKDIR/2nd_Align.submit;
    fi

    $QBATCH submit -T $LTHREADS -W -t $RUNTIME -p $CTHREADS -R $ALIGNMEM -n Align_2.submit $WORKDIR/2nd_Align.submit $WORKDIR/BATCH_2nd_Align;

    if [[ "$(grep -o 'done' ${WORKDIR}/BATCH_2nd_Align/2nd_Align.submit.e* | wc -l)" == "$(wc -l <${WORKDIR}/2nd_Align.submit)" ]]; then
	touch $WORKDIR/2nd_Align.done;
    else
	error 1 "Aligning split transcriptome files failed; please identify error and restart";
    fi
fi


###############################################################################
# Finish:
###############################################################################

# Merge second round bam files
# ============================================================  

printf "[%s] Merging final bam files \n" `basename $0` >&2;

if [[ -f "$WORKDIR/Final_Bam.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    ls $WORKDIR/Split_Fastq/*[0-2]_[0-9][0-9][0-9][0-9] | rev | cut -c8- | cut -d'/' -f1 | rev | sort | uniq | \
    while read sample; do 
	$SAMTOOLS merge -f -@ $STAR_THREADS $WORKDIR/Final_Bam/${sample}.bam $WORKDIR/RNA_Mapping/2nd_Mapping/${sample}*Aligned.sortedByCoord.out.bam; 
    done;
fi

if [[ "$(find ${WORKDIR}/Final_Bam/* -maxdepth 0 -type f -name '*bam' -size +1c | wc -l)" == "$(ls ${RNADIR}/*[0-1]${FQFORMAT} | wc -l)" ]]; then
    touch $WORKDIR/Final_Bam.done;
    printf "[%s] Finished %s %s %s %s %s %s\n" `basename $0` `date` >&2;
else
    error 1 "Merging final bam files failed; please identify error and restart";
fi

exit 0;
