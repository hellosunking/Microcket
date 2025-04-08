#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 2 ]
then
	echo "Usage: $0 <genome.fa> <genome.id> [aligner=all|bwa|star] [thread=all]" >/dev/stderr
	echo "INFO: You may consider removing the unplaced contigs using 'util/clean.genome.pl' script." >/dev/stderr
	exit 2
fi

ver=1.4
currSHELL=`readlink -f $0`
PRG=`dirname $currSHELL`
index=$PRG/../index
BIN=$PRG/../bin
anno=$PRG/../anno

fasta=$1
gid=$2
aligner=${3:-all}
thread=${4:-0}

## to lowercase
aligner=${aligner,,}

if [ $thread -eq 0 ]
then
	thread=`cat /proc/cpuinfo | grep processor | wc -l`
	echo "INFO: $thread threads available."
fi

## STAR
doSTAR=no
doBWA=no

if [ "$aligner" == "all" ]
then
	doSTAR=yes
	doBWA=yes
elif [ "$aligner" == "star" ]
then
	doSTAR=yes
	echo "INFO: I will build index for STAR only."
elif [ "$aligner" == "bwa" ]  ## bwa
then
	doBWA=yes
	echo "INFO: I will build index for BWA only."
else
	echo "ERROR: unknown aligner! MUST be all, bwa, or star!"
	exit 1
fi

if [ "$doBWA" == "yes" ]
then
	echo "======== Build index for BWA  ========"
	mkdir -p $index/$gid/BWA/
	$BIN/bwa index -p $index/$gid/BWA/$gid $fasta
fi

if [ "$doSTAR" == "yes" ]
then
	echo "======== Build index for STAR ========"
	mkdir -p $index/$gid/STAR/
	$BIN/STAR --runThreadN $thread --runMode genomeGenerate \
		--genomeDir $index/$gid/STAR \
		--genomeFastaFiles $fasta
fi

## generate sam header
echo "======== Generating annotation files ========"
perl $PRG/make.sam.header.pl $fasta $gid $anno
echo "======== Job done: now you can use '-g $gid' in microcket ========"

