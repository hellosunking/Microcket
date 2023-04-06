#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

## Distiller, use BWA + pairtools
## this sccript is adapted from https://hichip.readthedocs.io/en/latest/fastq_to_bam.html

if [ $# -lt 2 ]
then
	echo "Usage: $0 <in.fq.list> <sid> [thread=16]" > /dev/stderr
	exit 2
fi

fqlist=$1
sid=$2
THREAD=${3:-16}

BWAREF=/mnt/software/bwa-0.7.17/index/hg38p13
genomeinfo=/mnt/Genomes/bowtie2.index/hg38.info

let proc_in=$THREAD/2
let proc_out=$THREAD-$proc_in

touch $sid.Distiller.start

firstFQ=1
while read R1 R2 extra
do
	if [ $firstFQ == 1 ]
	then
		bwa mem -5SP -T0 -v 2 -t $THREAD $BWAREF $R1 $R2
		firstFQ=0
	else
		bwa mem -5SP -T0 -v 2 -t $THREAD $BWAREF $R1 $R2 | grep -v '^@'
	fi
done < $fqlist >$sid.raw.sam

pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $proc_in --nproc-out $proc_out --chroms-path $genomeinfo $sid.raw.sam | \
pairtools sort --tmpdir=/tmp --nproc $THREAD | \
pairtools dedup --nproc-in $proc_in --nproc-out $proc_out --mark-dups --output-stats $sid.stats.txt | \
pairtools split --nproc-in $proc_in --nproc-out $proc_out --output-pairs $sid.final.pairs --output-sam $sid.final.sam

touch $sid.Distiller.finished

## files to check running time (within the HiC-Pro.$sid/ directory):
# 1. Read alignment: $sid.Distiller.start to $sid.raw.sam
# 2. interaction extraction: $sid.raw.sam to $sid.final.pairs
