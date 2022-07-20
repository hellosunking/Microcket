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
	echo "Usage: $0 <genome.fa> <genome.id> [thread=all]" > /dev/stderr
	exit 2
fi

ver=1.0.0

PRG=`dirname $0`
index=$PRG/../index
BIN=$PRG/../bin
anno=$PRG/../anno

fasta=$1
gid=$2
thread=${3:-0}

if [ $thread == 0 ]
then
	thread=`cat /proc/cpuinfo | grep processor | wc -l`
fi

## STAR
echo "======== Build index for STAR ========"
mkdir -p $index/STAR/$gid
$BIN/STAR --runThreadN $thread --runMode genomeGenerate --genomeDir $index/STAR/$gid --genomeFastaFiles $fasta

## bwa
echo "======== Build index for BWA  ========"
$BIN/bwa index -p $index/BWA/$gid $fasta

## generate sam header
echo "======== Generating annotation files ========"
cp $index/STAR/$gid/chrNameLength.txt $anno/$gid.anno
echo -e "@HD\tVN:1.0\tSO:coordinate" > $anno/$gid.sam.header
cat $index/STAR/$gid/chrNameLength.txt | perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]"' >> $anno/$gid.sam.header
echo -e "@PG\tID:Microcket\tPN:Microcket\tVN:$ver\tDS:$gid" >> $anno/$gid.sam.header

echo "======== Job done: now you can use '-g $gid' in microcket ========"

