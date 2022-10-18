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
	echo "Usage: $0 <genome.fa> <genome.id> [thread=all] [target=STAR|both]"
	echo "Note: gzip-file is NOT acceptable."
	exit 2
fi >/dev/stderr

ver=1.1.0

PRG=`dirname $0`
index=$PRG/../index
BIN=$PRG/../bin
anno=$PRG/../anno

fasta=$1
gid=$2
thread=${3:-0}
target=${4:-STAR}
BWA=0

if [ $thread == 0 ]
then
	thread=`cat /proc/cpuinfo | grep processor | wc -l`
fi

if [ $target == "both" ] || [ $target == "Both" ] || [ $target == "BOTH" ]
then
	BWA=1
elif [ $target == "star" ] || [ $target == "Star" ] || [ $target == "STAR" ]
then
	BWA=0
else
	echo "ERROR: Unknown parameter! Must be both or STAR!" >/dev/stderr
	exit 1
fi

## STAR
echo "======== Build index for STAR ========"
mkdir -p $index/STAR/$gid
$BIN/STAR --runThreadN $thread --runMode genomeGenerate --genomeDir $index/STAR/$gid --genomeFastaFiles $fasta

## bwa
if [ $BWA == 1 ]
then
	echo "======== Build index for BWA  ========"
	$BIN/bwa index -p $index/BWA/$gid $fasta
fi

## generate sam header
echo "======== Generating annotation files ========"
cp $index/STAR/$gid/chrNameLength.txt $anno/$gid.anno
echo -e "@HD\tVN:1.0\tSO:coordinate" > $anno/$gid.sam.header
cat $index/STAR/$gid/chrNameLength.txt | perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]"' >> $anno/$gid.sam.header
echo -e "@PG\tID:Microcket\tPN:Microcket\tVN:$ver\tDS:$gid" >> $anno/$gid.sam.header

echo "======== Job done: now you can use '-g $gid' in microcket ========"

