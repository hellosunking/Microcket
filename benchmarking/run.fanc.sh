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
	echo "Usage: $0 <in.fq.list> <sid> [aligner=bwa|bowtie2] [Enzyme=HindIII]" > /dev/stderr
	exit 2
fi

export PATH=$PATH:/lustre/home/ksun/lib/Python.3.10.7.build/bin/

fqlist=$1
sid=$2
aligner=${3:-bwa}
ENZYME=${4:-HindIII}

THREAD=16

fanc=/lustre/home/ksun/lib/Python.3.10.7.build/bin/fanc
genome=/lustre/home/ksun/Genomes/hg38/hg38.p13.fa
if [ $aligner == "bowtie2" ]
then
	Index=/lustre/home/ksun/Genomes/bowtie2.index/hg38
else
	Index=/lustre/home/ksun/software/bwa-0.7.17/index/hg38p13
fi

files=""
while read R1 R2 extra
do
	files="$files $R1 $R2"
done < $fqlist

#echo $files
touch fanc.$sid.start

$fanc auto --no-hic -g $genome -t $THREAD -i $Index \
	-r $ENZYME -n $sid $files fanc.$sid

touch fanc.$sid.end

## files to check running time (within the fanc.$sid directory):
# 1. Read alignment: fanc.$sid.start to sam/XXX.bam (use the newest if there is more than 1 file)
# 2. interaction extraction: sam/XXX.bam to pairs/*.pairs

