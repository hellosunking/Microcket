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
	echo "Usage: $0 <in.fq.list> <sid> [Enzyme=HindIII]" > /dev/stderr
	exit 2
fi

fqlist=$1
sid=$2
ENZYME=${3:-HindIII}

THREAD=8
genome=/mnt/Genomes/hg38/hg38.fa
BWAREF=/mnt/software/bwa-0.7.17/index/hg38p13

files=""
while read R1 R2 extra
do
	files="$files $R1 $R2"
done < $fqlist

touch fanc.$sid.start
fanc auto --no-hic -g $genome -t $THREAD -i $BWAREF -r $ENZYME -n $sid $files fanc.$sid
touch fanc.$sid.end

## files to check running time (within the fanc.$sid directory):
# 1. Read alignment: fanc.$sid.start to sam/XXX.bam (use the newest if there is more than 1 file)
# 2. interaction extraction: sam/XXX.bam to pairs/*.pairs

