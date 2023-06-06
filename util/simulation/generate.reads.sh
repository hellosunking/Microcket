#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

## simulate Hi-C reads using sim3C
enzyme=MboI
number=10000000
cycle=150
insert=300
## need to manually update the path to genome file
genome=hg38.fa

rm -f sim3C.fastq sim3C.log profile.tsv
sim3C --convert --dist uniform -n $number -l $cycle --insert-mean $insert \
	-e $enzyme -m hic --machine-profile HiSeqXtruSeqL150 $genome sim3C.fastq

perl split.sim3C.pl sim3C.fastq sim3C

