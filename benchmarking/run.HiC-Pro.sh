#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 1 ]
then
	echo
	echo "Usage: $0 <sampleID>"
	echo
	echo "Note: Prepare a sampleID.HiC-Pro.conf file, and put all the reads in fastq/sampleID before running this script."
	echo
	exit 2
fi > /dev/stderr

HICPro=~/software/HiC-Pro-3.0.0/bin/HiC-Pro
sid=$1

touch HiC-Pro.$sid.start
$HICPro -i fastq -o HiC-Pro.$sid -c $sid.HiC-Pro.conf
touch HiC-Pro.$sid.end

## convert HiC-Pro's output to pairs
#cat HiC-Pro.allValidPairs | perl -lane 'print join("\t", $F[0], $F[1], $F[2], $F[4], $F[5], $F[3], $F[6])' >HiC-Pro.pairs

## files to check running time (within the HiC-Pro.$sid/ directory):
# 1. Read alignment: $sid.HiC-Pro.start to bowtie_results/bwt2/$sid/*.bwt2pairs.bam (use the newest if there is more than 1 file)
# 2. interaction extraction: bowtie_results/bwt2/$sid/*.bwt2pairs.bam to hic_results/data/$sid/$sid.allValidPairs
