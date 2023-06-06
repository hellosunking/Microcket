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
	echo "Usage: $0 <wk.dir> [genome=hg38] [enzyme=HindIII] [thread=16]" > /dev/stderr
	exit 2
fi

wkDIR=$1
genome=${2:-hg38}
enzyme=${3:-HindIII}
THREAD=${4:-16}

JUICER_HOME=/lustre/home/ksun/software/juicer-1.6
JUICER_EXE=$JUICER_HOME/scripts
## the FASTQ files should be under $wkDIR/fastq
## check the enzyme used during lib-prep, commonly ones: HindIII, Mbol, DpnII
## generate the restriction sites
## cd $JUICER/restriction_sites; python $JUICER/misc/generate_site_positions.py $enzyme hg38 $JUICER/references/$genome.fasta

cd $wkDIR
touch juicer.start
sh $JUICER_EXE/juicer.sh -t $THREAD -D $JUICER_HOME \
	-g $genome -z $JUICER_EXE/references/$genome.fasta \
	-p $JUICER_EXE/restriction_sites/$genome.chrom.sizes \
	-s $enzyme -y $JUICER_EXE/restriction_sites/${genome}_$enzyme.txt \
	-d $PWD
touch juicer.end

## set "-s none" and do not specify "-y" for Micro-C data

## generate filtered pairs
#perl filter.Juicer.pl -p $sid.filtered.pairs aligned/merged_nodups.txt

## files to check running time (within the $wkDIR directory):
# 1. Read alignment: $sid.juicer.start to splits/XXX.fastq.gz.sam (use the newest if there is more than 1 file)
# 2. interaction extraction: splits/XXX.fastq.gz.sam to aligned/merged_nodups.txt
