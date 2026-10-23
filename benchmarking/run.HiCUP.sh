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
	echo "Usage: $0 <in.conf>"
	exit 2
fi > /dev/stderr

conf=$1

HICUP=~/software/HiCUP-0.8.3
hicup2pairs=/lustre/home/ksun/kqsoft/Microcket/Microcket-1.3/benchmarking/hicup2pairs

## The OUTDIR should be identical to that in the configuration file
OUTDIR=`cat $conf | perl -ne 'next if /^\s*#/; print $1 if /Outdir:\s*(\S+)/i'`
if [ "$OUTDIR" == "" ]
then
	echo "ERROR: could not infer output directory from configuration file!"
	exit 1
else
	mkdir -p $OUTDIR
	echo "INFO: use output directory $OUTDIR"
fi

touch $OUTDIR/HiCUP.start
$HICUP/hicup --config $conf

cd $OUTDIR
touch HiCUP.aligned
## there could be more than 1 sam files
#cat *sam | $HICUP/Conversion/hicup2juicer -
$hicup2pairs *sam >HiCUP.pairs
touch HiCUP.end

## files to check running time:
# 1. Preprocessing: HiCUP.start to hicup_truncater_summary_XXX.txt
# 2. Alignment: hicup_truncater_summary_XXX.txt to HiCUP.aligned
# 3. Interaction extraction: HiCUP.aligned to $sid.HiCUP.end
# Note that HiCUP deletes the intermediate files and reports the alignments after filtering and de-duplication
