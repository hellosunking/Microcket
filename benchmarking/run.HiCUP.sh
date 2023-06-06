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
	echo "Usage: $0 <in.conf> <out.dir>"
	echo "The 'out.dir' should be identical to that in the configuration file."
	exit 2
fi > /dev/stderr

conf=$1
OUTDIR=$2

currSHELL=`readlink -f $0`
PRG=`dirname $currSHELL`
hicup2pairs=$PRG/hicup2pairs
HICUP=/mnt/software/HiCUP-0.8.3

mkdir -p $OUTDIR
touch $OUTDIR/HiCUP.start
$HICUP/hicup --config $conf

cd $OUTDIR
touch HiCUP.aligned
## convert to pairs format
perl $hicup2pairs *sam >HiCUP.pairs
touch HiCUP.end

## files to check running time:
# 1. Preprocessing: HiCUP.start to hicup_truncater_summary_XXX.txt
# 2. Alignment: hicup_truncater_summary_XXX.txt to HiCUP.aligned
# 3. Interaction extraction: HiCUP.aligned to $sid.HiCUP.end
# Note that HiCUP deletes the intermediate files and reports the alignments after filtering and de-duplication
