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

HICUP=/mnt/software/HiCUP-0.8.3

mkdir -p $OUTDIR
touch $OUTDIR/HiCUP.start
$HICUP/hicup --config $conf

cd $OUTDIR
touch HiCUP.aligned
## there could be more than 1 sam files
cat *sam | $HICUP/Conversion/hicup2juicer -
touch HiCUP.end

## files to check running time:
# 1. Preprocessing: HiCUP.start to hicup_truncater_summary_XXX.txt
# 2. Alignment: hicup_truncater_summary_XXX.txt to HiCUP.aligned
# 3. Interaction extraction: HiCUP.aligned to $sid.HiCUP.end
# Note that HiCUP deletes the intermediate files and reports the alignments after filtering and de-duplication
