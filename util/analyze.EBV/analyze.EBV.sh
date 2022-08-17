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
	echo "Usage: $0 <in.pairs[.gz]> <sid>" > /dev/stderr
	echo "The in.pairs is Microcket's output, where you should include EBV genome when building index."
	exit 2
fi

PRG=`dirname $0`
sid=$2

less $1 | grep -w chrEBV >$sid.EBV.pairs
perl $PRG/calc.inter.EBV.matrix.and.circos.pl $sid.EBV.pairs >$sid.EBV2hs.bedgraph 2>$sid.EBV2hs.links
perl $PRG/calc.loop2EBV.pl $sid.EBV.pairs | sort -k1,1 -k2,2n >$sid.hs2EBV.bedgraph

cp $PRG/../anno/4DN.DCIC.header $sid.chrEBV.intra.pairs 
cat $sid.EBV.pairs | perl -lane 'print if $F[0]!~/^#/ && $F[1] eq "chrEBV"' >>$sid.chrEBV.intra.pairs
java -jar $PRG/../bin/juicer_tools.jar pre --threads 4 -r 1000,500,200 $sid.chrEBV.intra.pairs $sid.EBV.hic $PRG/EBV.info

