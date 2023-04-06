#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

## this script is used to call loops from hESC Micro-C data: pooled high-depth, Microcket, and Distiller

juicer_tools=/mnt/software/Microcket-1.1/bin/juicer_tools.2.13.07.jar
RES=10000
THREAD=16

## VC_SQRT is a normalization methods existing in all the 3 files
java -jar $juicer_tools hiccups --cpu --threads $THREAD -r $RES -k VC_SQRT \
	4DNFI2TK7L2F.hic hiccups.highDepth

java -jar $juicer_tools hiccups --cpu --threads $THREAD -r $RES -k VC_SQRT \
	Microcket.hic hiccups.Microcket

java -jar $juicer_tools hiccups --cpu --threads $THREAD -r $RES -k VC_SQRT \
	Distiller.hic hiccups.Distiller

for T in highDepth Microcket Distiller
do
	cat hiccups.$T/postprocessed_pixels_*.bedpe | \
	perl -lane 'next if /^#/ || $F[16]>=0.01; $F[0]=~s/^chr//; $F[3]=~s/^chr//; print "$F[0]:$F[1]-$F[3]:$F[4]"' >hiccups.$T/loop
done

R --slave < plot.loop.overlap.R

