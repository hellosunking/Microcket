#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for plotting overlaps among loops from different datasets 
#

options( stringsAsFactors=F )
library(VennDiagram)

microcket=read.table("hiccups.Microcket/loop")
juicer=read.table("hiccups.Juicer/loop")
hicpro=read.table("hiccups.HiC-Pro/loop")
highdepth=read.table("hiccups.highDepth/loop")

venn.diagram( list(A=microcket$V1, B=highdepth$V1, C=juicer$V1, D=hicpro$V1),
			 height=2000, width=2000, imagetype="png", 
			 category.names=c("Microcket", "High-depth", "Juicer", "HiC-Pro"),
			 "loop.venn.png" )

