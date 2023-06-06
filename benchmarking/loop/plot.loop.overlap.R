#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for plotting overlaps among loops from different datasets 
#

options( stringsAsFactors=F );
library(VennDiagram)

microcket=read.table("hiccups.Microcket/loop")
distiller=read.table("hiccups.Distiller/loop")
highdepth=read.table("hiccups.highDepth/loop")

venn.diagram( list(A=microcket$V1, B=distiller$V1, C=highdepth$V1),
			 height=2000, width=2000, imagetype="png", 
			 category.names=c("Microcket", "Distiller", "High-depth"),
			 "loop.venn.png" )

