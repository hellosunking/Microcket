#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.pairs[.gz]> [EBV.bin=200] [hs.bin=1000000] [hotspot.loop.ratio=0.01] >intraEBV.bedgraph 2>interEBV.links\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $target = 'chrEBV';
my $bin = $ARGV[1] || 200;
my $bin_hg38 = $ARGV[2] || 1000000;
my $hotspot_loop_ratio = $ARGV[3] || 0.01;

## param for circos
my $efactor = 5000;	## enlarge factor for EBV genome for circos plots
my $log2 = log(2);

my %m;
my %link;
open IN, "less $ARGV[0] |" or die( "$!" );
my $total = 0;
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##A00153:1085:HGWCWDSX3:2:2146:28872:35712	chr1	87209	chrEBV	165386	-	+
	next unless $l[1]=~/chr\d+$/ && $l[3] eq $target;

	my $i = int($l[4]/$bin);
	$m{$i} ++;

	## for circos
	$l[1]=~s/chr/hs/;
	$l[2] = int($l[2]/$bin_hg38);
	$link{"$i\t$l[1]\t$l[2]"} ++;
	++ $total;
}
close IN;

## count table for bedgraph
my @k = sort {$a<=>$b} keys %m;
my $max = $k[-1];
foreach my $i ( 0..$max ) {
	print join("\t", $target, $i*$bin, $i*$bin+$bin, $m{$i}||0), "\n";
}

## determine the cutoff
## a dynamic cutoff is used in v2
my %loop_counts;
foreach my $k ( sort keys %link ) {
	my $v = $link{$k};
	$loop_counts{$v} ++;
}
my @index  = sort {$b<=>$a} keys %loop_counts;
my $cutoff = int($total * $hotspot_loop_ratio);
#print STDERR "Cutoff\t$cutoff\n";

my $min_loop;
my ($i, $accum) = (0, 0);
foreach $i ( @index ) {
	$accum += $i * $loop_counts{$i};
#	print STDERR "$i\t$accum\n";
	if( $accum > $cutoff ) {
		$min_loop = $i + 1;
		last;
	}
}
$min_loop = $index[0] if $min_loop > $index[0];	## the first one is too large
#print STDERR "Min.loop\t$min_loop\n";

## links for circos
foreach my $k ( sort keys %link ) {
	my $v = $link{$k};
	next if $v < $min_loop;
	$v = int(log($v)/$log2);
	$v .= "p";

	my ($i, $c, $j) = split /\t/, $k;
	$i = $i * $bin * $efactor;
	$j = $j * $bin_hg38;
	print STDERR join("\t", "hs0", $i, $i+$bin*$efactor, $c, $j, $j+$bin_hg38, "thickness=$v"), "\n";
}

