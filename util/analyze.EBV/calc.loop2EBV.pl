#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.pairs[.gz]> [bin=1000] [min.loop=5]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $target = 'chrEBV';
my $bin = $ARGV[1] || 1000;
my $min_loop = $ARGV[2] || 5;

my %m;
my %link;
open IN, "less $ARGV[0] |" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##A00153:1085:HGWCWDSX3:2:2146:28872:35712	chr1	87209	chrEBV	165386	-	+
	next unless $l[1]=~/chr\d+$/ && $l[3] eq $target;

	$l[2] = int($l[2]/$bin);
	$m{"$l[1]\t$l[2]"} ++;
}
close IN;

## count table for bedgraph
foreach my $k ( keys %m ) {
	next if $m{$k} < $min_loop;
	my ($chr, $pos) = split /\t/, $k;
	print join("\t", $chr, $pos*$bin, $pos*$bin+$bin, $m{$k}), "\n";
}

