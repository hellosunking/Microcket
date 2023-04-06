#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <sample.mpairstat> <allValidPairs.mergestat>\n";
	print STDERR "There 2 files should be under \"HiC-Pro.output/hic_results/stats/sample/\".\n\n";
	exit 2;
}

my ($total, $unmapped) = ( 0, 0 );
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	if( $l[0] eq 'Total_pairs_processed' ) {
		$total = $l[1];
	} elsif ( $l[0] eq 'Unmapped_pairs' || $l[0] eq 'Pairs_with_singleton' ) {
		$unmapped += $l[1];
	}
}
close IN;
my $mapped = $total-$unmapped;

my ($valid, $reported) = (0,0);
open IN, "$ARGV[1]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	if( $l[0] eq 'valid_interaction' ) {
		$valid = $l[1];
	} elsif( $l[0] eq 'valid_interaction_rmdup' ) {
		$reported = $l[1];
	}
}
close IN;
my $filtered = $mapped - $valid;
my $dup = $valid - $reported;

print   "Total\t$total\n",
		"Low-quality\t0\n",
		"Unmapped\t", $unmapped, "\n",
		"Filtered\t$filtered\n",
		"Duplicate\t$dup\n",
		"Reported\t$reported\n";

