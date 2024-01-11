#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <Distiller.stats.txt>\n\n";
	exit 2;
}

my ($total, $unmapped, $filtered, $dup, $reported) = ( 0,0,0,0,0 );
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	if( $l[0] eq 'total' ) {
		$total = $l[1];
	} elsif ( $l[0] eq 'total_unmapped' ) {
		$unmapped = $l[1];
	} elsif ( $l[0] eq 'total_single_sided_mapped' ) {
		$filtered = $l[1];
	} elsif ( $l[0] eq 'total_dups' ) {
		$dup = $l[1];
	} elsif ( $l[0] eq 'total_nodups' ) {
		$reported = $l[1];
	}
}
close IN;

$total /= 1e6;
$unmapped /= 1e6;
$filtered /= 1e6;
$dup /= 1e6;
$reported /= 1e6;

print   "Total\t$total\n",
		"Low-quality\t0\n",
		"Unmapped\t$unmapped\n",
		"Filtered\t$filtered\n",
		"Duplicate\t$dup\n",
		"Reported\t$reported\n";

