#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <Juicer.work.dir>\n\n";
	exit 2;
}

my $dir = shift;
## NOTE that Juicer (v1.6)'s statistics are incorrect when deal with multiple lanes of data

my $total = 0;
open IN, "cat $dir/splits/*.res.txt |" or die( "$!" );
while( <IN> ) {
	my @l = split /\s+/;
	$total += $l[0];
}
close IN;

my $filtered = 0;
my $reported = 0;
open IN, "$dir/aligned/inter_30.txt" or die( "$!" );
while( <IN> ) {
	s/,//g;
	if( /Intra-fragment.*:\s+(\d+)/ ) {
		$filtered += $1;
	} elsif( /Below MAPQ.*:\s+(\d+)/ ) {
		$filtered += $1;
	} elsif( /Hi-C.*:\s+(\d+)/ ) {
		$reported = $1;
	}
}
close IN;

my $dup = 0;
open IN, "$dir/aligned/dups.txt" or die( "$!" );
while( <IN> ) {
	++ $dup;
}
close IN;

## adjust to Millions
$total /= 1e6;
$filtered /= 1e6;
$reported /= 1e6;
my $non_dup = $filtered+$reported;
$dup /= 1e6;

print	"Total\t$total\n",
		"Low-quality\t0\n",
		"Unmapped\t", $total-$non_dup-$dup, "\n",
		"Filtered\t$filtered\n",
		"Duplicate\t$dup\n",
		"Reported\t$reported\n";

