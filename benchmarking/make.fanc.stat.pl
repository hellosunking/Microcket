#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "Usage: $0 <Total.read.number> <sample.stat> [sample.stat ...]\n";
	print STDERR "\nI cannot figure out how to calculate the total read number from FANC's statistics files, ",
				 "please specify this number manually.\n",
				 "The stat files are under \"FANC.output/plots/stats/\"; add all these files to the parameter list.\n";
	exit 2;
}

my $total = shift;
my ( $unmapped, $dup, $reported) = (0,0,0);

foreach my $file ( @ARGV ) {
	open IN, "$file" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/;

		if( $l[0] eq 'valid' ) {
			$reported += $l[1];
		} elsif( $l[0] eq 'PCR duplicates' ) {
			$dup += $l[1];
		} elsif( $l[0] eq 'unmappable' ) {
			$unmapped += $l[1];
		}
	}
	close IN;
}

my $filtered = $total - $unmapped - $dup - $reported;
$total /= 1e6;
$unmapped /= 1e6;
$filtered /= 1e6;
$dup /= 1e6;
$reported /= 1e6;

print	"Total\t$total\n",
		"Low-quality\t0\n",
		"Unmapped\t$unmapped\n",
		"Filtered\t$filtered\n",
		"Duplicate\t$dup\n",
		"Reported\t$reported\n";

