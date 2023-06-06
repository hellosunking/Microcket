#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <flash.tab.out> <out.prefix> [cut.tail=10] [min.size=36]\n\n";
	exit 2;
}

my $cutTail = $ARGV[2] || 10;
my $minSize = $ARGV[3] || 36;
$minSize += $cutTail;

open EXT, ">$ARGV[1].ext.fq" or die( "$!" );
open UN1, ">$ARGV[1].cut.read1.fq" or die( "$!" );
open UN2, ">$ARGV[1].cut.read2.fq" or die( "$!" );

my ($ext, $unc, $pass) = ( 0, 0 );
open IN, "$ARGV[0]" or die( "$!" );
my @l;
while( <IN> ) {
	chomp;
	@l = split /\t/;	## sid seq1 qual1 [seq2 qual2]

	if( $#l == 2 ) {	## combined
		++ $ext;
		print EXT "\@$l[0]\n$l[1]\n+\n$l[2]\n";
	} else {	## uncombined
		++ $unc;
		if( length($l[1]) >= $minSize ) {
			++ $pass;
			$l[1] = substr( $l[1], 0, length($l[1])-$cutTail );
			$l[2] = substr( $l[2], 0, length($l[2])-$cutTail );
			$l[3] = substr( $l[3], 0, length($l[3])-$cutTail );
			$l[4] = substr( $l[4], 0, length($l[4])-$cutTail );

			print UN1 "\@$l[0]\n$l[1]\n+\n$l[2]\n";
			print UN2 "\@$l[0]\n$l[3]\n+\n$l[4]\n";
		}
	}
}
close IN;

close EXT;
close UN1;
close UN2;

# statistics
open OUT, ">$ARGV[1].stitch.stat" or die( "$!" );
print OUT "Combined\t$ext\tUncombined\t$unc\tPass\t$pass\n";
close OUT;

