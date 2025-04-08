#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.fa[.gz]>";
	print STDERR "\nThis program is designed to remove non-primary contigs thing in the genome.\n\n";
	exit 2;
}

if( $ARGV[0] =~ /\.gz$/ ) {
	open IN, "zcat $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}

my $keepflag = 0;
my $r;
while( $r = <IN> ) {
	if( $r =~ /^>(\S+)/ ) {
		my $chr = $1;
		if( $chr!~/^chr/ || $chr=~/_/ || $chr=~/Un/ || $chr=~/rand/ ) {
			$keepflag = 0;
		} else {
			$keepflag = 1;
			print $r;
		}
	} else {
		print $r if $keepflag;
	}
}
close IN;

