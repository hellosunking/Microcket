#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <genome.fa> <id> <anno.dir>\n\n";
	exit 2;
}
my $ver = '1.4';

my %genome;
my $chr = 'NA';
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	if( /^>(\S+)/ ) {
		$chr = $1;
	} else {
		chomp;
		$genome{$chr} += length $_;
	}
}
close IN;

open INFO,   ">$ARGV[2]/$ARGV[1].info" or die( "$!" );
open HEADER, ">$ARGV[2]/$ARGV[1].sam.header" or die( "$!" );

print HEADER "\@HD\tVN:1.0\tSO:coordinate\n";
foreach my $chr ( sort keys %genome ) {
	my $len = $genome{$chr};
	print INFO "$chr\t$len\n";
	print HEADER "\@SQ\tSN:$chr\tLN:$len\n";
}
print HEADER "\@PG\tID:Microcket\tPN:Microcket\tVN:$ver\tDS:$ARGV[1]\n";

close INFO;
close HEADER;

