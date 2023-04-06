#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <file1> <file2> [file3...]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $lastfile = shift;
my @start = stat( $lastfile );

foreach my $file ( @ARGV ) {
	my @here = stat( $file );
	my $runtime = $here[9] - $start[9];
	my $hour = $runtime/3600;
	print "$lastfile -> $file\t$runtime\t$hour\n";

	@start = @here;
	$lastfile = $file;
}

