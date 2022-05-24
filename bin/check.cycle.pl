#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.fq[.gz]> [read.num=10000] [fraction=0.25]\n\n";
	exit 2;
}

my $readNum  = $ARGV[1] || 10000;	## only analyze the top reads
my $fraction = $ARGV[2] || 0.25;	## report the 1/4th longest cycle, in case the input reads are nonstandard

if( $readNum<=0 || $fraction>1 ) {
	print STDERR "ERROR: incorrect parameter!\n";
	print -1;
	exit 10;
}

if( $ARGV[0] =~ /\.gz$/ ) {
	open IN, "zcat $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}
## skip the heading 1000 reads, which usually have poor quality
my $curr = 0;
while( <IN> ) {
	<IN>;
	<IN>;
	<IN>;

	++$curr;
	last if $curr == 1000;
}

my @size;
$curr = 0;
while( <IN> ) {
	my $seq = <IN>;
	chomp( $seq );
	<IN>;
	<IN>;

	push @size, length($seq);
	++ $curr;
	last if $curr == $readNum;
}
close IN;

if( $curr == 0 ) {	## no valid reads
	print STDERR "ERROR: NO enough reads!\n";
	print 0;
	exit 1;
}

if( $curr != $readNum ) {	## read not enough
	print STDERR "WARNING: NO enough reads!\n";
}

my @sorted = sort {$a<=>$b} @size;
-- $curr;
my $pos = int($curr * $fraction);
print $sorted[$pos];

