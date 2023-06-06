#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <fq.list> [read.num=10000] [fraction=0.25]\n\n";
	exit 2;
}

my $readNum  = $ARGV[1] || 10000;	## only analyze the top reads
my $fraction = $ARGV[2] || 0.25;	## report the 1/4th longest cycle, in case the input reads are nonstandard

if( $readNum<=0 || $fraction>1 ) {
	print STDERR "ERROR: incorrect parameter!\n";
	print -1;
	exit 10;
}

my @files;
my $err = 0;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	## R1 R2 extra
	my @R1 = split /,/, $l[0];
	my @R2 = split /,/, $l[1];
	if( $#R1 != $#R2 ) {
		print STDERR "ERROR: file number does not match!\n";
		$err = 1;
		last;
	}

	foreach my $i ( @R1, @R2 ) {
		if( -s $i ) {	## check whether file exist
			push @files, $i;
		} else {
			print STDERR "ERROR: file $i does not exist!\n";
			$err = 1;
			last;
		}
	}
}
close IN;

exit 1 if $err;

if( $files[0] =~ /\.gz$/ ) {
	open IN, "zcat $files[0] |" or die( "$!" );
} else {
	open IN, "$files[0]" or die( "$!" );
}

## skip the heading 1000 reads, which usually have poor quality
my $curr = 0;
while( <IN> ) {
	<IN>;<IN>;<IN>;

	++ $curr;
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
	print STDERR "ERROR: NO reads loaded!\n";
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

