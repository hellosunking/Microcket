#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <a.pairs> <b.pairs> [dist=200] [inconsistent=/dev/null]\n\n",
					"NOTE: Numbers are reported in millions.\nThis program will use LOTS OF memory.\n",
					"'a.pairs' will be hashed, use the smaller one is recommended.\n\n";
	exit 2;
}

my $maxDist = $ARGV[2] || 200;
print "Using Distance: $maxDist\n";

my $writeInconsistent = 0;
if( $#ARGV >= 3 ) {
	$writeInconsistent = 1;
	if( $ARGV[3] =~ /\.gz$/ ) {
		open WIRED, "| gzip >$ARGV[3]" or die( "$!" );
	} else {
		open WIRED, ">$ARGV[3]" or die( "$!" );
	}
}

my %pairs;
my ($A_trans, $A_cis0, $A_cis1K, $A_cis10K) = ( 0,0,0,0 );
print "Loading file A: $ARGV[0] ...\n";
if( $ARGV[0] =~ /\.gz$/ ) {
	open IN, "pigz -p 4 -cd $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##K00284:94:HNJCLBBXX:7:1117:26768:13728 chr1 181408 chr1 181477 + -
	$pairs{$l[0]} = join("\t", $l[1], $l[2], $l[3], $l[4]);

	if( $l[1] ne $l[3] ) {
		++ $A_trans;
	} else {
		my $d = abs( $l[4]-$l[2] );
		if( $d > 10000 ) {
			++ $A_cis10K;
		} elsif( $d < 1000 ) {
			++ $A_cis0;
		} else {
			++ $A_cis1K;
		}
	}
}
close IN;

my ($consistent, $dis, $a_only, $b_only) = ( 0, 0, 0, 0 );

print "Loading file B: $ARGV[1] ...\n\n";
if( $ARGV[1] =~ /\.gz$/ ) {
	open IN, "pigz -p 4 -cd $ARGV[1] |" or die( "$!" );
} else {
	open IN, "$ARGV[1]" or die( "$!" );
}
my ($B_trans, $B_cis0, $B_cis1K, $B_cis10K) = ( 0,0,0,0 );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @pb = split /\t/;
	if( $pb[1] ne $pb[3] ) {
		++ $B_trans;
	} else {
		my $d = abs( $pb[4]-$pb[2] );
		if( $d > 10000 ) {
			++ $B_cis10K;
		} elsif( $d < 1000 ) {
			++ $B_cis0;
		} else {
			++ $B_cis1K;
		}
	}

	my $id = shift @pb;
	if( exists $pairs{$id} ) {
		my @pa = split /\t/, $pairs{$id};

		if( ($pa[0] eq $pb[0] && abs($pa[1]-$pb[1])<$maxDist && 
			 $pa[2] eq $pb[2] && abs($pa[3]-$pb[3])<$maxDist) || 
			($pa[0] eq $pb[2] && abs($pa[1]-$pb[3])<$maxDist &&
			 $pa[2] eq $pb[0] && abs($pa[3]-$pb[1])<$maxDist) ) {
			++ $consistent;
		} else {
			++ $dis;
			$#pb = 3;	## remove the tailing strands
			print WIRED join("\t", "I", $id, $pairs{$id}, "vs", @pb), "\n" if $writeInconsistent;
		}
		delete $pairs{$id};
	} else {
		++ $b_only;
		$#pb = 3;	## remove the tailing strands
		print WIRED join("\t", 'B', $id, @pb), "\n" if $writeInconsistent;
	}
}
$a_only = keys %pairs;

if( $writeInconsistent ) {
	foreach my $id ( keys %pairs ) {
		print WIRED join("\t", 'A', $id, $pairs{$id}), "\n";
	}
	close WIRED;
}

print "Overalp report:\n",
	  "A_only\t", ($a_only+$dis)/1e6, "\n",
	  "B_only\t", ($b_only+$dis)/1e6, "\n",
	  "Consistent\t", $consistent/1e6, "\n",
	  "#Disconcordant\t", $dis/1e6, "\n\n",

	  "Pairs distribution report (file A):\n",
	  "Cis (<1 K)\t", $A_cis0/1e6, "\n",
	  "Cis (1-10 K)\t", $A_cis1K/1e6, "\n",
	  "Cis (>10 K)\t", $A_cis10K/1e6, "\n",
	  "Trans\t", $A_trans/1e6, "\n\n",

	  "Pairs distribution report (file B):\n",
	  "Cis (<1 K)\t", $B_cis0/1e6, "\n",
	  "Cis (1-10 K)\t", $B_cis1K/1e6, "\n",
	  "Cis (>10 K)\t", $B_cis10K/1e6, "\n",
	  "Trans\t", $B_trans/1e6, "\n\n";

## The Disconcordant are included into A-only and B-only

