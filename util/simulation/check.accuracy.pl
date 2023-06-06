#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "Usage: $0 <in.pairs[.gz]>\n";
	print STDERR "This program is designed to check the alignment accuracy on sim3C reads.\n";
	exit 2;
}

my $maxDist = 500;	## this is due to enzymatic cutting limit
my ($pass, $err) = ( 0, 0 );

open IN, "less $ARGV[0] |" or die( "$!" );
while( my $r = <IN> ) {
	next if $r =~ /^#/;

#	print STDERR "Checking $r";
	my @l = split /\t/, $r;
	$l[0] =~ s/.*#//;
	my ($c1, $p1, $c2, $p2, $extra);
	if( $l[0] =~ /-/) { ##chr1:38369-chr1:211787796	chr1	38116	chr1	38373	+	-
		($c1, $p1, $c2, $p2, $extra) = split /[:-]/, $l[0];
	} else { ##chr1:10007..10439:F	chr1	10008	chr1	10209	+	-
		$l[0] =~ s/\.\./-/;
		($c1, $p1, $p2, $extra) = split /[:-]/, $l[0];
		$c2 = $c1;
	}
#	print STDERR "$c1:$p1 + $c2:$p2 vs $l[1]:$l[2] + $l[3]:$l[4]\n";
	if( $c1 eq $c2 ) {
		if( $c1 eq $l[1] && $c2 eq $l[3] ) {
			if( abs($l[2]-$p1)<$maxDist && abs($l[4]-$p2)<$maxDist ) {
				++ $pass;
			} elsif( abs($l[2]-$p2)<$maxDist && abs($l[4]-$p1)<$maxDist ) {
				++ $pass;
			} elsif( abs($l[2]-$p1)<$maxDist && abs($l[4]-$p1)<$maxDist ) {
				## sometimes the 3D connection loci could not be fully recovered and ONLY 1 part is correct
				++ $pass;
			} elsif( abs($l[2]-$p2)<$maxDist && abs($l[4]-$p2)<$maxDist ) {
				++ $pass;
			} else {
				++ $err;
#				print STDERR "$r";
			}
		} else {
			++ $err;
#			print STDERR "$r";
		}
	} else {
		if( $c1 eq $l[1] && $c2 eq $l[3] ) {
			if( abs($l[2]-$p1)<$maxDist && abs($l[4]-$p2)<$maxDist ) {
				++ $pass;
			} else {
				++ $err;
#				print STDERR "$r";
			}
		} elsif( $c1 eq $l[3] && $c2 eq $l[1] ) {
			if( abs($l[2]-$p2)<$maxDist && abs($l[4]-$p1)<$maxDist ) {
				++ $pass;
			} else {
				++ $err;
#				print STDERR "$r";
			}
		} elsif ( $c1 eq $l[3] && $c1 eq $l[1] && abs($l[2]-$p1)<$maxDist && abs($l[4]-$p1)<$maxDist ) {
			## sometimes the 3D connection loci could not be fully recovered and ONLY 1 part is correct
			++ $pass;
		} elsif ( $c2 eq $l[3] && $c2 eq $l[1] && abs($l[2]-$p2)<$maxDist && abs($l[4]-$p2)<$maxDist ) {
			++ $pass;
		} else {
			++ $err;
#			print STDERR "$r";
		}
	}
}
close IN;

my $total = $pass + $err;
my $pr = $pass/$total*100;
my $er = $err/$total*100;
print "File\t$ARGV[0]\nTotal\t$total\t100%\nCorrect\t$pass\t$pr%\nError\t$err\t$er%\n";
 
