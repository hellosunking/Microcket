#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;
#use KSLIB::Digitalize qw/digitalize/;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <sid> <concat=yes|no>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $sid    = $ARGV[0];
my $concat = $ARGV[1];

print "#Category\tCount\tFraction(%)\n";
## preprocess
my %trim;
open IN, "cat $sid.trim.log |" or die( "$!" );	## in case there are many files in biological replicates
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##Total   60802850
	$trim{$l[0]} += $l[1];
}
close IN;

my %rmdup;
open IN, "$sid.rmdup.log" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	$rmdup{$l[0]} += $l[1];
}
close IN;

print "## Preprocessing and alignment\n";
printf "Total\t%s\t100.0\nKtrim\t%s\t%.1f\nUnique\t%s\t%.1f\n",
		d($trim{"Total"}),
		d($rmdup{"Total"}), $rmdup{"Total"}/$trim{"Total"}*100,
		d($rmdup{"Uniq"}), $rmdup{"Uniq"}/$rmdup{"Total"}*100;

## flash
my $prealign;
if( $concat eq "yes" ) {
	my ($cat, $unc, $cut) = ( 0, 0, 0 );

	if( -s "$sid.flash.log" ) {	## old version
		open IN, "$sid.flash.log" or die( "$!" );
		while( <IN> ) {
			chomp;
			if( /\sCombined pairs:\s+(\d+)/ ) {
				$cat = $1;
				last;
			}
		}
		close IN;

		if( -s "$sid.cut.log" ) {
			open IN, "$sid.cut.log" or die( "$!" );
			while( <IN> ) {
				chomp;
				if( /Total\s+(\d+)/ ) {$unc = $1;}
				if( /Pass\s+(\d+)/  ) {$cut = $1;}
			}
			close IN;
		} else {
			$unc = $rmdup{"Uniq"} - $cat;
			$cut = $unc;
		}
	} else {	## new version
		open IN, "$sid.stitch.stat" or die( "$!" );
		my $line = <IN>;chomp($line);my @l = split /\t/, $line;
		$cat = $l[1];
		$unc = $l[3];
		$cut = $l[5];
		close IN;
	}

	printf "Stitched\t%s\t%.1f\nUnstitched\t%s\t%.1f\n  Discarded(too-short)\t%s\t%.1f\n",
			d($cat), $cat/$rmdup{"Uniq"}*100,
			d($cut), $cut/$rmdup{"Uniq"}*100,
			d($unc-$cut), ($unc-$cut)/$rmdup{"Uniq"}*100;

	$prealign = $cat + $cut;
} else {
	$prealign = $rmdup{"Uniq"};
}

## alignment
my %align;
if( $concat eq "yes" ) {
	open IN, "$sid.flash2pairs.log" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/;
		$align{$l[0]} += $l[1];
		$align{"all"} += $l[1];
	}
	close IN;
}
open IN, "$sid.unc2pairs.log" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	$align{$l[0]} += $l[1];
	$align{"all"} += $l[1];
}
close IN;

printf "Mappable\t%s\t%.1f\n", d($align{"all"}), $align{"all"}/$prealign*100;

print "## Interactions\n";
my $uncalled = $align{"lowMap"}+$align{"manyHits"}+$align{"unpaired"}+$align{"selfCircle"};
printf "Uncalled\t%s\t%.1f\n", d($uncalled), $uncalled/$align{"all"}*100;
printf "  Incomplete-mapping\t%s\t%.1f\n", d($align{"lowMap"}), $align{"lowMap"}/$align{"all"}*100;
printf "  Too-many-segments\t%s\t%.1f\n", d($align{"manyHits"}), $align{"manyHits"}/$align{"all"}*100;
printf "  Unpairable\t%s\t%.1f\n",  d($align{"unpaired"}), $align{"unpaired"}/$align{"all"}*100;
printf "  Self-circle\t%s\t%.1f\n", d($align{"selfCircle"}), $align{"selfCircle"}/$align{"all"}*100;

my $valid = $align{"trans"} + $align{"cis10K"} + $align{"cis1K"} + $align{"cis0"};
printf "Valid-pair\t%s\t%.1f\n", d($valid), $valid/$align{"all"}*100;
printf "  Cis(<1K)\t%s\t%.1f\n", d($align{"cis0"}), $align{"cis0"}/$valid*100;
printf "  Cis(1-10K)\t%s\t%.1f\n", d($align{"cis1K"}), $align{"cis1K"}/$valid*100;
printf "  Cis(>=10K)\t%s\t%.1f\n", d($align{"cis10K"}), $align{"cis10K"}/$valid*100;
printf "  Trans\t%s\t%.1f\n", d($align{"trans"}), $align{"trans"}/$valid*100;

## digitalize a given number
sub d {
	my $v = shift;
	while($v =~ s/(\d)(\d{3})((:?,\d\d\d)*)$/$1,$2$3/){};
	return $v;
}

