#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <chromap.stderr>\n\n";
	exit 2;
}

my %stat;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	if( /Mapped all reads in ([\d\.]+)s/ ) {
		$stat{AlignmentTime} = $1 / 3600;
	} elsif( /Sorted, deduped and outputed mappings in ([\d\.]+)s/ ) {
		$stat{ReportTime} = $1 / 3600;
	} elsif( /Number of reads: (\d+)/ ) {
		$stat{Total} = $1 / 2e6;
	} elsif( /Number of mapped reads: (\d+)/ ) {
		$stat{mapped} = $1 / 2e6;
	} elsif( /Number of output mappings \(passed filters\): (\d+)/ ) {
		$stat{pass} = $1 / 1e6;
	} elsif( /, total: (\d+)/ ) {
		$stat{uniq} = $1 / 1e6;
	}
}
close IN;

print "## Running time\nRead alignment\t", $stat{AlignmentTime}, "\n",
	"Interaction-extraction\t", $stat{ReportTime}, "\n\n";

print "## Statistics\n",
	"Total\t", $stat{Total}, "\n",
    "Low-quality\t0\n",
	"Unmapped\t", $stat{Total}-$stat{mapped}, "\n",
	"Filtered\t", $stat{uniq}-$stat{pass}, "\n",
	"Duplicate\t", $stat{mapped}-$stat{uniq}, "\n",
	"Reported\t", $stat{pass}, "\n";

