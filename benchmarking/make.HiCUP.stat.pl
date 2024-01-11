#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <hicup_summary_report.txt>\n\n";
	exit 2;
}

my ($total, $lowQual, $mapped, $filtered, $reported) = (0,0,0,0,0);

open IN, "$ARGV[0]" or die( "$!" );
<IN>;	## header
#File	Total_Reads_1	Total_Reads_2	Not_Truncated_Reads_1	Not_Truncated_Reads_2	Truncated_Read_1	Truncated_Read_2	Average_Length_Truncated_1	Average_Length_Truncated_2	Too_Short_To_Map_Read_1	Too_Short_To_Map_Read_2	Unique_Alignments_Read_1	Unique_Alignments_Read_2	Multiple_Alignments_Read_1	Multiple_Alignments_Read_2	Failed_To_Align_Read_1	Failed_To_Align_Read_2	Paired_Read_1	Paired_Read_2	Valid_Pairs	Valid_Cis_Close	Valid_Cis_Far	Valid_Trans	Invalid_Pairs	Same_Circularised	Same_Dangling_Ends	Same_Fragment_Internal	Re_Ligation	Contiguous_Sequence	Wrong_Size	Deduplication_Read_Pairs_Uniques	Deduplication_Cis_Close_Uniques	Deduplication_Cis_Far_Uniques	Deduplication_Trans_Uniques	Percentage_Mapped	Percentage_Valid	Percentage_Uniques	Percentage_Unique_Trans	Percentage_Ditags_Passed_Through_HiCUP
while( <IN> ) {
	chomp;
	my @l = split /\s+/;

	$total += $l[1];
	$lowQual += ($l[9] > $l[10]) ? $l[9] : $l[10];
	$mapped += $l[17];
	$filtered += $l[23];
	$reported += $l[30];
}
close IN;

$total /= 1e6;
$lowQual /= 1e6;
$mapped /= 1e6;
$filtered /= 1e6;
$reported /= 1e6;

print   "Total\t$total\n",
		"Low-quality\t", $lowQual, "\n",
		"Unmapped\t", $total-$lowQual-$mapped, "\n",
		"Filtered\t$filtered\n",
		"Duplicate\t", $mapped-$filtered-$reported, "\n",
		"Reported\t$reported\n";

