#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <sim3C.fastq> <out.prefix>\n\n";
	exit 2;
}

#@SIM3C:1684747613:3C:1:1:1:1 1:Y:18:1 HIC chr7:72956331 chr7:66024672
#CGGAGGTTGCAGTAAGCTGAGATTGTGCCACCGCACTCCAGCCTGGGGGATAGAGCAAGCTAGCTAGCTTTTGAAGTTTATCAGCTCTAGAACTGACTCTAAAAT
#+
#AAFFFFKKKKKKKKKKKKKKKKKKKKKKKK,KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKFAKKKKKKKKKKKFKKKKKKKKKKKKKKKKKKKKKKKKKKKKF
#@SIM3C:1684747613:3C:1:1:1:1 2:Y:18:1 HIC chr7:72956331 chr7:66024672
#AAATCTTTTCTTGGAGTTCCAGTTTGTATGACCTATCCCATCACCCCCACTCCCAGAATTCAGACCAAATTGGATGAAAAGTAGGGTTGACAAGTTAAAGTACTT
#+
#AAF,FKKKKKKAKKKKKKKKKKKFKKKKKKKKKKK(KKK,KKKKKKKK,KFK7KKKKK,KKKKKKKKK7KKKAKKKF,KFKKKKKKKKKKKKK<KF<FKKAKKAK

my $fastq  = shift;
my $prefix = shift;

open R1, "| gzip >${prefix}_R1.fastq.gz" or die( "$!" );
open R2, "| gzip >${prefix}_R2.fastq.gz" or die( "$!" );

my $cnt = 0;
open IN, "$fastq" or die( "$!" );
while( my $r1 = <IN> ) {
	my $s1 = <IN>;
	<IN>;
	my $q1 = <IN>;

	my $r2 = <IN>;
	my $s2 = <IN>;
	<IN>;
	my $q2 = <IN>;

	chomp( $r1 );
	my @l = split /\s+/, $r1, 4;
	$l[3] =~ s/\s+/-/g;

	++$cnt;
	print R1 "\@$cnt#$l[3]\n$s1+\n$q1";
	print R2 "\@$cnt#$l[3]\n$s2+\n$q2";
}
close IN;

close R1;
close R2;

