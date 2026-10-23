#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.file.list> <PWD>\n\n";
	exit 2;
}

my $pwd = $ARGV[1];

my @path;
foreach my $file ( split /,/, $ARGV[0] ) {
	if( $file =~ /^\// ) {
		push @path, $file;
	} else {
		push @path, "$pwd/$file";
	}
}
print join(",", @path);

