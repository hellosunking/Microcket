#!/usr/bin/perl

## This proram is adapted from "statistics.pl" in Juicer software (v1.6)
## Designed to filter "merged_nodups.txt", convert the passed data to ".pairs" format, and add some statistics

use File::Basename;
use POSIX;
use List::Util qw/min max/;
use Getopt::Std;
use vars qw/$opt_s $opt_l $opt_d $opt_o $opt_q $opt_p $opt_h/;

# Check arguments
getopts('s:l:o:q:p:h');

my $site_file = "";
my $ligation_junction = "GATCGATC";
my $stats_file = "stats.txt";
my $mapq_threshold = 30;
my $output = '/dev/stdout';

if ($opt_s) {
  $site_file = $opt_s;
}
if ($opt_l) {
  $ligation_junction = $opt_l;
}
if ($opt_o) {
  $stats_file = $opt_o;
}
if ($opt_q) {
  $mapq_threshold = $opt_q;
}
if ($opt_p) {
	$output = $opt_p;
}

if (scalar(@ARGV)==0 || $opt_h) {
    print STDERR "Usage: statistics.pl -s [site file] -l [ligation] -o [stats file] -q [mapq threshold] -p [out pairs] <infile>\n";
    print STDERR " <infile>: file in intermediate format to calculate statistics on, set '-' to use stdin stream\n";
    print STDERR " [site file]: list of HindIII restriction sites, one line per chromosome (default NULL)\n";
    print STDERR " [ligation]: ligation junction (default $ligation_junction)\n";
    print STDERR " [stats file]: output file containing total reads, for library complexity (default $stats_file)\n";
    print STDERR " [mapq threshold]: mapping quality threshold, do not consider reads < threshold (default $mapq_threshold)\n";
    print STDERR " [out pairs]: filtered interactions in pairs format (default stdout)\n";
    exit 2;
}

# Remove parenthesis
$ligation_junction =~ s/\(//;
$ligation_junction =~ s/\)//;

# with OR | symbol, this won't work, need to explicitly fix
my $dangling_junction = substr $ligation_junction, length($ligation_junction)/2;

# Global variables for calculating statistics
my %chromosomes;
my %hindIII;
my %mapQ;
my %mapQ_inter;
my %mapQ_intra;
my %innerM;
my %outerM;
my %rightM;
my %leftM;
my $three_prime_end=0;
my $five_prime_end=0;
my $total = 0;
my $dangling = 0;
my $ligation = 0;
my $inner = 0;
my $outer = 0;
my $left = 0;
my $right = 0;
my $inter = 0;
my $intra = 0;
my $small = 0;
my $large = 0;
my $short0K = 0;
my $inter10K = 0;
my $large10K = 0;
my $very_small = 0;
my $very_small_dangling = 0;
my $small_dangling = 0;
my $large_dangling = 0;
my $inter_dangling = 0;
my $true_dangling_intra_small = 0;
my $true_dangling_intra_large = 0;
my $true_dangling_inter = 0;
my $total_current = 0;
my $under_mapq = 0;
my $intra_fragment = 0;
my $unique = 0;
# logspace bins
my @bins = (10,12,15,19,23,28,35,43,53,66,81,100,123,152,187,231,285,351,433,534,658,811,1000,1233,1520,1874,2310,2848,3511,4329,5337,6579,8111,10000,12328,15199,18738,23101,28480,35112,43288,53367,65793,81113,100000,123285,151991,187382,231013,284804,351119,432876,533670,657933,811131,1000000,1232847,1519911,1873817,2310130,2848036,3511192,4328761,5336699,6579332,8111308,10000000,12328467,15199111,18738174,23101297,28480359,35111917,43287613,53366992,65793322,81113083,100000000,123284674,151991108,187381742,231012970,284803587,351119173,432876128,533669923,657933225,811130831,1000000000,1232846739,1519911083,1873817423,2310129700,2848035868,3511191734,4328761281,5336699231,6579332247,8111308308,10000000000);

if(length($site_file==0) || index($site_file, "none") != -1) {
   #no restriction enzyme, no need for RE distance
} else {
  # read in restriction site file and store as multidimensional array
  open FILE, "$site_file" or die $!;
  while (<FILE>) {
    my @locs = split;
    my $key = shift(@locs);
    my $ref = \@locs;
    $chromosomes{$key} = $ref;
  }
  close(FILE);
}
# read in infile and calculate statistics
#open FILE, $infile or die $!;

open PAIRS, ">$output" or die $!;
print PAIRS "## pairs format v1.0\n#columns: readID chr1 position1 chr2 position2 strand1 strand2\n";

while (<>) {
  $unique++;
	my @record = split;	##0 chr1 10013 0 16 chr1 10442 0 0 101M TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCAAACCCAAACC 0 17S36M1D48M CCCCCAACCCCAACCCCAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC SRR1658610.sra.16 SRR1658610.sra.16
	my $num_records = scalar(@record);
  # don't count as Hi-C contact if fails mapq or intra fragment test
  my $countme = 1;

  if (($record[1] eq $record[5]) && $record[3] == $record[7]) {
    $intra_fragment++;
    $countme = 0;
  }
	elsif ($num_records > 8) {
		my $mapq_val = min($record[8],$record[11]);
    if ($mapq_val < $mapq_threshold) {
      $under_mapq++;
      $countme = 0;
    }
  }


  if ($countme) {
    $total_current ++;
	print PAIRS join("\t", $record[14], $record[1],$record[2], $record[5],$record[6], ($record[0] & 16)?'-':'+',($record[4] & 16)?'-':'+'), "\n";

    # position distance
    my $pos_dist = abs($record[2] - $record[6]);
    
    my $hist_dist = &bsearch($pos_dist,\@bins); 
    
    my $is_dangling = 0;
    # one part of read pair has unligated end
    if ($num_records > 8 && ($record[10] =~ m/^$dangling_junction/ || $record[13] =~ m/^$dangling_junction/)) {
      $dangling++;
      $is_dangling=1;
    }
    # look at chromosomes
    if ($record[1] eq $record[5]) {	## cis-
      $intra++;
	  if( $pos_dist > 10000 ) {
		  ++ $large10K;
	  } elsif( $pos_dist < 1000 ) {
		  ++ $short0K;
	  } else {
		  ++ $inter10K;
	  }
      # determine right/left/inner/outer ordering of chromosomes/strands
      if ($record[0] == $record[4]) {
        if ($record[0] == 0) {
          if ($pos_dist >= 20000) {
            $right++;
          }
          $rightM{$hist_dist}++;
        }else {
          if ($pos_dist >= 20000) {
            $left++;
          }
          $leftM{$hist_dist}++;
        }
      }else {
        if ($record[0] == 0) {
          if ($record[2] < $record[6]) {
            if ($pos_dist >= 20000) {
              $inner++;
            }
            $innerM{$hist_dist}++;
          } else {
            if ($pos_dist >= 20000) {
              $outer++;
            }
            $outerM{$hist_dist}++;
          }
        } else {
          if ($record[2] < $record[6]) {
            if ($pos_dist >= 20000) {
              $outer++;
            }
            $outerM{$hist_dist}++;
          } else {
            if ($pos_dist >= 20000) {
              $inner++;
            }
            $innerM{$hist_dist}++;
          }
        }
      }
      # intra reads less than 20KB apart
      if ($pos_dist < 10) {
        $very_small++;
        if ($is_dangling) {
          $very_small_dangling++;
        }
      } elsif ($pos_dist < 20000) {
        $small++;
        if ($is_dangling) {
          $small_dangling++;
        }
      } else {
        $large++;
        if ($is_dangling) {
          $large_dangling++;
        }
      }
    } else {
      $inter++;
      if ($is_dangling) {
        $inter_dangling++;
      }
    }
    if ($num_records > 8) {
      my $mapq_val = min($record[8],$record[11]);
      if ($mapq_val <= 200) {
        $mapQ{$mapq_val}++;
        if ($record[1] eq $record[5]) {
          $mapQ_intra{$mapq_val}++;
        }
        else {
          $mapQ_inter{$mapq_val}++;
        }
      }
      # read pair contains ligation junction
      if ($record[10] =~ m/($ligation_junction)/ || $record[13] =~ m/($ligation_junction)/) {
        $ligation++;
      }
    }
    # determine distance from nearest HindIII site, add to histogram
    if (index($site_file, "none") == -1) {
      my $report = (($record[1] != $record[5]) || ($pos_dist >= 20000));
      my $dist = &distHindIII($record[0], $record[1], $record[2], $record[3], $report);
      if ($dist <= 2000) {
        $hindIII{$dist}++;
      }

      $dist = &distHindIII($record[4], $record[5], $record[6], $record[7], $report);
      if ($dist <= 2000) {
        $hindIII{$dist}++;
      }
    }
    if ($is_dangling) {
      if ($record[10] =~ m/^$dangling_junction/) {
        $dist = &distHindIII($record[0], $record[1], $record[2], $record[3], 1);
      }
      else { #	$record[13] =~ m/^$dangling_junction/) 
        $dist = &distHindIII($record[4], $record[5], $record[6], $record[7], 1);
      }
      if ($dist == 1) {
        if ($record[1] == $record[5]) {
          if ($pos_dist < 20000) {
            $true_dangling_intra_small++;
          }
          else {
            $true_dangling_intra_large++;
          }
        }
        else {
          $true_dangling_inter++;
        }
      }
    }
  }
 }
close PAIRS;

my($statsfilename, $directories, $suffix)= fileparse($stats_file, qr/\.[^.]*/);
my $histsfile = $directories . $statsfilename . "_hists.m";

my $seq=0;
if (-e $stats_file) {
  open FILE, $stats_file or die $!;
  while (<FILE>) {
    if (/Sequenced/g) {
      ($label, $reads) = split(':', $_);
      $seq=1;
      $reads =~ s/,//g;
      $reads =~ s/ //g;
    }
  } 
  close FILE;
}
open FILE, " >> $stats_file", or die $!;
if ($unique==0) {
	$unique=1;
}

print FILE "Intra-fragment Reads: " . commify($intra_fragment);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $intra_fragment*100/$reads; 
} else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $intra_fragment*100/$unique; 

print FILE "Below MAPQ Threshold: " . commify($under_mapq);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $under_mapq*100/$reads; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $under_mapq*100/$unique; 

print FILE "Hi-C Contacts: " . commify($total_current);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $total_current*100/$reads; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $total_current*100/$unique; 

printf FILE " Ligation Motif Present: %s ", commify($ligation);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $ligation*100/$reads; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $ligation*100/$unique; 

if ($five_prime_end + $three_prime_end > 0) {
  my $f1 = $three_prime_end*100/($five_prime_end + $three_prime_end);
  my $f2 = $five_prime_end*100/($five_prime_end + $three_prime_end);
  printf FILE " 3' Bias (Long Range): %0.0f\%", $f1;
  printf FILE " - %0.0f\%\n", $f2;
}else {
  print FILE " 3' Bias (Long Range): 0\% \- 0\%\n";
}
if ($large > 0) {
  printf FILE " Pair Type %(L-I-O-R): %0.0f\%", $left*100/$large;
  printf FILE " - %0.0f\%", $inner*100/$large;
  printf FILE " - %0.0f\%", $outer*100/$large;
  printf FILE " - %0.0f\%\n", $right*100/$large;
}else {
	print FILE " Pair Type %(L-I-O-R): 0\% - 0\% - 0\% - 0\%\n";
}

printf FILE "Inter-chromosomal: %s ", commify($inter);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $inter*100/$reads; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $inter*100/$unique; 

printf FILE "Intra-chromosomal: %s ", commify($intra);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $intra*100/$reads; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $intra*100/$unique; 

printf FILE "-Short Range (<20Kb): %s ", commify($small);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $small*100/$intra; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $small*100/$unique; 

printf FILE "-Long Range (>20Kb): %s ", commify($large);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $large*100/$intra; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $large*100/$intra;

## added distance statistics
printf FILE "=Short Range (<1Kb): %s ", commify($short0K);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $short0K*100/$intra;
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $short0K*100/$intra;

printf FILE "=Intermediate Range (1-10Kb): %s ", commify($inter10K);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $inter10K*100/$intra; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $inter10K*100/$intra;

printf FILE "=Long Range (>10Kb): %s ", commify($large10K);
if ($seq == 1) {
  printf FILE " (%0.2f\% / ", $large10K*100/$intra; 
}else {
  print FILE "(";
}
printf FILE "%0.2f\%)\n", $large10K*100/$intra;

close FILE;

open FILE, "> $histsfile", or die $!; 
print FILE "A = [\n";
for (my $i=1; $i <= 2000; $i++) {
	my $tmp =	 $hindIII{$i} || 0;
	print FILE "$tmp ";
}
print FILE "\n];\n";
print FILE "B = [\n";
for (my $i=0; $i <= 200; $i++) {
	my $tmp = $mapQ{$i} || 0;
	my $tmp2 = $mapQ_intra{$i} || 0;
	my $tmp3 = $mapQ_inter{$i} || 0;
	print FILE "$tmp $tmp2 $tmp3\n ";
}
print FILE "\n];\n";
print FILE "D = [\n";
for (my $i=0; $i < scalar(@bins); $i++) {
	my $tmp = $innerM{$i} || 0;
	print FILE "$tmp ";
	$tmp = $outerM{$i} || 0;
	print FILE "$tmp ";
	$tmp = $rightM{$i} || 0;
	print FILE "$tmp ";
	$tmp = $leftM{$i} || 0;
	print FILE "$tmp\n";
}

print FILE "\n];";
print FILE "x = [\n";
for (my $i=0; $i < scalar(@bins); $i++) {
	print FILE "$bins[$i] ";
}
print FILE "\n];\n";
close FILE;

# Find distance to nearest HindIII restriction site
sub distHindIII {
	# find upper index of position in sites array via binary search
	my $index = $_[3];
	# get distance to each end of HindIII fragment
	my $dist1;
	if ($index == 0) {
		# first fragment, distance is position
		$dist1 =	$_[2];	
	}
	else {
		$dist1 = abs($_[2] - $chromosomes{$_[1]}[$index-1]);
	}
	my $dist2 = abs($_[2] - $chromosomes{$_[1]}[$index]);
	
	# get minimum value -- if (dist1 <= dist2), it's dist1, else dist2
	my $retval = $dist1 <= $dist2 ? $dist1 : $dist2; 
	# get which end of the fragment this is, 3' or 5' (depends on strand)
	if ($retval == $dist1 && $_[4]) {
		$_[0] == 0 ? $five_prime_end++ : $three_prime_end++;
	}
	elsif ($retval == $dist2 && $_[4]) {
		$_[0] == 16 ? $five_prime_end++ : $three_prime_end++;
	}
	return $retval;
}

# Binary search, array passed by reference
# search array of integers a for given integer x
# return index where found or upper index if not found
sub bsearch {
	my ($x, $a) = @_;		 # search for x in array a
	my ($l, $u) = (0, @$a - 1);	 # lower, upper end of search interval
	my $i;	      		 # index of probe
	while ($l <= $u) {
		$i = int(($l + $u)/2);
		if ($a->[$i] < $x) {
	    $l = $i+1;
		}
		elsif ($a->[$i] > $x) {
	    $u = $i-1;
		}
		else {
	    return $i; # found
		}
	}
	return $l;				 # not found, return upper
}
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
