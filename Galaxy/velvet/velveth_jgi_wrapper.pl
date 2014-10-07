#!/usr/bin/env perl

use strict;
use warnings;
my $start=time;
my $outfile=shift @ARGV;
my $outdir=shift @ARGV;

#my $option=shift @ARGV;

#print "ARG @ARGV \n";

#if ($option eq "file_input")
#{
#	open (FH, $ARGV[2]) or die "Can't open '$$ARGV[2]' : $!";
#	$ARGV[2] = <FH>;
#	close FH;
#	print $ARGV[2];
#}


my $kmer=$ARGV[2];
#die ("USER ERROR: Hash length (kmer) must be odd!\n") unless $kmer % 2;
my $tot_reads=0;
#print "@ARGV 2>&1|\n";
open (VELVETH, "@ARGV 2>&1|") or die("Unable to run velveth: $!\n");
open(OUT, ">$outfile") or die($!);
while (<VELVETH>) {
    print OUT $_;
    if (/^\[\d+\.\d+\] (\d+) sequences found/) {
        $tot_reads += $1;
    }
}
close VELVETH;
close OUT;
die("No reads found\n") unless $tot_reads;

my $sec=time-$start;
my $min=int($sec/60);
$sec -= ($min*60);
my $hr=int($min/60);
$min -= ($hr*60);
print "$tot_reads processed in";
print " $hr hr" if $hr;
print " $min min" if $min;
print " $sec sec\n";
exit

