#!/usr/bin/perl 
use strict;

#my $end_kmer = $ARGV[0];

#open FH, $ARGV[0];

#my $end_kmer = <FH>;

#close FH;

#my $step_size = $ARGV[1];


#open FH, $ARGV[2];


#my $best_kmer = <FH>;

#close FH;

#my $start_kmer = $best_kmer - 10 * $step_size;

my $option = $ARGV[0];

my $start_kmer;
my $end_kmer;
my $step_size;
my $input;
my @args;





if ($option eq "manual")
{
	$start_kmer = $ARGV[1];
	
	$end_kmer = $ARGV[2];
	
	$step_size = $ARGV[3];
	
	open FH, $ARGV[4];
	
	my $input = <FH>;
	@args = split ( /,/, $input);
}
if ($option eq "offset")
{	
	my $offset= $ARGV[1];
	
	open FH, $ARGV[2];
	
	my $input = <FH>;
	@args = split ( /,/, $input);
	my $best= @args[0];

	
	$start_kmer = $best - $offset;
	$end_kmer = $best + $offset;
	$step_size = $args[1];
}


if ($start_kmer % 2 == 0)
{
	$start_kmer = $start_kmer - 1;
}

if ($end_kmer % 2 == 0)
{
	$end_kmer = $end_kmer - 1;
}
if ($start_kmer <= 0)
{
	
	die " The first kmer $start_kmer is less then or equal to 0";
}

if ($end_kmer >= $args[2])
{
	die " The last kmer out of bounds";
}


print $start_kmer;
print ",";
print $end_kmer;
print ",";
print $step_size;
