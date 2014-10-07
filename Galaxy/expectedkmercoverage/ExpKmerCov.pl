#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);

# Bug: 2481
# Input: [read lengths, num reads], kmer (val or range), 
# OR: input reads files, kmer (val or range), 
# AND: est. genome size or contigs file
# Output: exp_cov values for each input kmer.

my $options = {};
GetOptions ($options,
    'read_length|l=i@',
    'num_reads|n=i@',
    'reads_file|f=s@',
    'kmer_range|r=s@',
    'kmer_list|k=s@',
    'genome_size|g=i',
    'contig_file|f=s',
    );

# Currently only works with read lengths, num reads, and genome size (e.g. no file parsing).

# From a set of kmer ranges of form e.g. 41,59,2
# and kmer values of form 31,39,47,...
# create a single array of kmer values to use.
sub get_kmer_list
{
    my %uniq_kmers = ();
    for my $kr (@{$options->{kmer_range}}) {
        my @fields = split(/,/, $kr);
        if (scalar @fields == 3) {
            my ($min, $max, $inc) = @fields;
            if ($min % 2 == 1 and $max % 2 == 1 and $inc % 2 == 0 and $min <= $max and $min > 0) {
                for (my $i = $min; $i <= $max; $i += $inc) {
                    $uniq_kmers{$i}++;
                }
            } else {
                die "kmer range $kr is invalid - require odd min,max, even increment, and min<=max\n";
            }
        } else {
            die "kmer range $kr is invalid - must have three entries: min,max,inc\n";
        }
    }
    for my $kl (@{$options->{kmer_list}}) {
        my @fields = split(/,/,$kl);
        foreach my $kval (@fields) {
            if ($kval =~ /^\d+$/ and $kval %2 == 1 and $kval > 0) {
                $uniq_kmers{$kval}++;
            } else {
                print "omitting non-integer kmer '$kval'\n";
            }
        }
    }
    # return numerically sorted list of unique kmer values
    my $kmer_list = [sort {$a <=> $b} keys %uniq_kmers];
    return $kmer_list;
}
            

sub calc_common_stats
{
    my @rls = @{$options->{read_length}};
    my @nrs = @{$options->{num_reads}};
    my ($total_cov, $avg_readlen) = ('','');
    if (scalar (@rls) != scalar (@nrs)) {
        die "Error: number of values differs between read lengths and num reads\n";
    } else {
        my $num_files = scalar @rls;
        my @products = ();
        for (my $i=0; $i<$num_files; $i++) {
            push (@products, ($rls[$i] * $nrs[$i]));
        }
        my $total_bases = sum (@products);
        $avg_readlen = $total_bases / sum (@nrs);
        $total_cov = $total_bases / $options->{genome_size};
    }
    return ($total_cov, $avg_readlen);
}

sub calc_exp_cov
{
    my $total_cov = shift;
    my $avg_readlen = shift;
    my $kmer = shift;
    my $exp_cov_float = -0.5;
    if ($avg_readlen) {
        $exp_cov_float = $total_cov * ($avg_readlen - $kmer + 1) / $avg_readlen;
    }
    my $exp_cov_rounded_int = int($exp_cov_float + 0.5);
    return $exp_cov_rounded_int;
}
        
sub get_all_exp_cov
{
    my $kmer_list = get_kmer_list;
    #my @exp_cov_list = ();
    my ($total_cov, $avg_readlen) = calc_common_stats;
    for my $kmer (@$kmer_list) {
        my $exp_cov = calc_exp_cov ($total_cov, $avg_readlen, $kmer);
        print $kmer . "\t" . $exp_cov . "\n";
        #push (@exp_cov_list, $exp_cov);
    }
    #print join(",",@exp_cov_list);
}

get_all_exp_cov;
