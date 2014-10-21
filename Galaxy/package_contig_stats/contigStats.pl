#!/usr/bin/env perl


use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $options = {};
GetOptions($options,
        'contigs_filename|c=s'
        );

sub calc_contig_stats
{
    my $contigs_filename = shift;
    my $contigs = Bio::SeqIO->new(-file => $contigs_filename, -format => "Fasta");

    my @seqlen = ();
    while (my $seq = $contigs->next_seq()) {
        push(@seqlen, length($seq->seq));
    }

    @seqlen = sort { $b <=> $a } @seqlen;
    my $num_contigs = scalar @seqlen;

    my $max_contig_len = $seqlen[0];
    my $min_contig_len = $seqlen[$num_contigs-1];

    my $mid_idx = int ($num_contigs - 0.5) / 2; # 0-base index
    my $median_contig_len = '';
    if ($num_contigs % 2 == 0) {
        $median_contig_len = ($seqlen[$mid_idx] + $seqlen[$mid_idx+1])/2;
    } else {
        $median_contig_len = $seqlen [$mid_idx];
    }
    
    return [$num_contigs, $min_contig_len, $median_contig_len, $max_contig_len];
}

if ($options->{contigs_filename}) {
    my $stats = calc_contig_stats($options->{contigs_filename});
    my ($count, $min, $median, $max) = @$stats;
    print "Num contigs: " . $count . "\n";
    print "Min contig length: " . $min . "\n";
    print "Median contig length: " . $median . "\n";
    print "Max contig length: " . $max . "\n";
}

