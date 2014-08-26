#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use YAML qw(DumpFile LoadFile);
use Bio::SeqIO;

my ($quast_report, $contigs_file, $velvet_log, $yaml_out) = ('','','', '');
GetOptions(
    'quast_report|q=s' => \$quast_report,
    'contigs_file|c=s' => \$contigs_file,
    'velvet_log|v=s' => \$velvet_log,
    'yaml_out|y=s' => \$yaml_out
    );

my $genome_stats = {};

sub fill_yaml
{
    $genome_stats->{assembly_sample_types} = "";
    $genome_stats->{estimated_genome_length} = "";
    $genome_stats->{kingdom} = "";
    $genome_stats->{reads_used} = "";
    $genome_stats->{release_kmer} = "";
    $genome_stats->{species_dir} = "";
}

sub parse_quast_stats
{
    open (QUAST, '<', $quast_report) or die "Error: couldn't open file $quast_report\n";
    while (my $line = <QUAST>) {
        if ($line =~ /N50\s+([0-9]+)/) {
            $genome_stats->{N50} = $1;
        }
        if ($line =~ /Total length\s+([0-9]+)/) {
            $genome_stats->{total_length} = $1;
        }
        if ($line =~ /Largest contig\s+([0-9]+)/) {
            $genome_stats->{max_contig} = $1;
        }
        if ($line =~ /# contigs.*\b0 bp\)\s+([0-9]+)/) {
            $genome_stats->{num_contigs} = $1;
        }
    }
    close (QUAST);
}

sub parse_velvet_stats
{
    open (VLOG, '<', $velvet_stats) or die "Error: couldn't open file $velvet_stats.\n";
    while (my $line = <VLOG>) {
        if ($line =~ /Median coverage depth = ([0-9\.]+)/) {
            $genome_stats->{velvet_median_coverage} = $1;
        }
        if ($line =~ /using ([0-9]+\/[0-9]+) reads/) {
            $genome_stats->{reads_used} = $1;
        }
    }
    close (VLOG);
}

sub calc_contig_stats
{
    my $contigs = Bio::SeqIO->new(-file => $contigs_file, -format => "Fasta");

    my @seqlen = ();
    while (my $seq = $contigs->next_seq()) {
        push(@seqlen, length($seq->seq));
    }

    @seqlen = sort { $b <=> $a } @seqlen;
    my $num_contigs = scalar @seqlen;

    # my $max_contig_len = $seqlen[0];
    # my $min_contig_len = $seqlen[$num_contigs-1];
    $genome_stats->{min_contig_len} = $seqlen[$num_contigs-1];

    my $mid_idx = int ($num_contigs - 0.5) / 2; # 0-base index
    my $median_contig_len = '';
    if ($num_contigs % 2 == 0) {
        $median_contig_len = ($seqlen[$mid_idx] + $seqlen[$mid_idx+1])/2;
    } else {
        $median_contig_len = $seqlen [$mid_idx];
    }
    $genome_stats->{median_contig_len} = $median_contig_len;
    # return [$num_contigs, $min_contig_len, $median_contig_len, $max_contig_len];
}

fill_yaml;
parse_quast_stats;
if ($velvet_log) {
    parse_velvet_stats;
}
calc_contig_stats;
DumpFile($yaml_out, $genome_stats);
