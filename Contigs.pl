#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use Bio::SeqIO;

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/12_velvet_stats.yml
        verbose 0
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{yaml_in}) {
        die "Usage: $0 -i <input yaml file>
            Optional:
                --verbose
                --stats_outfile <filename>
                ";
    }
}

sub gather_opts
{
    GetOptions($options,
        'yaml_in|i=s',
        'verbose',
        );
    set_default_opts;
    check_opts;
}

sub get_kmer_range
{
    my $rec = shift;
    my $tr = shift;
    my $kmin = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "min_kmer"]);
    my $kmax = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "max_kmer"]);
    my $krange = [];
    if ($kmin =~ /^\d+$/ and $kmax =~ /^\d+$/) {
        for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
            push (@$krange, $kmer);
        }
    }
    return $krange;
}

sub get_contigs_stats
{
    my $kdir = shift;
    my $contigs_path = $kdir . "/contigs.fa";
    my $seq_in = Bio::SeqIO->new(-file => $contigs_path, -format => 'FASTA');
    my @contig_counts = ();
    while (my $seq = $seq_in->next_seq) {
        push (@contig_counts, length($seq->seq));
    }
    # sort descending
    my @sorted_counts = sort {$b <=> $a} @contig_counts;
    return \@sorted_counts;
}   

sub write_contig_counts
{
    my $count_arr = shift;
    my $fname = ($options->{count_file} ? $options->{count_file} : '');
    if ($fname) {
        open (FCOUNT, '>', $fname) or die "Error: couldn't open ouptut counts file $fname\n";
        print FCOUNT join ("\n", @$count_arr) . "\n";
        close (FCOUNT);
    }
}

sub parse_contigs
{
    my $records = shift;
    for my $species (keys %$records) {
        for my $strain (keys %${records->{$species}->{DNA}}) {
            for my $trimraw (qw(trim raw)) {
                my $kmer_range = get_kmer_range($rec, $trimraw);
                for my $kmer (@$kmer_range) {
                    my $kdir = Assembly::Utils::get_check_records($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer, "kmer_dir"]);
                    get_contigs_stats($kdir);
                }
            }
        }
    }
}

gather_opts;
$records = LoadFile($options->{yaml_in});



# If no counts file, create
if ($options->{input_fasta}

if ($options->{input_fasta} and $options->{count_file}) {
    count_contigs;
}

