#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use Assembly::Utils;

# Keep all files for best kmer (trim and raw)
# Remove all but contigs.fa, Log, and Roadmaps files for others.
# Do not touch folders with files less than one day old - may still be processing
# Future: compression of Roadmaps/contigs for non-best.
# Future: removal of all but Log file for poor (e.g. <2000) N50's
# Future: support for inclusion of other 'bests' e.g. max contig length.

# Note: we assume all of the stats are in the input yaml file and 
# do not explicitly parse the VelvetStats.tab file.

my $options = {};
my @del_files = qw(CnyUnifiedSeq.names Roadmaps stats.txt LastGraph);
# push (@del_files, "CnyUnifiesSeq");
# push (@del_files, "Graph2")
# push (@del_files, "PreGraph");

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/12_velvet_stats.yml
        verbose 1
        run 0
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
                --run
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'verbose',
        'run',
        );
    set_default_opts;
    check_opts;
}

sub print_verbose
{
    if ($options->{verbose}) {
        print (@_);
    }
}

sub remove_files
{
    my $rec = shift;
    my $kmer = shift;
    my $kdir = Assembly::Utils::get_check_record($rec, ["kmer", $kmer, "kmer_dir"]);
    if ($kdir and -e $kdir) {
        for my $fname (@del_files) {
            my $fpath = $kdir . "/" . $fname;
            if (-e $fpath and (-M $fpath >= 1)) {
                print_verbose("Removing file $fpath\n");
                if ($options->{run}) {
                    system ("rm $fpath");
                }
            } else {
                print_verbose("No such file or file modified within the past 24 hours: $fpath\n");
            }
        }
    }
}

sub remove_range
{
    my $rec = shift;
    my $kmin = Assembly::Utils::get_check_record($rec, ["min_kmer"]);
    my $kmax = Assembly::Utils::get_check_record($rec, ["max_kmer"]);
    my $best_kmer = Assembly::Utils::get_check_record($rec, ["max_n50_kmer"]);
    if ($best_kmer and $kmin and $kmax) {
        for (my $k = $kmin; $k <= $kmax; $k += 2) {
            if ($k != $best_kmer) {
                remove_files($rec, $k);
            }
        }
    } else {
        print_verbose "Could not determine best/min/max kmers for this record";
    }
} 

sub remove_all
{
    my $records = shift;
    for my $species (keys %$records) {
        for my $strain (keys %{$records->{$species}->{DNA}}) {
            for my $trimraw (qw(trim raw)) {
                print "species $species strain $strain trimraw $trimraw\n";
                my $rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw]);
                if ($rec) {
                    remove_range($rec);
                } else {
                    print_verbose "Could not get YAML record for species $species strain $strain trim/raw: $trimraw\n";
                }
            }
        }
    }
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
remove_all($records);


















