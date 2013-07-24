#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my $records = {};
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
my @sample_list = ();
my $tr;
my $trdata;

# output files created by velvetg
my $vg_outfiles = [qw(Graph2 LastGraph PreGraph stats.txt)];

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/10_velvetg_qsub.yml
        yaml_out yaml_files/11_velvet_stats.yml
        verbose 0
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    $options->{qsub_opts} = $options->{qsub_opts} . ""; # . -N velvet_hg ";
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
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
        );
    set_default_opts;
    check_opts;
}

sub get_check_record
{
    my $ref = shift;
    my $kref = shift;
    for my $key (@$kref) {
        if (defined ($ref->{$key})) {
            $ref = $ref->{$key};
        } else {
            return '';
        }
    }
    return $ref;
}

sub set_check_record
{
    my $ref = shift;
    my $kref = shift;
    my $last_key = shift;
    my $value = shift;
    for my $key (@$kref) {
        if (defined ($ref->{$key})) {
            $ref = $ref->{$key};
        } else {
            $ref->{$key} = {};
            $ref = $ref->{$key};
        }
    }
    $ref->{$last_key} = $value;
}

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $tr = shift;
    my $kmin = get_check_record($rec, ["velvet", $tr, "min_kmer"]);
    my $kmax = get_check_record($rec, ["velvet", $tr, "max_kmer"]);
    my $krange = [];
    if ($kmin =~ /^\d+$/ and $kmax =~ /^\d+$/) {
        for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
            push (@$krange, $kmer);
        }
    }
    return $krange;
}

gather_opts;
$records = LoadFile($options->{yaml_in});
for my $sample (keys %$records) {
    print $sample . "\n";
    my $rec = $records->{$sample};
    for my $tr (qw(trim raw)) {
        my $krange = get_kmer_range($rec, $tr);
        if (scalar @$krange > 0) {
            my %max;
            $max{kmer} = $krange->[0];
            $max{n50} = 0;
            for my $kmer (@$krange) {
                #print "!!\n";
                my $kmer_dir = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "kmer_dir"]);
                #print "got kmer dir " . $kmer_dir . "\n";
                my $log_file = $kmer_dir . "/Log";
                my $contig_file = $kmer_dir . "/contigs.fa";
                if (-e $log_file and -e $contig_file and -s $contig_file > 0) {
                    open (FIN, '<', $log_file) or die "Error: couldn't open file ${log_file}\n";
                    my $n50 = '';
                    while (<FIN>) { 
                        if (/n50\s+of\s+(\d+)/i) {
                            $n50 = $1;
                            print "kmer $kmer - found n50 " . $n50 . "\n";
                        }
                    }
                    if ($n50) {
                        set_check_record($rec, ["velvet", $tr, "kmer", $kmer], "N50", $n50);
                        if ($n50 > $max{n50}) {
                            $max{n50} = $n50;
                            $max{kmer} = $kmer;
                        }
                    }
                }
            }
            if ($max{n50} > 0) {
                print "For sample $sample $tr got max n50 kmer " . $max{kmer} . " value " . $max{n50} . "\n";
                set_check_record($rec, ["velvet", $tr], "max_n50_kmer", $max{kmer});
                set_check_record($rec, ["velvet", $tr], "max_n50_value", $max{n50});
            }
        }
    }
}
DumpFile($options->{yaml_out}, $records);








