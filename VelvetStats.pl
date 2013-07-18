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
        yaml_in yaml_files/08_velvetg_qsub.yml
        yaml_out yaml_files/09_velvet_n50.yml
        trim 1
        verbose 0
        submit 0
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    if ($options->{raw}) {
        $tr = "raw";
        $trdata = "rawdata";
    } else {
        $tr = "trim";
        $trdata = "trimdata";
    }
    $options->{qsub_opts} = $options->{qsub_opts} . ""; # . -N velvet_hg ";
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --sample_list <ID1,ID2,...,IDX (no spaces)>
                --qsub_array_script <filename>
                --memory
                --trim
                --raw
                --submit
                --submit_max <max # qsubs to perform>
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'qsub_script=s',
        'verbose',
        'sample_list=s',
        'memory_requirements|m=s',
        'trim',
        'raw',
        'submit',
        'submit_max=s',
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