#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS;
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

# G
