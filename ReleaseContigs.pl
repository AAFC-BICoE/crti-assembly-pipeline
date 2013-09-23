#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $options = {};

sub set_default_opts
{
    my %defaults = qw(
        contigs_in Ac_LMG_21996_R01_contigs.fa
        contigs_out Ac_contigs_out.fa
        run 0
        verbose 1
    );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{contigs_in} and $options->{contigs_out}) {
        die "Usage: $0 -i <input contigs file> -o <output contigs file>
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
        'contigs_in|i=s',
        'contigs_out|o=s',
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

sub rename_contigs
{
    my $contigs_in = ($options->{contigs_in} ? $options->{contigs_in} : '');
    my $contigs_out = ($options->{contigs_out} ? $options->{contigs_out} : '');
    if ($contigs_in and $contigs_out) {
        my $prefix = '';
        if ($contigs_in =~ /([^\/]+)_contigs.fa/) {
            $prefix = $1;
        }
        open (CIN, '<', $contigs_in) or die "Error: couldn't open file $contigs_in\n";
        open (COUT, '>', $contigs_out) or die "Error: couldn't open file $contigs_out\n";
        my $node_count = 1;
        while (my $line = <CIN>) {
            if ($line =~ /^>(\S+)/) {
                my $old_id = $1;
                my $new_id = $prefix . "_" . $node_count;
                print COUT ">" . $new_id . " " . $old_id . "\n";
                $node_count++;
            } else {
                print COUT $line;
            }
        }
        close (CIN);
        close (COUT);
    }
}

gather_opts;
rename_contigs;











