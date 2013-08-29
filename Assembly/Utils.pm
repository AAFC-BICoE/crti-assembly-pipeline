#!/usr/bin/env perl
# Assembly::Utils.pm
# Author: Jeff Cullis
# Date: August 01, 2013

package Assembly::Utils;

use strict;
use warnings;

sub get_check_record
{
    my $ref = shift;
    my $kref = shift;
    for my $key (@$kref) {
        if ($ref and defined ($ref->{$key})) {
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

sub format_species_key
{
	my $species = shift;
	my $flag = 0;
	$species =~ s/^\s+|\s+$//g;
	$species =~ s/ /_/g;
	if ($species =~ /^([A-Z][a-z]+_[a-z]+)/) {
		$species = $1;
	} else {
		die "Error: malformed species name: $species. Aborting.\n";
	}
	return $species;
}

sub format_strain_key
{
    my $strain = shift;
    $strain =~ s/^\s+|\s+$//g;
    $strain =~ s/[^0-9A-Za-z \_\-]//g;
    $strain =~ s/\bstrain\b\s*//i;
    $strain =~ s/\s+/_/g;
    unless ($strain =~ /\S/) {
        $strain = "unknown";
    }
    return $strain;
}

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $trimraw = shift;
    my $kmin = get_check_record($rec, ["velvet", $trimraw, "min_kmer"]);
    my $kmax = get_check_record($rec, ["velvet", $trimraw, "max_kmer"]);
    my $krange = [];
    if ($kmin =~ /^\d+$/ and $kmax =~ /^\d+$/) {
        for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
            push (@$krange, $kmer);
        }
    }
    return $krange;
}

sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

return 1;