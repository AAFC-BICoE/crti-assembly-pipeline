#!/usr/bin/env perl
# Assembly::Velvet.pm
# Author: Jeff Cullis
# Date: August 11, 2013

package Assembly::Velvet;

use strict;
use warnings;
use Assembly::Utils;

my $velvet_bin_dir = "/opt/bio/velvet";
my @kbins = (0, 31, 63, 127, 145); 

sub get_assembly_outdir
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $species_dir = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "PE", "species_dir"];
    my $assembly_outdir = $species_dir . "DNA/assemblies/" . $strain . "/velvet/";
    unless (-e $assembly_outdir) {
        mkpath $assembly_outdir;
    }
    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "assembly_outdir", $assembly_outdir);
    return $assembly_outdir;
}

sub create_kmer_range
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $use_velvetk = 1;
    my $vk_radius = 6;
    my $min_kmer = 21;
    my $max_kmer = 101;
    my $kmer_range = [];
    if ($use_velvetk and $vk_radius) {
        my $velvetk_best = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_best_kmer"]);
        my $rad = $vk_radius * 2;
        if ($velvetk_best =~ /^\s*\d+\s*$/ and $rad =~ /^\s*\d+\s*$/) {
            if ($velvetk_best - $rad < 21) {
                $min_kmer = 21;
            } else {
                $min_kmer = $velvetk_best - $rad;
            }
            if ($velvetk_best + $rad > 101) {
                $max_kmer = 101;
            } else {
                $max_kmer = $velvetk_best + $rad;
            }
        } else {
            # If no velvetk value found, don't create any commands.
            $min_kmer = 1;
            $max_kmer = 0;
        }
    }
    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "min_kmer", $min_kmer);
    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "max_kmer", $max_kmer);
    for (my $i = $min_kmer; $i <= $max_kmer; $i = $i+2) {
        push (@$kmer_range, $i);
    }
    return $kmer_range;
}

sub get_kmer_range
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $kmin = get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "min_kmer"]);
    my $kmax = get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "max_kmer"]);
    my $krange = [];
    if ($kmin =~ /^\d+$/ and $kmax =~ /^\d+$/) {
        for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
            push (@$krange, $kmer);
        }
    }
    return $krange;
}

# Given a kmer, determine which velvet binary (defined in @kbins, top) should be used to run it.
sub get_kmer_bin
{
    my $kmer = shift;
    my $bin = '';
    if ($kmer < $kbins[0]) { 
        print "Error: kmer value $kmer must be a positive integer!\n";
    } elsif ($kmer > $kbins[$#kbins]) {
        print "kmer value $kmer not supported. Recompile velvet for higher kmer.\n";
    } else {
        for (my $i = 0; $i < $#kbins; $i++) {
            if ($kbins[$i] < $kmer and $kmer <= $kbins[$i+1]) {
                $bin = $kbins[$i+1];
            }
        }
    }
    return $bin;
}

sub get_coverage_vars
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $trdata = $trimraw . "data";
    
    
    my $var = {};
    $var->{R1_nreads} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R1", $trdata, "num_reads"]);
    $var->{R1_readlen} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R1", $trdata, "read_length"]);
    $var->{R2_nreads} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R2", $trdata, "num_reads"]);
    $var->{R2_readlen} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R2", $trdata, "read_length"]);
    $var->{genome_length} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "related_genome_length", "RG_Est_Genome_Length"]);
    my $pass = 1;
    for my $key (keys %$var) {
        if ($var->{$key} !~ /^\s*\d+\s*$/) {
            print "Got bad var " . $var->{$key} . " at key " . $key . "\n";
            $pass = 0;
        }
    }
    if ($pass) {
        $var->{total_coverage} = ($var->{R1_nreads} * $var->{R1_readlen} + $var->{R2_nreads} * $var->{R2_readlen}) / $var->{genome_length};
        $var->{avg_readlen} = ($var->{R1_readlen} + $var->{R2_readlen}) / 2;
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "total_coverage", $var->{total_coverage});
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "average_read_length", $var->{avg_readlen});
        return $var;
    } else {
        return '';
    }
}

sub calc_exp_cov
{
    my $var = shift;
    my $kmer = shift;
    my $exp_cov_float = $var->{total_coverage} * ($var->{avg_readlen} - $kmer + 1) / $var->{avg_readlen};
    my $exp_cov_rounded_int = int($exp_cov_float + 0.5);
    return $exp_cov_rounded_int;
} 

return 1;  
