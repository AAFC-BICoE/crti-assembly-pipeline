#!/usr/bin/env perl
# Assembly::Velvet.pm
# Author: Jeff Cullis
# Date: August 11, 2013

package Assembly::Velvet;

use strict;
use warnings;
use Assembly::Utils;
use File::Path;

my $velvet_bin_dir = "/opt/bio/velvet";
my @kbins = (0, 31, 63, 127, 145); 

sub get_assembly_outdir
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $species_dir = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "PE", "species_dir"]);
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
    my $read_total = 0;
    my $rl_sum = 0;
    my $rl_count = 0;
    for my $sample_type (qw(PE PER MP MP3 MP8)) {
        if ($sample_type =~ /MP/) {
            #$trdata = "rev" . $trdata; # don't do this yet - we don't yet get the reads stats for reversed files.
        }
        for my $rval (qw(R1 R2)) {
            my $numreads = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "data_stats", $rval, $trdata, "num_reads"]);
            my $readlen = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "data_stats", $rval, $trdata, "read_length"]);
            if ($numreads =~ /\d+/ and $readlen =~ /\d+/) {
                $read_total += ($numreads * $readlen);
                $rl_sum += $readlen;
                $rl_count++;
            }
        }
    }
    
    my $avg_readlen = 0;
    if ($rl_count) {
        $avg_readlen = $rl_sum/$rl_count;
    }
    my $genome_length = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "related_genome_length", "RG_Est_Genome_Length"]);
    my $total_coverage = 0;
    if ($genome_length) {
        $total_coverage = $read_total / $genome_length;
    }
    return ($avg_readlen, $total_coverage);
}

sub calc_exp_cov
{
    my $total_coverage = shift;
    my $avg_readlen = shift;
    my $kmer = shift;
    my $exp_cov_float = -0.5;
    if ($avg_readlen) {
        $exp_cov_float = $total_coverage * ($avg_readlen - $kmer + 1) / $avg_readlen;
    }
    my $exp_cov_rounded_int = int($exp_cov_float + 0.5);
    return $exp_cov_rounded_int;
} 

return 1;  
