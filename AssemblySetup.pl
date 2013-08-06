#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS;
use Assembly::Util;

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/05_readinfo.yml
        yaml_out yaml_files/06_assembly_setup.yml
        verbose 0
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    #$options->{qsub_opts} = $options->{qsub_opts} . " -p -500 "; # . -N velvet_hg ";
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

sub split_sample_list
{
    my $records = shift;
    my $genomic_samples = [];
    my $rnaseq_samples = [];
    for my $sample (%$records) {
        my $rec = $records->{sample};
        if ($rec->{bio_type} eq "DNA") {
            push (@$genomic_samples, $sample);
        } elsif ($rec->{bio_type} eq "RNA") {
            push (@$rnaseq_samples, $sample);
        }
    }
    return ($genomic_samples, $rnaseq_samples);
}

sub genomic_super_rec
{
    my $records = shift;
    my $genomic_samples = shift;
    my $genomic_super = {};
    for my $sample (@$genomic_samples) {
        my $rec = $records->{$sample};
        my $species = $rec->{species};
        my $strain = get_check_record($rec, ["sequence_metadata", "Genotype"]);
        my $ga_key = get_genomic_key($species, $strain);
        set_check_record($genomic_super, [$ga_key, "samples"], $sample, $rec);
        #$genomic_super{$ga_key}{samples}{$sample} = $rec;
    }
    return $genomic_super;
}

sub get_assembly_ids
{
    my $species = $rec->{
    

my $records = LoadFile($options->{yaml_in});
my ($genomic_samples, $rnaseq_samples) = split_sample_list($records);
# Add Transcriptome_Assemblies and Genome_Assemblies super-records
# Genome_Assemblies: Make new keys species_sample with approp. sample elements below.
my $assemblies = setup_assemblies($records, $sample_list);


