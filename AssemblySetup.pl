#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use Assembly::Utils;

my $options = {};

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


sub alter_super_records
{
    my $records = shift;
    my $super_records = {};
    for my $sample (keys %$records) {
        my $rec = $records->{$sample};
        my $bio_type = Assembly::Utils::get_check_record($rec, ["bio_type"]);
        my $species = Assembly::Utils::get_check_record($rec, ["species"]);
        my $species_key = Assembly::Utils::format_species_key($species);
        my $strain = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Genotype"]);
        my $strain_key = Assembly::Utils::format_strain_key($strain);
        Assembly::Utils::set_check_record ($super_records, [$species_key, $bio_type, $strain_key, "samples"], $sample, $rec);
    }
    return $super_records;
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
my $super_records = alter_super_records($records);
DumpFile($options->{yaml_out}, $super_records);



