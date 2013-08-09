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
        table_out output_files/AssemblySetupOut.tab
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
                --table_out <filename>
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'table_out=s',
        'verbose',
        );
    set_default_opts;
    check_opts;
}

sub get_sample_type
{
    my $rec = shift;
    my $biomaterial = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Biomaterial"]);
    my $species_long = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Organism"]);
    my $notes = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Notes"]);
    my $type = "PE";
    if ($biomaterial =~ /\b3\s*kb\b/i) {
        $type = "MP3";
    } elsif ($biomaterial =~ /\b8\s*kb\b/i) {
        $type = "MP8";
    } elsif ($biomaterial =~ /\bMP\b/) {
        $type = "MP";
    } elsif ($species_long =~ /resequencing/i or $notes =~ /resequencing/i) {
        $type = "PER";
    }
    return $type;
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
        my $sample_type = get_sample_type($rec);
        my $type_filled = Assembly::Utils::get_check_record($super_records, [$species_key, $bio_type, $strain_key, $sample_type]);
        if ($type_filled) {
            my $prev_sample = Assembly::Utils::get_check_record($super_records, [$species_key, $bio_type, $strain_key, $sample_type, "sample"]);
            #die "Error: Trying to assign sample $sample as species: $species, bio_type: $bio_type, " .
                #"strain: $strain_key, sample type: $sample_type\nbut another sample already exists there: $prev_sample\n";
        } else {
            Assembly::Utils::set_check_record ($super_records, [$species_key, $bio_type, $strain_key], $sample_type, $rec);
        }
    }
    return $super_records;
}

sub write_output_file
{
    my $super_records = shift;
    my @headings_start = qw(Species Strain);
    my @type_list = qw(PE PER MP MP3 MP8);
    my $fname = ($options->{table_out} ? $options->{table_out} : '');
    if ($fname) {
        open (FOUT, '>', $fname) or die "Error: couldn't open file $fname\n";
        print FOUT join("\t", @headings_start) . "\t" . join("\t", @type_list) . "\n";
        for my $species (keys %$super_records) {
            my $ref = $super_records->{$species}->{DNA};
            for my $strain (keys %$ref) {
                print FOUT $species . "\t" . $strain . "\t";
                for my $sample_type (@type_list) {
                    my $sample_val = Assembly::Utils::get_check_record($super_records, [$species, "DNA", $strain, $sample_type, "sample"]);
                    print FOUT $sample_val . "\t";
                }
                print FOUT "\n";
            }
        }
        close (FOUT);
    }
}
                    

gather_opts;
my $records = LoadFile($options->{yaml_in});
my $super_records = alter_super_records($records);
write_output_file($super_records);
DumpFile($options->{yaml_out}, $super_records);



