#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS;
use POSIX qw(strftime);
use Assembly::Utils;

my $options = {};

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/09_velvet_advisor.yml
        yaml_out yaml_files/10_velvetk.yml
        release_table input_data/Release.tab
        verbose 0
    );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file> -r <release table>
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
        'release_table|r=s',
        'verbose',
        );
}

sub get_software_versions
{
    my $versions = {};
    $versions->{fastqx} = `fastqc --version`;
    $versions->{fastx_trimmer} = `fastx_trimmer -h | grep FASTX`;
    $versions->{velvetk} = "2012\n";
    $versions->{velveth} = `velveth_31 | grep Version`;
    $versions->{velvetg} = `velvetg_31 | grep Version`;
    return $versions;
}

# Pull records such as input raw reads file location
# Any trimming work done
# Sequencing metadata (and read lengths, estimated genome size, etc.)
sub get_sample_records
{
    my $rec = shift;
    my $trimraw = shift;
    
    my $sr = {};
    $sr->{sample_id} = Assembly::Utils::get_check_record($rec, ["sample
}

sub create_release
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $read_types = shift;
    my $trimraw = shift;
    my $best_kmer = shift;
    
    # What we'll be handing back in our release yaml
    my $release_records = {};
    
    # Get the release part
    my $release = {};
    $release->{species} = $rec{Species};
    $release->{strain} = $rec{Strain};
    $release->{version} = "R01";
    $release->{date} = strftime "%m/%d/%Y", localtime;
    $release->{read_types} = $rec{Read_Types};
    $release_records->{release} = $release;
    
    # For each of PE, MP3 MP8, get location of base reads file
    # and the metadata for the file
    # and any trim commands used to generate the file.
    for my $sample_type (split(/,/, $rec{Read_Types})) {
        my $rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type]);
        $release_records->{$sample_type} = get_sample_records($rec, $trimraw);
    }
    
    # Now get the sequence of velvet commands used to build the file
    my $velvet_ref = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet"]);
    my $release_records->{pipeline} = get_pipeline_info($velvet_ref);
    return $release_records;
}

sub parse_release_table
{
    my $fname = ($options->{release_table} ? $options->{release_table} : '');
    unless ($fname) { return; }
    my @headers = qw (Species Strain Read_Types Trim_Raw Best_Kmer Best_N50);
    open (RTAB, '<', $fname) or die "Error: couldn't open release table file $fname\n";
    <RTAB>; # skip header
    while (my $line = <RTAB>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        my $num_headers = scalar @headers;
        my $num_fields = scalar @fields;
        my %rec;
        for (my $i=0; $i<$num_headers; $i++) { 
            my $key = $headers[$i];
            my $val = ($fields[$i] ? $fields[$i] : '');
            $rec{$key} = $val; 
        }
        my $species = Assembly::Utils::format_species_key($rec{Species});
        my $strain = Assembly::Utils::format_strain_key($rec{Strain});
        my $trimraw = $rec{Trim_Raw};
        create_release($records, $species, $strain, $rec{Read_Types}, $rec{Trim_Raw}, $rec{Best_Kmer});
    }
}



# Open release spreadsheet
# Given hte species/strain/trim/raw/kmer
# Build the release info

# Q: what should the output file look like?

gather_opts;
my $records = LoadFile($options->{yaml_in});
my $versions = get_software_versions;










