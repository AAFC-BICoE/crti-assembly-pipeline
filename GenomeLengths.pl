#!/usr/bin/env perl
use strict;
use warnings;
use YAML::XS qw(LoadFile DumpFile);
use File::Path;
use Getopt::Long;

# From the records, pull just the genomic data
# Read in genome lengths from file
# Write out all the velveth/g commands per sample per kmer

# Ideally: have separate folders within velvet that include
# e.g. the different genome lengths, QC trims etc. that we try?

my $options = {};
my $records = {};
my @genomic_samples = ();
my @glen_col_headers = qw(IL_Species IL_Genotype IL_Biomaterial IL_Biomaterial_type 
    IL_Sample_ID IL_Kingdom IL_Resident_Expert RG_Species RG_Strain RG_Est_Genome_Length
     RG_Source RG_Citation RG_URL RG_Notes);

sub set_default_opts
{
    my %defaults = qw(
            yaml_in yaml_files/05_readinfo.yml
            yaml_out yaml_files/06_genome_lengths.yml
            genome_length_file input_data/GenomeLengthEst.tab
            verbose 1
            );
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file>
            Optional:
                --genome_length_file <filename>
                --verbose
                ";
    }
}


sub gather_opts
{
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'genome_length_file|g=s',,
            'verbose',
            );
    set_default_opts;
    check_opts;
}

sub parse_genome_lengths
{
    my $fname = ($options->{genome_length_file} ? $options->{genome_length_file} : '');
    if ($fname) {
        open (FGLEN, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FGLEN>; <FGLEN>; # Skip first two col header lines.
        while (my $line = <FGLEN>) {
            chomp $line;
            my @fields = split(/\t/, $line);
            my %fh = map { $glen_col_headers[$_] => ($fields[$_] ? $fields[$_] : '') } (0..$#glen_col_headers);
            my $sample = $fh{IL_Sample_ID};
            $sample =~ s/\s+//g;
            my $rec = ($records->{$sample} ? $records->{$sample} : '');
            if ($rec) {
                $rec->{related_genome_length} = {};
                my $rgl_ref = $rec->{related_genome_length};
                for my $ch (@glen_col_headers) {
                    if ($ch =~ /^RG/) {
                        $rgl_ref->{$ch} = $fh{$ch};
                    }
                }
            } else {
                print "Warning: in file $fname, sample $sample was not found in input yaml file.\n";
            }
        }
        close (FGLEN);
    }
}

gather_opts;
$records = LoadFile($options->{yaml_in});
parse_genome_lengths;
DumpFile($options->{yaml_out}, $records);

