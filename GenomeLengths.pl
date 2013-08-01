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
            genome_length_infile input_data/GenomeLengthEst.tab
            genome_length_outfile output_files/GenomeLengthEstOut.tab
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
                --genome_length_infile <filename>
                --genome_length_outfile <filename>
                --verbose
                ";
    }
}


sub gather_opts
{
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'genome_length_infile=s',
            'genome_length_outfile=s',
            'verbose',
            );
    set_default_opts;
    check_opts;
}

sub get_check_record
{
    my $ref = shift;
    my $kref = shift;
    for my $key (@$kref) {
        if (defined ($ref->{$key})) {
            $ref = $ref->{$key};
        } else {
            return '';
        }
    }
    return $ref;
}

sub get_genomic_samples
{
    my $sample_list = [];
    for my $sample (keys %$records) {
        if ($records->{$sample}->{bio_type} =~ /DNA/) {
            push (@$sample_list, $sample);
        }
    }
    return $sample_list;
}

sub parse_genome_lengths
{
    my $fname = ($options->{genome_length_infile} ? $options->{genome_length_infile} : '');
    if ($fname) {
        open (FGLEN, '<', $fname) or die "Error: couldn't open genome lengths input file $fname\n";
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

sub write_table
{
    my $fname = ($options->{genome_length_outfile} ? $options->{genome_length_outfile} : '');
    if ($fname) {
        open (FGOUT, '>', $fname) or die "Error: couldn't open genome lengths output file: $fname\n";
        print FGOUT join ("\t", @glen_col_headers) . "\n";
        my $genomic_samples = get_genomic_samples;
        my @glen_related = grep (/^RG/, @glen_col_headers);
        for my $sample (@$genomic_samples) {
            my $rec = $records->{$sample};
            my $species = get_check_record($rec, ["species"]);
            my $genotype  = '';
            my $biomaterial = '';
            my $biomaterial_type = get_check_record($rec, ["bio_type"]);
            my $kingdom = '';
            my $res_expert = '';
            my @print_line = ($species, $genotype, $biomaterial, $biomaterial_type, 
                    $sample, $kingdom, $res_expert);
            for my $key (@glen_related) {
                my $val = get_check_record($rec, ["related_genome_length", $key]);
                push (@print_line, $val);
            }
            
            print FGOUT join ("\t", @print_line) . "\n";
        }
        close (FGOUT);
    }
}

#my @glen_col_headers = qw(IL_Species IL_Genotype IL_Biomaterial IL_Biomaterial_type 
#    IL_Sample_ID IL_Kingdom IL_Resident_Expert RG_Species RG_Strain RG_Est_Genome_Length
#     RG_Source RG_Citation RG_URL RG_Notes);

gather_opts;
$records = LoadFile($options->{yaml_in});
parse_genome_lengths;
#write_table;
DumpFile($options->{yaml_out}, $records);

