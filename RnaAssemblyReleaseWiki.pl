#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS;

my $options = {};
my @colheaders = ("Species", "Strain", "Sample_ID", "Trim_Raw", "Reference_Strain", "Reference_Metafile", 
        "Output_Release_Dir", "Output_Release_Prefix");

sub set_default_opts
{
    my %defaults = qw(
        rna_assembly_table input_data/RnaAssemblyTable.tab
        genome_lengths_table input_data/GenomeLengths.tab
        wiki_table_file output_files/wiki_transcript_table.txt
        specimen_dir ../../processing_test2
        yaml_in yaml_files/15_rna_setup.yml
        yaml_out yaml_files/16_rna_release.yml
        create_release 0
        );
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
}

sub check_options
{
        unless ($options->{rna_assembly_table} and $options->{yaml_in} ) {
                die "Usage: $0 -a <RNA assembly table> -i <yaml input file> -o <yaml output file>
                        Optional:
                                --genome_lengths_table <filename>
                                --wiki_table_file <filename>
                                --testing
                                --verbose
                                --create_release
                                --specimen_dir <path to specimen/ dir>
                                ";
        }
}

sub gather_options
{
        GetOptions($options,
                'rna_assembly_table|a=s',
                'genome_lengths_table|g=s',
                'wiki_table_file|w=s',
                'testing|t',
                'verbose|v',
                'create_release|r',
                'specimen_dir|p=s',
                'yaml_in|i=s',
                'yaml_out|o=s',
                );
        set_default_opts;
        check_options;
}               

sub print_verbose
{
    if ($options->{verbose}) {
        print (@_);
    }
}

# Parse the input table 
sub parse_assembly_table
{
    my $table_recs = {};
    my $fname = ($options->{rna_assembly_table} ? $options->{rna_assembly_table} : '');
    if ($fname and -s $fname) {
        open (FRNA, '<', $fname) or die "Error: couldn't open file $fname\n";
        while (my $line = <FRNA>) {
            chomp $line;
            my @fields = split (/\t/, $line);
            if (scalar (@fields) == scalar (@colheaders)) {
                my $line_record = {};
                foreach my $i (0..$#colheaders) {
                    # Add hash entry key=column header, value=parsed field value
                    $line_record->{$colheaders[$i]} = $fields[$i];
                }
                my $species_key = Assembly::Utils::format_species_key($line_record->{Species});
                $line_record->{Species} = $species_key; # Reset value to standard access format.
                my $strain_key = Assembly::Utils::format_strain_key($line_record->{Strain});
                $line_record->{Strain} = $strain_key; # Reset value to standard access format.
                my $reference_strain_key = Assembly::Utils::format_strain_key($line_record->{Reference_Strain});
                $line_record->{Reference_Strain} = $reference_strain_key;
                my $sample = $line_record->{"Sample_ID"};
                $table_recs->{$sample} = $line_record;
            } else {
                print_verbose ("Error on line $. of rna assembly table file - incorrect number of cols.\n");
            }
        }
        close (FRNA);
    }
    return $table_recs;
}

# Parse the species, kingdom from the genome lengths table
sub parse_genome_lengths_table
{
    my $s2k = {};
    my $fname = ($options->{genome_lengths_table} ? $options->{genome_lengths_table} : '');
    if ($fname and -s $fname) {
        open (FGL, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FGL>; <FGL>; # skip top two (col headers) rows.
        while (my $line = <FGL>) {
            chomp $line;
            my @fields = split (/\t/, $line);
            if (scalar @fields > 1) {
                my $species = Assembly::Utils::format_species_key ($fields[0]);
                my $kingdom = $fields[1];
                $s2k->{$species} = $kingdom;
            }
        }
        close (FGL);
    }
    return $s2k;
}

sub get_key
{
    my ($species, $strain, $release) = @_;
    return join(":", ($species, $strain, $release));
}

sub get_release_recs
{
    my $table_records = shift;   
    # Get each unique combo of reference species, reference strain, release.
    my %release_recs = {};
    
    for my $sample (@$table_records) {
        my $species = $table_records->{$sample}->{"Species"};
        
        my $strain = $table_records->{$sample}->{"Reference_Strain"};
        my $release = $table_records->{$sample}->{"Output_Release_Prefix"};
        my $release_dir = dirname ($table_records->{$sample}->{"Output_Release_Dir"};
        my $key = get_key ($species, $strain, $release);
        unless ($release_recs{$key}) {
            $release_recs->{$key} = [$species, $strain, $release, $release_dir];
        }
    }
    return $release_recs;
}

sub get_release_link
{
    my $release_dir = shift;
    http://biocluster/project_data/CRTI-09S-462RD/specimen/G_rostochiensis/release
    my $link_prefix = "http://biocluster/project_data/CRTI-09S-462RD/specimen";
    $release_dir =~ s/.*\/specimen\///;
    my $link_dir = $link_prefix . "/" . $release_dir;
    my $link = "[[" . $link_dir . "][Release]]";
    return $link;
}
    

sub create_wiki_table
{
    my $yaml_records = shift;
    my $table_records = shift;
    my $spec_to_king = shift;
    my $release_recs = get_release_recs ($table_records);
    my $fname = ($options->{wiki_table_file} ? $options->{wiki_table_file} : '');
    if ($fname) {
        open (FWIKI, '>', $fname) or die "Error: couldn't open file for output: $fname\n";
        foreach my $rkey (keys %$release_recs) {
            my $rec = $release_recs->{$rkey};
            my ($species, $strain, $release, $release_dir) = @$rec;
            my $kingdom = $spec_to_king->{$species};
            my $link = get_release_link ($release_dir);
            $species =~ s/\_/ /g;
            $strain =~ s/\_/ /g;
            my $wiki_line = "|" . join ("|", ($kingdom, $species, $strain, $release, $link)) . "|\n";
            print FWIKI $wiki_line;
        }
        close (FWIKI);
    }
}

my $yaml_records = LoadFile ($options->{yaml_in});
my $table_records = parse_assembly_table;
my $spec_to_king = parse_genome_lengths_table;
create_wiki_table ($yaml_records, $table_records, $spec_to_king);
#DumpFile ($options->{yaml_out});
















