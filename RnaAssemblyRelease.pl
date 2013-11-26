#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Assembly::Utils;
use YAML::XS qw(LoadFile DumpFile);

# Open and parse the table as in prev script...????

# Open our metadata files and add the release section, then write to 
# new output file.

my $options = {};
my @colheaders = ("Species", "Strain", "Sample_ID", "Reference_Strain", "Reference_Metafile", 
        "Output_Release_Dir", "Output_Release_Prefix");

sub set_default_opts
{
    my %defaults = qw(
        rna_assembly_table input_data/RnaAssemblyTable2.tab
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
    }
    return $table_recs;
}

sub get_bowtie_stanza
{
    my $yrec = shift;
    #my $bowtie_cmd = Assembly::Utils::get_check_record($
}

sub get_tophat_stanza
{

}

sub get_cufflinks_stanza
{

}

sub update_release
{
    my ($input_release_filename, $output_release_filename, $species, $sample, $strain, $yrec) = @_;
    my $rrin = LoadFile ($input_release_filename);
    print "Got rrin\n";
    my $pipeline_ref = $rrin->[4]->{pipeline};
    my $bowtie_stanza = get_bowtie_stanza ($yrec);
    my $tophat_stanza = get_tophat_stanza ($yrec);
    my $cufflinks_stanza = get_cufflinks_stanza ($yrec);
    push (@$pipeline_ref, $bowtie_stanza);
    push (@$pipeline_ref, $tophat_stanza);
    push (@$pipeline_ref, $cufflinks_stanza);
    #DumpFile ($output_release_filename, $rrin);
}
    

sub create_release
{
    # Open the yaml metadata
    # Pull all of the commands into the pipeline
    # Add a section for stats about the transcriptome???
    my $yaml_records = shift;
    my $table_records = shift;

    for my $sample (keys %$table_records) {
    #for my $sample ("S001274") {
        my $species = $table_records->{$sample}->{"Species"};
        my $strain = $table_records->{$sample}->{"Strain"};
        my $input_release_filename = $table_records->{$sample}->{"Reference_Metafile"};
        my $output_release_filename = $table_records->{$sample}->{"Output_Release_Dir"} . 
            $table_records->{$sample}->{"Output_Release_Prefix"} . "_metadata.yml";
        my $yrec = Assembly::Utils::get_check_record ($yaml_records, [$species, "RNA", $strain, $sample, "RNA_assembly"]);
        update_release ($input_release_filename, $output_release_filename, $species, $sample, $strain, $yrec);
    }
}

gather_options;
my $yaml_records = LoadFile ($options->{yaml_in});
my $table_records = parse_assembly_table;
create_release ($yaml_records, $table_records);
























