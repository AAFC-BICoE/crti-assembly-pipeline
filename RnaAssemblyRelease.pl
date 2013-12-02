#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Assembly::Utils;
use YAML::XS qw(LoadFile DumpFile);
use Cwd;
use File::Path;

# Open and parse the table as in prev script...????

# Open our metadata files and add the release section, then write to 
# new output file.

my $options = {};
my @colheaders = ("Species", "Strain", "Sample_ID", "Trim_Raw", "Reference_Strain", "Reference_Metafile", 
        "Output_Release_Dir", "Output_Release_Prefix");

sub set_default_opts
{
    my %defaults = qw(
        rna_assembly_table input_data/RnaAssemblyTable.tab
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

sub create_dir
{
    my $dirname = shift;
    unless ($options->{testing}) {
        unless (-e $dirname) {
            print_verbose ("Creating directory $dirname\n");
            mkpath $dirname;
        }
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
    }
    return $table_recs;
}

sub get_sample_stanza
{
    my $strain_rec = shift;
    my $sample_type = shift;
    my $trimraw = shift;
    my $release_prefix = shift;
    my $release_dir = shift;
    my $sample_rec = Assembly::Utils::get_check_record($strain_rec, [$sample_type]);
    my $outrec = {};
    $outrec->{sequencing_metadata} = Assembly::Utils::get_check_record($sample_rec, ["sequencing_metadata"]);
    # $outrec->{species_abbr} = Assembly::Utils::get_check_record($sample_rec, ["species_abbr"]);
    $outrec->{read_data} = {};
    $outrec->{read_data}->{sample_dir} = Assembly::Utils::get_check_record($sample_rec, ["sample_dir"]);
    $outrec->{read_data}->{R1} = {};
    $outrec->{read_data}->{R1}->{raw} = {};
    $outrec->{read_data}->{R1}->{raw}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R1", "rawdata"]);
    $outrec->{read_data}->{R1}->{raw}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "rawdata", "num_reads"]);
    $outrec->{read_data}->{R1}->{raw}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "rawdata", "read_length"]);
    $outrec->{read_data}->{R2} = {};
    $outrec->{read_data}->{R2}->{raw} = {};
    $outrec->{read_data}->{R2}->{raw}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R2", "rawdata"]);
    $outrec->{read_data}->{R2}->{raw}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "rawdata", "num_reads"]);
    $outrec->{read_data}->{R2}->{raw}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "rawdata", "read_length"]);

    # get info for output file path for release raw reads
    my $bio_type = Assembly::Utils::get_check_record($sample_rec, ["bio_type"]);
    my $sample_id = Assembly::Utils::get_check_record($sample_rec, ["sample"]);
    #my $release_prefix = Assembly::Utils::get_check_record($strain_rec, ["release", "prefix"]);
    #my $release_dir = Assembly::Utils::get_check_record($strain_rec, ["release", "release_dir"]);
    my $fpath_base = $release_dir . "/" . $release_prefix . "_" . $bio_type . 
            "_" . $sample_id . "_";

    my $r1_in = Assembly::Utils::get_check_record($sample_rec, ["R1", "rawdata"]);
    my $r1_out = $fpath_base . "R1.fq";
    my $r2_in = Assembly::Utils::get_check_record($sample_rec, ["R2", "rawdata"]);
    my $r2_out = $fpath_base . "R2.fq";

    $outrec->{release} = [];
    $outrec->{release}->[0] = {};
    $outrec->{release}->[0]->{input_file} = $r1_in;
    $outrec->{release}->[0]->{output_file} = $r1_out;
    $outrec->{release}->[1] = {};
    $outrec->{release}->[1]->{input_file} = $r2_in;
    $outrec->{release}->[1]->{output_file} = $r2_out;

    # copy the files here??
    unless (-e $r1_out and (-s $r1_in) == (-s $r1_out)) {
        print_verbose "cp $r1_in $r1_out\n";
        if ($options->{run}) {
            system("cp $r1_in $r1_out");
        }
    }
    unless (-e $r2_out and (-s $r2_in) == (-s $r2_out)) {
        print_verbose "cp $r2_in $r2_out\n";
        if ($options->{run}) {
            system("cp $r2_in $r2_out");
        }
    }

    if ($trimraw =~ /trim/) {
        $outrec->{read_data}->{R1}->{trim} = {};
        $outrec->{read_data}->{R1}->{trim}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R1", "trimdata"]);
        $outrec->{read_data}->{R1}->{trim}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "trimdata", "num_reads"]);
        $outrec->{read_data}->{R1}->{trim}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "trimdata", "read_length"]);
        $outrec->{read_data}->{R2}->{trim} = {};
        $outrec->{read_data}->{R2}->{trim}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R2", "trimdata"]);
        $outrec->{read_data}->{R2}->{trim}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "trimdata", "num_reads"]);
        $outrec->{read_data}->{R2}->{trim}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "trimdata", "read_length"]);

        my $tr1_in = Assembly::Utils::get_check_record($sample_rec, ["R1", "trimdata"]);
        my $tr1_out = $fpath_base . "trim_R1.fq";
        my $tr2_in = Assembly::Utils::get_check_record($sample_rec, ["R2", "trimdata"]);
        my $tr2_out = $fpath_base . "trim_R2.fq";
        
        $outrec->{release}->[2] = {};
        $outrec->{release}->[2]->{input_file} = $tr1_in;
        $outrec->{release}->[2]->{output_file} = $tr1_out;
        $outrec->{release}->[3] = {};
        $outrec->{release}->[3]->{input_file} = $tr2_in;
        $outrec->{release}->[3]->{output_file} = $tr2_out;
        
        unless (-e $tr1_out and (-s $tr1_in) == (-s $tr1_out)) {
            print_verbose "cp $tr1_in $tr1_out\n";
            if ($options->{run}) {
                system("cp $tr1_in $tr1_out");
            }
        }
        unless (-e $tr2_out and (-s $tr2_in) == (-s $tr2_out)) {
            print_verbose "cp $tr2_in $tr2_out\n";
            if ($options->{run}) {
                system("cp $tr2_in $tr2_out");
            }
        }
    }

    # my $fastqc_version = `fastqc --version`; 
    my $fastqc_version = "FastQC v0.10.1";
    # my $fastx_version = `fastx_trimmer -h | grep FASTX`;
    my $fastx_version = "Part of FASTX Toolkit 0.0.13.2 by A. Gordon (gordon\@cshl.edu)";

    # Now add the fastqc info to the pipeline
    my $qcr1 = {};
    $qcr1->{description} = "FastQC";
    $qcr1->{version} = $fastqc_version;
    chomp $qcr1->{version};
    $qcr1->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "raw_cmd"]);
    $qcr1->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "raw_qsub"]);
    $qcr1->{run_dir} = getcwd;
    
    my $qcr2 = {};
    $qcr2->{description} = "FastQC";
    $qcr2->{version} = $fastqc_version;
    chomp $qcr2->{version};
    $qcr2->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "raw_cmd"]);
    $qcr2->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "raw_qsub"]);
    $qcr2->{run_dir} = getcwd;
    
    my $qtr1 = {};
    $qtr1->{description} = "Fastx_trimmer";
    $qtr1->{version} = $fastx_version;
    chomp $qtr1->{version};
    $qtr1->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R1", "trim_cmd"]);
    $qtr1->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R1", "qsub_cmd"]);
    $qtr1->{run_dir} = getcwd;
    
    my $qtr2 = {};
    $qtr2->{description} = "Fastx_trimmer";
    $qtr2->{version} = $fastx_version;
    chomp $qtr2->{version};
    $qtr2->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R2", "trim_cmd"]);
    $qtr2->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R2", "qsub_cmd"]);
    $qtr2->{run_dir} = getcwd;    
    
    # Now add the fastqc info to the pipeline
    my $qct1 = {};
    $qct1->{description} = "FastQC";
    $qct1->{version} = $fastqc_version;
    chomp $qct1->{version};
    $qct1->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "trim_cmd"]);
    $qct1->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "trim_qsub"]);
    $qct1->{run_dir} = getcwd;
    #$qct1->{release} = [];
    #$qct1->{release}->[0] = {};
    #$qct1->{release}->[0]->{input_file} = Assembly::Utils::get_check_record($sample_rec, ["R1", "trimdata"]);
    #$qct1->{release}->[0]->{output_file} = 
    # this is handled in the location above... should it be??
    
    my $qct2 = {};
    $qct2->{description} = "FastQC";
    $qct2->{version} = $fastqc_version;
    chomp $qct2->{version};
    $qct2->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "trim_cmd"]);
    $qct2->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "trim_qsub"]);
    $qct2->{run_dir} = getcwd;
    
    my $qc =[];
    push (@$qc, $qcr1);
    push (@$qc, $qcr2);
    if ($trimraw =~ /trim/) {
        push (@$qc, $qtr1);
        push (@$qc, $qtr2);
        push (@$qc, $qct1);
        push (@$qc, $qct2);
    }
    
    return ($outrec, $qc);
}



sub create_pipeline_stanza
{
    my $asm_rec = shift;
    my $name_prefix = shift;
    my $description = shift;
    my $release_stanza = {};
    my %lookup = ();
    for my $key (qw(cmd qsub_cmd version)) {
        $lookup{$key} = $name_prefix . "_" . $key;
    }
    print Assembly::Utils::get_check_record($asm_rec, ["bowtie_cmd"]) . "\n";
    $release_stanza->{"command"} = Assembly::Utils::get_check_record($asm_rec, [$lookup{"cmd"}]);
    $release_stanza->{"qsub_cmd"} = Assembly::Utils::get_check_record($asm_rec, [$lookup{"qsub_cmd"}]);
    $release_stanza->{"version"} = Assembly::Utils::get_check_record($asm_rec, [$lookup{"version"}]);
    $release_stanza->{"run_dir"} = getcwd;
    $release_stanza->{"description"} = $description;
    return $release_stanza;
}

sub update_release
{
    my ($release_rec, $yaml_records, $table_records, $sample, $strain, $trimraw) = @_;
        
    my $species = $table_records->{$sample}->{"Species"};
    my $strain_rec = Assembly::Utils::get_check_record ($yaml_records, [$species, "RNA", $strain]); 
    my $output_release_prefix = $table_records->{$sample}->{"Output_Release_Prefix"};
    my $output_release_dir = $table_records->{$sample}->{"Output_Release_Dir"};
    my ($sample_stanza, $qc_cmds) = get_sample_stanza ($strain_rec, $sample, $trimraw, $output_release_prefix, $output_release_dir);   
    push (@{$release_rec->[2]->{samples}}, $sample_stanza);
    my $pipeline_ref = $release_rec->[3]->{pipeline};
        for my $rec (@$qc_cmds) {
        push (@$pipeline_ref, $rec);
    }
    
    my $yrec = Assembly::Utils::get_check_record ($yaml_records, [$species, "RNA", $strain, $sample, "RNA_assembly", $trimraw]);
    my $bowtie_stanza = create_pipeline_stanza ($yrec, "bowtie", "bowtie2-build indexer");
    my $tophat_stanza = create_pipeline_stanza ($yrec, "tophat", "tophat alignment");
    my $cufflinks_stanza = create_pipeline_stanza ($yrec, "cufflinks", "cufflinks assembly");
    # Add the input/output file info for the released transcripts.gtf file
    # Copy the transcripts file to the release dir.
    my $cufflinks_dir = Assembly::Utils::get_check_record ($yrec, ["cufflinks_dir"]);
    my $transcripts_in = $cufflinks_dir . "/transcripts.gtf";
    my $transcripts_out = $output_release_dir . "/" . $output_release_prefix . "_" . $sample . "_transcripts.gtf";
    
    push (@$pipeline_ref, $bowtie_stanza);
    push (@$pipeline_ref, $tophat_stanza);
    push (@$pipeline_ref, $cufflinks_stanza);
}
    
sub samples_by_strain
{
    my $table_records = shift;
    my $table_strains = {};
    for my $sample (keys %$table_records) {
        my $strain = $table_records->{$sample}->{"Strain"};
        push (@{$table_strains->{$strain}}, $sample);
    }
    return $table_strains;
}

sub get_output_filename
{
    my $table_records = shift;
    my $sample = shift;
    my $output_release_prefix = $table_records->{$sample}->{"Output_Release_Prefix"};
    my $output_release_dir = $table_records->{$sample}->{"Output_Release_Dir"};
    create_dir ($output_release_dir);
    my $output_release_filename = $output_release_dir . "/" .
        $output_release_prefix . "_metadata.yml";
    return $output_release_filename;
}

sub create_release
{
    my $yaml_records = shift;
    my $table_records = shift;
    my $table_strains = samples_by_strain ($table_records);
    
    for my $strain (keys %$table_strains) {
        my @sample_list = @{$table_strains->{$strain}};
        my $release_rec = '';
        my $last_sample = '';
        for my $sample (@sample_list) {
            print "Working on sample $sample\n";
            my $input_release_filename = $table_records->{$sample}->{"Reference_Metafile"};
            unless ($release_rec) {
                $release_rec = LoadFile ($input_release_filename);
            }
            my $trimraw = $table_records->{$sample}->{"Trim_Raw"};
            update_release ($release_rec, $yaml_records, $table_records, $sample, $strain, $trimraw);
            $last_sample = $sample;
        }
        print "sample $last_sample\n";
        my $output_release_filename = get_output_filename ($table_records, $last_sample);
        print "Writing file $output_release_filename\n";
        DumpFile ($output_release_filename, $release_rec);
    }
}

gather_options;
my $yaml_records = LoadFile ($options->{yaml_in});
my $table_records = parse_assembly_table;
create_release ($yaml_records, $table_records);
#DumpFile ($options->{yaml_out}, $yaml_records);
























