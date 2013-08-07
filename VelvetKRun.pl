#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);
use Assembly::Utils;

my $options = {};
my $velvetk_bin = "./velvetk.pl";

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/07_velvet_advisor.yml
        yaml_out yaml_files/08_velvetk.yml
        trim 1
        raw 0
        verbose 0
        run 0
        velvetk_infile input_data/VelvetKBest.tab
        velvetk_outfile output_files/VelvetKBestOut.tab
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --sample_list <ID1,ID2,...,IDX (no spaces)>
                --trim
                --raw
                --run
                --velvetk_infile <filename>
                --velvetk_outfile <filename>
                ";
    }
}

sub gather_opts
{
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'verbose',
        'sample_list=s',
        'trim',
        'raw',
        'velvetk_infile=s',
        'velvetk_outfile=s',
        'run',
        );
    set_default_opts;
    check_opts;
}

sub get_sample_list
{
    my $records = shift;
    my $sample_list = [];
    if ($options->{sample_list}) {
        $sample_list = [split(/,/, $options->{sample_list})];
    } else {
        for my $sample (keys %$records) {
            if ($records->{$sample}->{bio_type} =~ /DNA/) {
                push (@$sample_list, $sample);
            }
        }
    }
    return $sample_list;
}

sub print_verbose
{
    if ($options->{verbose}) {
        print (@_);
    }
}

sub parse_input_table
{
    my $records = shift;
    my $fname = ($options->{velvetk_infile} ? $options->{velvetk_infile} : '');
    if ($fname) {
        open (FTAB, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FTAB>; #skip first line.
        while (my $line = <FTAB>) {
            chomp $line;
            my @fields = split(/\t/, $line);
            if (scalar @fields == 7) {
                my ($species, $strain, $bio_type, $sample, $trimraw, $best_kmer, $velvetk_cmd) = @fields;
                $species = Assembly::Utils::format_species_key($species);
                $strain = Assembly::Utils::format_strain_key($strain);
                print "species $species strain $strain sample $sample\n";
                my $sample_list = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples"]);
                if ($sample_list and scalar (keys (%$sample_list)) == 1) {
                    print "setting\n";
                    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_cmd", $velvetk_cmd);
                    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer", $best_kmer);
                }
            }
        }
    }
}

sub write_output_table
{
    my $records = shift;
    my $sample = '';
    my $fname = ($options->{velvetk_outfile} ? $options->{velvetk_outfile} : '');
    if ($fname) {
        open (FTAB, '>', $fname) or die "Error: couldn't open file $fname\n";
        print FTAB join("\t", qw(Species Strain Trim/Raw VK_Best_Kmer VK_Command)) . "\n";
        for my $species (keys %$records) {
            my $ss = $records->{$species}->{DNA};
            for my $strain (keys %$ss) {
                for my $trimraw (qw(trim raw)) {
                    my $vk_cmd = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_cmd"]);
                    my $vk_best = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_best_kmer"]);
                    print FTAB join("\t", ($species, $strain,$trimraw, $vk_best, $vk_cmd)) . "\n";
                }
            }
        }
    }
}
            

sub get_velvetk_cmd
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $velvetk_cmd = '';
    my $trimdata = $trimraw . "data";
    my $genome_len = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "related_genome_length", "RG_Est_Genome_Length"]);
    my $sample_hash = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples"]);
    if ($sample_hash and scalar (keys (%$sample_hash)) == 1) {
        my @sample_list = keys %$sample_hash;
        my $sample = $sample_list[0];
        my $r1data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "R1", $trimdata]);
        my $r2data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "R2", $trimdata]);
        $velvetk_cmd = $velvetk_bin . " --size " . $genome_len . " --best $r1data $r2data";
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_cmd", $velvetk_cmd);
    }
    return $velvetk_cmd;
}

sub get_velvetk_sample
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    
    my $have_best = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_best_kmer"]);
    unless ($have_best) {
        my $vk_cmd = get_velvetk_cmd($records, $species, $strain, $trimraw);
        print_verbose "Running command:\n" . $vk_cmd . "\n";
        #my $best = 0;
        if ($options->{run} and $vk_cmd) {
            my $best = `$vk_cmd`;
            chomp $best;
            print_verbose "velvetk.pl found best kmer: " . $best . "\n";
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer", $best);
        }
    } else {
        print_verbose "Already found best kmer for species $species strain $strain trim/raw $trimraw. Best is: " . $have_best . "\n";
    }
} 

sub run_velvetk
{
    my $records = shift;
    #my $sample_list = shift;
    for my $species (keys %$records) {
        my $sr = $records->{$species}->{DNA};
        for my $strain (keys %$sr) {
            my $sl = $sr->{$strain};
            for my $trimraw (qw(trim raw)) {
                if ($options->{$trimraw}) {
                    get_velvetk_sample($records, $species, $strain, $trimraw);
                }
            }   
        }
    }
}  

sub run_all
{
    gather_opts;
    my $records = LoadFile($options->{yaml_in});
    parse_input_table($records);
    #my $sample_list = get_sample_list($records);
    run_velvetk($records);
    write_output_table($records);
    DumpFile($options->{yaml_out}, $records);
}

run_all;


