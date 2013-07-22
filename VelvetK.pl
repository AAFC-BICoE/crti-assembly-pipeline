#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my $velvetk_bin = "./velvetk.pl";
my $tr;
my $trdata;

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/07_velveth_qsub.yml
        yaml_out yaml_files/08_velvetg_qsub.yml
        trim 1
        verbose 0
        velvetk_infile input_data/VelvetKBest.tab
        velvetk_outfile input_data/VelvetKBestOut.tab
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    if ($options->{raw}) {
        $tr = "raw";
        $trdata = "rawdata";
    } else {
        $tr = "trim";
        $trdata = "trimdata";
    }
    $options->{qsub_opts} = $options->{qsub_opts} . ""; # . -N velvet_hg ";
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
                --velvetk_infile <filename>
                --velvetk_outfile <filename>
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
        'sample_list=s',
        'trim',
        'raw',
        'velvetk_infile=s',
        'velvetk_outfile=s',
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

sub set_check_record
{
    my $ref = shift;
    my $kref = shift;
    my $last_key = shift;
    my $value = shift;
    for my $key (@$kref) {
        if (defined ($ref->{$key})) {
            $ref = $ref->{$key};
        } else {
            $ref->{$key} = {};
            $ref = $ref->{$key};
        }
    }
    $ref->{$last_key} = $value;
}

sub get_sample_list
{
    my $records = shift;
    my $sample_list = [];
    if ($options->{sample_list}) {
        $sample_list = [split(/,/, $options->{sample_list})];
    } else {
        $sample_list = [keys %$records];
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
            if (scalar @fields == 4) {
                my ($sample, $tr, $velvetk_cmd, $best_kmer) = @fields;
                my $rec = $records{$sample};
                set_check_record($rec, ["velvet", $tr], "velvetk_cmd", $velvetk_cmd);
                set_check_record($rec, ["velvet", $tr], "velvetk_best_kmer", $best_kmer);
            }
        }
    }
}

sub get_velvetk_cmd
{
    my $rec = shift;
    my $genome_len = get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
    my $r1data = get_check_record($rec, ["R1", $trdata]);
    my $r2data = get_check_record($rec, ["R2", $trdata]);
    my $velvetk_cmd = $velvetk_bin . " --size " . $genome_len . " --best $r1data $r2data";
    set_check_record($rec, ["velvet", $tr], "velvetk_cmd", $velvetk_cmd);
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
my $sample_list = get_sample_list($records);
for my $sample (@$sample_list) {
    my $rec = $records->{$sample};
    my $have_best = get_check_record($rec, ["velvet", $tr, "velvetk_best_kmer"]);
    unless ($have_best) {
        my $vk_cmd = get_velvetk_cmd($rec);
        print_verbose "Running command:\n" . $vk_cmd . "\n";
        #my $best = 0;
        my $best = `$vk_cmd`;
        chomp $best;
        print_verbose "velvetk.pl found best kmer: " . $best . "\n";
        set_check_record($rec, ["velvet", $tr], "velvetk_best_kmer", $best);
    } else {
        print_verbose "Already found best kmer for sample $sample trim/raw $tr. Best is: " . $have_best . "\n";
    }   
}
DumpFile($options->{yaml_out}, $records);


