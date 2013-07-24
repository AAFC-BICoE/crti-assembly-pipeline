#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my $velvetk_bin = "./velvetk.pl";
#my $tr;
#my $trdata;

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/07_velveth_qsub.yml
        yaml_out yaml_files/08_velvetg_qsub.yml
        trim 1
        verbose 0
        advisor_infile input_data/VelvetAdvisorBest.tab
        velvetk_outfile output_files/VelvetAdvisorBestOut.tab
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    #if ($options->{raw}) {
    #    $tr = "raw";
    #    $trdata = "rawdata";
    #} else {
    #    $tr = "trim";
    #    $trdata = "trimdata";
    #}
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
                --advisor_infile <filename>
                --advisor_outfile <filename>
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
        'advisor_infile=s',
        'advisor_outfile=s',
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
            if (scalar @fields == 3) {
                my ($sample, $trimraw, $advisor_best_kmer) = @fields;
                my $rec = $records->{$sample};
                set_check_record($rec, ["velvet", $trimraw], "velvet_advisor_best_kmer", $advisor_best_kmer);
            }
        }
    }
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
my $sample_list = get_sample_list($records);
parse_input_table($records);

print join ("\t", (qw`Sample Trim/raw Num_reads_(M) Avg_readlen_(bp) Est_genome_size_(Mbp) Advised_kmer`)) . "\n";
for my $sample (@$sample_list) {
    my $rec = $records->{$sample};
    for my $tr (qw(trim raw)) {
        my $trdata = $tr . "data";
        my $advised_kmer = get_check_record($rec, ["velvet", $tr, "velvet_advisor_best_kmer"]);
        unless ($advised_kmer) {
            my $r1_numreads = get_check_record($rec, ["data_stats", "R1", $trdata, "num_reads"]);
            my $r2_numreads = get_check_record($rec, ["data_stats", "R2", $trdata, "num_reads"]);
            my $total_numreads = '';
            if ($r1_numreads and $r2_numreads) {
                $total_numreads = ($r1_numreads + $r2_numreads) / 1000000;
            }
            my $avg_readlen = get_check_record($rec, ["velvet", $tr, "average_read_length"]);
            my $genome_len = get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
            if ($genome_len) {
                $genome_len = $genome_len / 1000000;
            }
            #print "Sample: $sample\n";
            #print "Trim/raw: $tr\n";
            #print "Number of reads (millions): " . $total_numreads . "\n";
            #print "Average read length (bp): " . $avg_readlen . "\n";
            #print "Est genome size (millions of bp): " . $genom_len . "\n";
            print join ("\t", ($sample, $tr, $total_numreads, $avg_readlen, $genome_len)) . "\n";
        }
    }
}



