#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my @col_headers = ("Sample", "Trim/raw", "Num reads (M)", "Avg read length (bp)", 
        "Est. Genome Size (Mbp)", "Advisor k-mer");

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/06_genome_lengths.yml
        yaml_out yaml_files/07_velvet_advisor.yml
        trim 1
        verbose 0
        advisor_infile input_data/VelvetAdvisorBest.tab
        advisor_outfile output_files/VelvetAdvisorBestOut.tab
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
                --advisor_infile <filename>
                --advisor_outfile <filename>
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
    my $fname = ($options->{advisor_infile} ? $options->{advisor_infile} : '');
    if ($fname) {
        open (FTAB, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FTAB>; #skip first line.
        while (my $line = <FTAB>) {
            chomp $line;
            print "Line $.\n";
            my @fields = split(/\t/, $line);
            if (scalar @fields == 6) {
                print "Parsing!\n";
                my ($sample, $trimraw, $num_reads, $avg_readlen, $est_gen_len, $advisor_best_kmer) = @fields;
                my $rec = $records->{$sample};
                set_check_record($rec, ["velvet", $trimraw], "velvet_advisor_best_kmer", $advisor_best_kmer);
            } else {
                print_verbose "Couldn't parse input line $.:\n$line\n";
            }
        }
    }
}

sub pull_genome_data
{
    my $rec = shift;
    my $sample = shift;
    my $tr = shift;
    my $trdata = $tr . "data";
    my $advisor_kmer = get_check_record($rec, ["velvet", $tr, "velvet_advisor_best_kmer"]);
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
    my $out_str = join ("\t", ($sample, $tr, $total_numreads, $avg_readlen, $genome_len, $advisor_kmer));
    return $out_str;
}

sub write_genome_data
{
    my $records = shift;
    my $sample_list = shift;
    my $fname = ($options->{advisor_outfile} ? $options->{advisor_outfile} : '');
    if ($fname) {
        open (FOUT, '>', $fname) or die "Error: couldn't open output file $fname.\n";
        print FOUT join ("\t", @col_headers) . "\n";
        for my $sample (@$sample_list) {
            my $rec = $records->{$sample};
            for my $tr (qw(trim raw)) {
                my $out_line = pull_genome_data($rec, $sample, $tr);
                print FOUT $out_line . "\n";
            }
        }
    }
    close (FOUT);
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
my $sample_list = get_sample_list($records);
parse_input_table($records);
write_genome_data;
DumpFile($options->{yaml_out}, $records);

