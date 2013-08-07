#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);
use Assembly::Utils;

my $options = {};
my @col_headers = ("Species", "Strain", "Bio_type", "Sample", "Trim/raw", "Num reads (M)", "Avg read length (bp)", 
        "Est. Genome Size (Mbp)", "Advisor k-mer");

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/07_genome_lengths.yml
        yaml_out yaml_files/07_velvet_advisor.yml
        trim 1
        verbose 0
        advisor_infile input_data/VelvetAdvisorKmers.tab
        advisor_outfile output_files/VelvetAdvisorKmerOut.tab
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
            my @fields = split(/\t/, $line);
            if (scalar @fields == 9) {
                my ($species, $strain, $bio_type, $sample, $trimraw, $num_reads, $avg_readlen, $est_gen_len, $advisor_best_kmer) = @fields;
                if ($species) {
                    $species = Assembly::Utils::format_species_key($species);
                }
                $strain = Assembly::Utils::format_strain_key($strain);
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvet_advisor_best_kmer", $advisor_best_kmer);
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
    my $advisor_kmer = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "velvet_advisor_best_kmer"]);
    my $r1_numreads = Assembly::Utils::get_check_record($rec, ["data_stats", "R1", $trdata, "num_reads"]);
    my $r2_numreads = Assembly::Utils::get_check_record($rec, ["data_stats", "R2", $trdata, "num_reads"]);
    my $total_numreads = '';
    if ($r1_numreads and $r2_numreads) {
        $total_numreads = ($r1_numreads + $r2_numreads) / 1000000;
    }
    my $avg_readlen = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "average_read_length"]);
    my $genome_len = Assembly::Utils::get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
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
        close (FOUT);
    }
}

sub run_all
{
    gather_opts;
    my $records = LoadFile($options->{yaml_in});
    #my $sample_list = get_sample_list($records);
    parse_input_table($records);
    #write_genome_data($records, $sample_list);
    DumpFile($options->{yaml_out}, $records);
}

run_all;
