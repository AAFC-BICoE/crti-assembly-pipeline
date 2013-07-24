#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
# output files created by velvetg
my $vg_outfiles = [qw(Graph2 LastGraph PreGraph stats.txt)];
my @col_headers = ("Sample", "Trim/raw", "Total num reads (M)", "Avg. read len (bp)",
        "Min kmer", "Max kmer", "Best kmer", "Best N50", "VelvetK kmer", "VelvetK N50",
        "Advisor kmer", "Advisor N50"); 

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/10_velvetg_qsub.yml
        yaml_out yaml_files/11_velvet_stats.yml
        stats_outfile output_files/VelvetStats.tab
        verbose 0
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    $options->{qsub_opts} = $options->{qsub_opts} . ""; # . -N velvet_hg ";
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --stats_outfile <filename>
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
        'stats_outfile=s',
        );
    set_default_opts;
    check_opts;
}

sub get_check_record
{
    my $ref = shift;
    my $kref = shift;
    for my $key (@$kref) {
        if (ref($ref) eq "HASH" and defined ($ref->{$key})) {
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

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $tr = shift;
    my $kmin = get_check_record($rec, ["velvet", $tr, "min_kmer"]);
    my $kmax = get_check_record($rec, ["velvet", $tr, "max_kmer"]);
    my $krange = [];
    if ($kmin =~ /^\d+$/ and $kmax =~ /^\d+$/) {
        for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
            push (@$krange, $kmer);
        }
    }
    return $krange;
}

sub parse_log_n50
{
    my $rec = shift;
    my $tr = shift;
    my $kmer = shift;
    my $kmer_dir = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "kmer_dir"]);
    my $log_file = $kmer_dir . "/Log";
    my $contig_file = $kmer_dir . "/contigs.fa";
    my $n50 = '';
    if (-e $log_file and -e $contig_file and -s $contig_file > 0) {
        open (FIN, '<', $log_file) or die "Error: couldn't open file ${log_file}\n";
        while (<FIN>) { 
            if (/n50\s+of\s+(\d+)/i) {
                $n50 = $1;
                print "kmer $kmer - found n50 " . $n50 . "\n";
            }
        }
        close (FIN);
    }
    return $n50;
}

sub get_max_kmer
{
    my $rec = shift;
    my $tr = shift;
    my $kmer_range = get_kmer_range($rec, $tr);
    if (scalar @$kmer_range > 0) {
        my %max;
        $max{kmer} = $kmer_range->[0];
        $max{n50} = 0;
        for my $kmer (@$kmer_range) {
            my $n50 = parse_log_n50($rec, $tr, $kmer);
            if ($n50) {
                set_check_record($rec, ["velvet", $tr, "kmer", $kmer], "N50", $n50);
                if ($n50 > $max{n50}) {
                    $max{n50} = $n50;
                    $max{kmer} = $kmer;
                }
            }
        }
        if ($max{n50} > 0) {
            my $sample = $rec->{sample};
            print "For sample $sample $tr got max n50 kmer " . $max{kmer} . " value " . $max{n50} . "\n";
            set_check_record($rec, ["velvet", $tr], "max_n50_kmer", $max{kmer});
            set_check_record($rec, ["velvet", $tr], "max_n50_value", $max{n50});
        }
    }
}

sub get_all_max
{
    my $records = shift;
    for my $sample (keys %$records) {
        print $sample . "\n";
        my $rec = $records->{$sample};
        for my $tr (qw(trim raw)) {
            get_max_kmer($rec, $tr);
        }
    }
}   

sub get_stats_line
{
    my $rec = shift;
    my $trimraw = shift;
    my $trimdata = $trimraw . "data";
    my $sample = ''; #$rec->{sample};
    
    my $r1_numreads = get_check_record($rec, ["data_stats", "R1", $trimdata, "num_reads"]);
    my $r2_numreads = get_check_record($rec, ["data_stats", "R2", $trimdata, "num_reads"]);
    
    my $total_numreads = '';
    if ($r1_numreads and $r2_numreads) {
        $total_numreads = ($r1_numreads + $r2_numreads) / 1000000;
    }
    my $avg_readlen = get_check_record($rec, ["velvet", $trimraw, "average_read_length"]);
    my $genome_len = get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
    if ($genome_len) {
        $genome_len = $genome_len / 1000000;
    }
    
    my $advisor_kmer = get_check_record($rec, ["velvet", $trimraw, "velvet_advisor_best_kmer"]);
    my $advisor_n50 ='';
    if ($advisor_kmer) {
        $advisor_n50 = get_check_record($rec, ["velvet", $trimraw, "kmer", $advisor_kmer, "n50"]);
    }

    my $vk_kmer = get_check_record($rec, ["velvet", $trimraw, "velvetk_best_kmer"]);
    my $vk_n50 = '';
    if ($vk_kmer) {
        $vk_n50 = get_check_record($rec, ["velvet", $trimraw, "kmer", $vk_kmer, "n50"]);
    }
    
    my $best_kmer = get_check_record($rec, ["velvet", $trimraw, "max_n50_kmer"]);
    my $best_n50 = get_check_record($rec, ["velvet", $trimraw, "max_n50_value"]);

    my $kmin = get_check_record($rec, ["velvet", $trimraw, "min_kmer"]);
    my $kmax = get_check_record($rec, ["velvet", $trimraw, "max_kmer"]);
    #my $missing_kmers = get_check_record($rec, [ ]);
    
    my $out_str = join ("\t", ($sample, $trimraw, $total_numreads, $avg_readlen, $kmin, $kmax,
            $best_kmer, $best_n50, $vk_kmer, $vk_n50, $advisor_kmer, $advisor_n50));
    
    return $out_str;
}

sub write_stats_file
{
    my $records = shift;
    my $fname = ($options->{stats_outfile} ? $options->{stats_outfile} : '');
    if ($fname) {
        open (FSTATS, '>', $fname) or die "Error: couldn't open output stats file $fname.\n";
        print FSTATS join ("\t", @col_headers) . "\n";
        for my $sample (keys %$records) {
            my $rec = $records->{sample};
            for my $trimraw (qw(trim raw)) {
                my $stats_line = get_stats_line ($rec, $trimraw);
                print FSTATS $stats_line . "\n";
            }
        }
        close (FSTATS);
    }
}
    
sub run_all
{          
    gather_opts;
    my $records = LoadFile($options->{yaml_in});
    get_all_max($records);
    write_stats_file($records);
    DumpFile($options->{yaml_out}, $records);
}

run_all;







