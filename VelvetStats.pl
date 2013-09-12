#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);
use Assembly::Utils;

my $options = {};
# output files created by velvetg
my $vg_outfiles = [qw(Graph2 LastGraph PreGraph stats.txt)];
my @col_headers = ("Species", "Strain", "Trim/raw", "VelvetK", "Best Kmer", "Best N50", "Max Contig", "Total Length", "Reads Used", "Missing Kmers");

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/11_velvet_cmds.yml
        yaml_out yaml_files/12_velvet_stats.yml
        stats_outfile output_files/VelvetStats.tab
        verbose 0
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
                --stats_outfile <filename>
                ";
    }
}

sub gather_opts
{
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'verbose',
        'stats_outfile=s',
        );
    set_default_opts;
    check_opts;
}

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $tr = shift;
    my $kmin = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "min_kmer"]);
    my $kmax = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "max_kmer"]);
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
    my $kmer_dir = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "kmer_dir"]);
    my $log_file = $kmer_dir . "/Log";
    my $contig_file = $kmer_dir . "/contigs.fa";
    my $n50 = '';
    my $max_len = '';
    my $total_len = '';
    my $reads_used = '';
    if (-e $log_file and -e $contig_file and -s $contig_file > 0) {
        open (FIN, '<', $log_file) or die "Error: couldn't open file ${log_file}\n";
        while (<FIN>) { 
            if (/n50\s+of\s+(\d+)/i) {
                $n50 = $1;
                #print "kmer $kmer - found n50 " . $n50 . "\n";
            }
            if (/n50 of (\d+), max (\d+), total (\d+), using ([0-9\/]+) reads/) {
                ($n50, $max_len, $total_len, $reads_used) = ($1, $2, $3, $4);
            }
        }
        close (FIN);
    }
    return ($n50, $max_len, $total_len, $reads_used);
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
            my ($n50, $max_contig_len, $total_len, $reads_used) = parse_log_n50($rec, $tr, $kmer);
            if ($n50) {
                Assembly::Utils::set_check_record($rec, ["velvet", $tr, "kmer", $kmer], "N50", $n50);
                if ($n50 > $max{n50}) {
                    $max{n50} = $n50;
                    $max{kmer} = $kmer;
                    $max{max_contig_len} = $max_contig_len;
                    $max{total_len} = $total_len;
                    $max{reads_used} = $reads_used;
                }
                #if ($n50 < 2000) {
                 #   print $rec->{species} . "\t" . $rec->{sample} . "\t$tr\t" . $kmer . "\t" . $n50 . "\t";
                 #   my $kdir = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "kmer_dir"]);
                    #print $kdir . "\t";
                    #my $du = `du -h $kdir`;
                    #print "Disk usage: $du\n";
                #}
            }
        }
        if ($max{n50} > 0) {
            Assembly::Utils::set_check_record($rec, ["velvet", $tr], "max_n50_kmer", $max{kmer});
            Assembly::Utils::set_check_record($rec, ["velvet", $tr], "max_n50_value", $max{n50});
            Assembly::Utils::set_check_record($rec, ["velvet", $tr], "max_n50_max_contig", $max{max_contig_len});
            Assembly::Utils::set_check_record($rec, ["velvet", $tr], "max_n50_total_len", $max{total_len});
            Assembly::Utils::set_check_record($rec, ["velvet", $tr], "max_n50_reads_used", $max{reads_used});
        }
    }
}

sub get_all_max
{
    my $records = shift;
    for my $species (keys %$records) {
        for my $strain (keys %{$records->{$species}->{DNA}}) {
            my $rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain]);
            for my $tr (qw(trim raw)) {
                get_max_kmer($rec, $tr);
            }
        }
    }
}   

sub get_missing_kmers
{
    my $rec = shift;
    my $trimraw = shift;
    my $kmin = shift;
    my $kmax = shift;
    my @missing = ();
    if ($kmin and $kmax) {
        for (my $k=$kmin; $k<=$kmax; $k+= 2) {
            my $n50 = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "kmer", $k, "N50"]);
            unless ($n50 =~ /\d/) {
                push (@missing, $k);
            }
        }
    }
    return join(",", @missing);
}

sub get_stats_line
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $trimdata = $trimraw . "data";

    my $rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain]);
    my $sample = Assembly::Utils::get_check_record($rec, ["PE", "sample"]);

    my $vk_kmer = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "velvetk_best_kmer"]);
    
    my $best_kmer = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "max_n50_kmer"]);
    my $best_n50 = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "max_n50_value"]);
    my $best_n50_max_contig_len = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "max_n50_max_contig"]);
    my $best_n50_total_len = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "max_n50_total_len"]);
    my $best_n50_reads_used = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "max_n50_reads_used"]);
    
    my $kmin = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "min_kmer"]);
    my $kmax = Assembly::Utils::get_check_record($rec, ["velvet", $trimraw, "max_kmer"]);
    my $missing_kmers = get_missing_kmers($rec, $trimraw, $kmin, $kmax);
    
    my $out_str = join ("\t", ($species, $strain, $sample, $trimraw, $vk_kmer, $best_kmer, $best_n50, $best_n50_max_contig_len, 
            $best_n50_total_len, $best_n50_reads_used, $missing_kmers));
    
    return $out_str;
}

sub write_stats_file
{
    my $records = shift;
    my $fname = ($options->{stats_outfile} ? $options->{stats_outfile} : '');
    if ($fname) {
        open (FSTATS, '>', $fname) or die "Error: couldn't open output stats file $fname.\n";
        print FSTATS join ("\t", @col_headers) . "\n";
        for my $species (keys %$records) {
            for my $strain (keys %{$records->{$species}->{DNA}}) {
                for my $trimraw (qw(trim raw)) {
                    my $stats_line = get_stats_line ($records, $species, $strain, $trimraw);
                    print FSTATS $stats_line . "\n";
                }
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







