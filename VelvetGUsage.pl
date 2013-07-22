#!/usr/bin/env perl
use strict;
use warnings;
use YAML::XS qw(LoadFile DumpFile);
use Getopt::Long;
use Date::Parse;
use Date::Calc;

my $yaml_in = "yaml_files/09_velvet_n50.yml";
my $stats_file = "output_files/velvetg_run_stats.txt";
my $tr = "trim";
my @completed_jobs = ();

GetOptions('yaml_in|i=s' => \$yaml_in,
        'stats_file|s=s' => \$stats_file);

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

sub get_kmer_range
{
    my $rec = shift;
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

sub print_stats
{
    my $rec = shift;
    my $kmer = shift;
    my $n50 = shift;
    my $sample = $rec->{sample};
    print "sample:\t$sample\ntrim/raw:\t$tr\nkmer:\t$kmer\nn50:\t$n50\n";
    my $glen = get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
    my $numreadsR1 = get_check_record($rec, ["data_stats", "R1", $tr . "data", "num_reads"]);
    my $numreadsR2 = get_check_record($rec, ["data_stats", "R2", $tr . "data", "num_reads"]);
    my $readlenR1 = get_check_record($rec, ["data_stats", "R1", $tr . "data", "read_length"]);
    my $readlenR2 = get_check_record($rec, ["data_stats", "R2", $tr . "data", "read_length"]);
    print "genome length:\t$glen:\nR1 numreads:\t$numreadsR1\nR1 read length:\t$readlenR1\n";
    print "R2 numreads:\t$numreadsR2\nR2 read length:\t$readlenR2\n";
}

sub print_time_difference
{
    my $start_time = shift;
    my $end_time = shift;
    if ($start_time and $end_time) {
        my $diff = ($end_time - $start_time);
        my $hours = int($diff/3600);
        my $minutes = int( ($diff - $hours * 3600)/60 );
        my $seconds = $diff - $hours * 3600 - $minutes * 60;
        print "Time difference: " . $hours . " hr " . $minutes . " min " . $seconds . " sec \n";
    }
}

sub print_qacct_stats
{
    my $qsub_jobid = shift;
    my $qacct_str = `qacct -j ${qsub_jobid}`;
    my @qacct_lines = split(/\n/, $qacct_str);
    my ($start_time, $end_time) = ('','');
    for my $line (@qacct_lines) {
        chomp $line;
        if ($line =~ /^(hostname|start_time|end_time|mem|maxvmem)\s+(.*)$/) {
            print $line . "\n";
            my $key = $1;
            my $val = $2;
            if ($key =~ /start_time/) {
                $start_time = str2time($val);
            }
            if ($key =~ /end_time/) {
                $end_time = str2time($val);
            }
        }
    }
    print_time_difference($start_time, $end_time);
    print "\n"; 
}     

sub print_sample_stats
{
    my $rec = shift;
    my $krange = get_kmer_range($rec);
    my $num_printed = 0;
    for my $kmer (@$krange) {
        my $n50 = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "N50"]);
        my $qsub_jobid = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "velvetg_qsub_jobid"]);
        if ($n50 and $qsub_jobid) {
            print_stats($rec, $kmer, $n50);
            print_qacct_stats($qsub_jobid);
            $num_printed++;
        }
    }
    return $num_printed;
}

my $records = {};
my $i = 0;
$records = LoadFile($yaml_in);
for my $sample (keys %$records) { 
    my $num_printed = print_sample_stats($records->{$sample});
    $i += $num_printed;
}
print "Number of completed records found: " . $i . "\n";
