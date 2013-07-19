#!/usr/bin/env perl
use strict;
use warnings;
use YAML::XS qw(LoadFile DumpFile);
use Getopt::Long;
use Date::Parse;
use Date::Calc;

my $yaml_in = "yaml_files/09_velvet_n50.yml";
my $tr = "trim";
my @completed_jobs = ();

GetOptions('yaml_in|i=s' => \$yaml_in);

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

my $i=0;
my $records = {};
$records = LoadFile($yaml_in);
for my $sample (keys %$records) {
    my $rec = $records->{$sample};
    my $krange = get_kmer_range($rec);
    for my $kmer (@$krange) {
        my $n50 = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "N50"]);
        my $qsub_jobid = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "velvetg_qsub_jobid"]);
        if ($n50 and $qsub_jobid) {
            $i++;
            print "sample\t$sample\ntrimraw\t$tr\nkmer\t$kmer\nn50\t$n50\n";
            my $glen = get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
            my $numreads = get_check_record($rec, ["data_stats", "R1", $tr . "data", "num_reads"]);
            my $readlen = get_check_record($rec, ["data_stats", "R1", $tr . "data", "read_length"]);
            print "genome length\t$glen\nnum reads (R1)\t$numreads\nread length (R1)\t$readlen\n";
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
            if ($start_time and $end_time) {
                my $diff = ($end_time - $start_time);
                my $hours = int($diff/3600);
                my $minutes = int( ($diff - $hours * 3600)/60 );
                my $seconds = $diff - $hours * 3600 - $minutes * 60;
                print "Time difference: " . $hours . " hr " . $minutes . " min " . $seconds . " sec \n";
            }
            print "\n";     
        }
    }
}
print "record count: $i\n";
#my $comp_rec = [];
#            push (@$comp_rec, [$sample, $kmer]);
#            push (@completed_jobs, $comp_rec);
#for my $arr (@completed_jobs) {
#    my ($sample, $kmer) = 