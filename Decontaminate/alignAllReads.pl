#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;

# Given a set of input reads, and an input genome,
# align reads to the genome and return the set of reads
# that aligned and the set of reads that did not align.

my $alignReads_path = getcwd . "/alignReads.pl";
my $qsub_path = "/opt/gridengine/bin/linux-x64/qsub";
my $qsub_script = "qsub_script.sh";

my @alignment_headers = qw(genome_file keep_reads reads_in_R1 reads_in_R2 reads_out_R1 
        reads_out_R2);
my $options = {};
GetOptions($options,
    'alignment_list|a=s',
    'testing|t',
    );

my $global_qsub_jobid = 11223300; 
# Only used for testing. Never tries to qsub with this value.

my %reads_jobids = (); # Stores the jobid that will create a pair of reads files. 

# from AssemblyPipeline/Assembly/Qsub.pm - bug
# fix so we 'use' the module above rather than copy the function.
sub get_jobid
{
    my $self = shift;
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

sub parse_list
{
    my $fname = shift;
    my $headers = shift;
    open (FIN, '<', $fname) or die "Error: couldn't open file $fname\n";
    my $records = [];
    <FIN>; # Skip headers
    while (my $line = <FIN>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        if (scalar @fields == scalar @$headers) {
            my $frec = {};
            for (my $i=0; $i<scalar @fields; $i++) {
                my $key = $headers->[$i];
                my $value = $fields[$i];
                $frec->{$key} = $value;
            }
            push (@$records, $frec);
        }
    }
    return $records;   
}

# Given an alignment record, check that the input files exist.
sub input_reads_exist
{
    my $arec = shift;
    my $inr1 = ($arec->{reads_in_R1} ? $arec->{reads_in_R1} : "");
    my $inr2 = ($arec->{reads_in_R2} ? $arec->{reads_in_R2} : "");
    if ($inr1 and $inr2 and -e $inr1 and -e $inr2) {
        return 1;
    }
    return 0;
}

sub get_reads_key
{
    my $reads_R1 = shift;
    my $reads_R2 = shift;
    my $key = $reads_R1 . "_" . $reads_R2;
    return $key;
}

# Set the jobid to wait on that generates the desired reads
sub set_reads_jobid
{
    my $arec = shift;
    my $qsub_jobid = shift;
    my $outr1 = ($arec->{reads_out_R1} ? $arec->{reads_out_R1} : "");
    my $outr2 = ($arec->{reads_out_R2} ? $arec->{reads_out_R2} : "");
    if ($outr1 and $outr2) {
        my $reads_key = get_reads_key($outr1, $outr2);
        $reads_jobids{$reads_key} = $qsub_jobid;
    } else {
        die "Error - could not get the reads output files from the alignment record.\n";
    }
}

sub submit_align_job
{
    my $arec = shift;
    my $parent_qsub_jobid = shift;
    my $alignReads_cmd = $alignReads_path . " --keep_" . $arec->{keep_reads} . "_reads " .
        " --genome " . $arec->{genome_file};
    for my $arg (qw(reads_in_R1 reads_in_R2 reads_out_R1 reads_out_R2)) {
        my $str = " --" . $arg . " " . $arec->{$arg} . " ";
        $alignReads_cmd .= $str;
    }
    print $alignReads_cmd . "\n";
    my $qsub_cmd = $qsub_path . " -N align_reads -hold_jid " . $parent_qsub_jobid . " " . 
            $qsub_script . " '" . $alignReads_cmd . "' ";
    my $qsub_jobid;
    if ($options->{testing}) {
        $qsub_jobid = $global_qsub_jobid;
        $global_qsub_jobid++;
        print "Job above has testing id $qsub_jobid\n";
    } else {
        my $qsub_str = `$qsub_cmd`;
        $qsub_jobid = get_jobid ($qsub_str);
    }
    return $qsub_jobid;
}   

sub run_all
{
    my $alignment_recs = parse_list ($options->{alignment_list}, \@alignment_headers);
    for my $arec (@$alignment_recs) {
        my $parent_qsub_jobid;
        my $qsub_jobid;
        if (input_reads_exist($arec)) {
            $parent_qsub_jobid = 1;
            my $qsub_jobid = submit_align_job($arec, $parent_qsub_jobid);
            set_reads_jobid($arec, $qsub_jobid);
        } else {
            $parent_qsub_jobid = get_parent_jobid($arec);
            my $qsub_jobid = submit_align_job($arec, $parent_qsub_jobid);
            set_reads_jobid($arec, $qsub_jobid);
        }
    }
}

run_all;

