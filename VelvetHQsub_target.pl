#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my $records = {};
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
my @sample_list = ();
my $tr;
my $trdata;


# output files created by velveth
my @vh_outputs = qw(CnyUnifiedSeq CnyUnifiedSeq.names Log Roadmaps);

my $vh_cmds = {};

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/06_velvet_cmds.yml
        yaml_out yaml_files/07_velveth_qsub_target.yml
        qsub_script qsub_script.sh
        trim 0
        raw 1
        verbose 0
        submit 0
        submit_max 0
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
    $options->{qsub_opts} = $options->{qsub_opts} . " -p -500 "; # . -N velvet_hg ";
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --sample_list <ID1,ID2,...,IDX (no spaces)>
                --qsub_script <filename>
                --memory
                --trim
                --raw
                --submit
                --submit_max <max # jobs to qsub>
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'qsub_script=s',
        'verbose',
        'sample_list=s',
        'trim',
        'raw',
        'submit',
        'submit_max=s',
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

# Get an array of all the kmer values
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

# Check if:
# 1. kmer dir exists
# 2. all output files have been created
# 3. No file size is 0.
sub check_vh_kmer_cmd
{
    my $rec = shift;
    my $kmer = shift;
    #my $kdir = get_check_record($rec, ["velvet", $tr, "kmer_dirs", $kmer]);
    #my $cmd = get_check_record($rec, ["velveth", $tr, "cmd", $kmer]);
    my $kdir = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "kmer_dir"]);
    my $cmd = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "velveth_cmd"]);
    
    my $outputs_exist = 1;
    for my $fname (@vh_outputs) {
        my $fpath = $kdir . "/" . $fname;
        if (!-e $fpath or !-s $fpath) { # file nonexistent or empty
            $outputs_exist = 0;
        }
    }
    if (!$outputs_exist) {
        if ($options->{verbose}) {
            print "sample: " . $rec->{sample} . ":\n";
            print "Found missing/empty output files, running the following command:\n$cmd\n";
        }
        # delete any files that do exist first?
        return $cmd;
    }
    if ($options->{verbose}) {
        print "Not running for this kmer - found valid output files for sample " . $rec->{sample} . " trim/raw=" . $tr . " kmer=" . $kmer . "\n";
    }
    return '';
}

sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

sub get_vh_sample_cmds
{
    my $rec = shift;
    my $sample = $rec->{sample};
    my $krange = get_kmer_range($rec);
    #my $vh_cmds = [];
    for my $kmer (@$krange) {
        my $cmd = check_vh_kmer_cmd($rec, $kmer);
        if ($cmd) {
            $vh_cmds->{$sample}->{$kmer} = $cmd;
        }
    }

}

sub get_processor_slots
{
    my $slots = [];
    for (my $i=0; $i<4; $i++) {
        push (@$slots, $i+4) for (1..2);
    }
    push (@$slots, 9) for (1..4);
    push (@$slots, 10) for (1..8);
    return $slots;
}

sub do_host_qsub
{
    my $host = shift;
    my $cmd = shift;
    my $hold_jid = shift;
    my $hold_param_str = ($hold_jid ? " -hold_jid " . $hold_jid . " " : "");
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N velveth_" . $tr . 
            " -l h=" . $host . " -p -400 " . $hold_param_str . $options->{qsub_script} . " '" . $cmd . "'";
    my ($jobid, $qsub_str) = ('','');
    print $qsub_cmd . "\n";
    if ($options->{submit}) {
        $qsub_str = `$qsub_cmd`;
        print $qsub_str . "\n";
        $jobid = get_jobid($qsub_str);
        print "Got jobid " . $jobid . "\n";
    }
    return ($qsub_cmd, $jobid);
}

sub qsub_cmds
{
    my $slots = shift;
    my $num_slots = scalar @$slots;
    my @holds = ();
    push (@holds, 72429) for (1..$num_slots);
    
    my $i = 0;
    for my $sample (keys %$vh_cmds) {
        my $vh_sample_cmds = $vh_cmds->{$sample};
        for my $kmer (keys %$vh_sample_cmds) {
            my $j = $i % $num_slots;
            my $cmd = $vh_sample_cmds->{$kmer};
            my $host = "biocomp-0-" . $slots->[$j];
            my ($qsub_cmd, $jobid) = do_host_qsub($host, $cmd, $holds[$j]);
            $holds[$j] = $jobid;
            set_check_record($records, [$sample, "velvet", $tr, "kmer", $kmer], "velveth_qsub_cmd", $qsub_cmd);
            set_check_record($records, [$sample, "velvet", $tr, "kmer", $kmer], "velveth_qsub_jobid", $jobid);
            $i++;
            last if ($options->{submit_max} and $i >= $options->{submit_max});
        }
        last if ($options->{submit_max} and $i >= $options->{submit_max});
    }
}

gather_opts;
$records = LoadFile($options->{yaml_in});
for my $sample (keys %$records) {
    my $rec = $records->{$sample};
    get_vh_sample_cmds($rec);
}
my $slots = get_processor_slots;
qsub_cmds($slots);
DumpFile($options->{yaml_out}, $records);


