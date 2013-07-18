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

# output files created by velveth/g
my $vh_outfiles = [qw(CnyUnifiedSeq CnyUnifiedSeq.names Log Roadmaps)];
my $vg_outfiles = [qw(Graph2 LastGraph PreGraph stats.txt)];

my $vg_cmds = {};


sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/07_velveth_qsub.yml
        yaml_out yaml_files/08_velvetg_qsub.yml
        qsub_script qsub_script.sh
        vg_batch_dir qsub_files/04_vg_cmds
        trim 1
        verbose 0
        submit 0
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
                --qsub_script <filename>
                --memory
                --trim
                --raw
                --submit
                --submit_max <max # qsubs to perform>
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
        'memory_requirements|m=s',
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
    

sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

sub get_sample_list
{
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

# check that a set of filenames exist in a given directory.
sub outfiles_exist
{
    my $kdir = shift;
    my $filenames = shift;
    my $outputs_exist = 1;
    for my $fname (@$filenames) {
        my $fpath = $kdir . "/" . $fname;
        if (!-e $fpath or !-s $fpath) { # file nonexistent or empty
            $outputs_exist = 0;
        }
    }
    return $outputs_exist;
}  

# Check if:
# 1. kmer dir exists
# 2. all output files have been created
# 3. No file size is 0.
sub get_vg_kmer_cmd
{
    my $rec = shift;
    my $kmer = shift;
    my $kdir = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "kmer_dir"]);
    my $cmd = get_check_record($rec, ["velvet", $tr, "kmer", $kmer, "velvetg_cmd"]);
    
    my $vh_files_exist = outfiles_exist($kdir, $vh_outfiles);
    my $vg_files_exist = outfiles_exist($kdir, $vg_outfiles);

    if ($vh_files_exist and !$vg_files_exist) {
        if ($cmd) {
            print_verbose "Found vh files but not vg files. Running the following command:\n$cmd\n";
        } else {
            print_verbose "No command found to run for kmer " . $kmer . " sample " . $rec->{sample} . "\n";
            $cmd = '';
        }
    } elsif (!$vh_files_exist) {
        print_verbose "Need to run velveth before it's possible to run velvetg for folder $kdir\n";
        $cmd = '';
    } elsif ($vh_files_exist) {
        print_verbose "Not running for this kmer - found valid velvetg output files for sample " . $rec->{sample} . " trim/raw=" . $tr . " kmer=" . $kmer . "\n";
        $cmd = '';
    }
    return $cmd;
}

sub get_vg_sample_cmds
{
    my $rec = shift;
    my $sample = $rec->{sample};
    my $krange = get_kmer_range($rec);
    for my $kmer (@$krange) {
        my $cmd = get_vg_kmer_cmd($rec, $kmer);
        if ($cmd) {
            #push (@vg_cmds, $cmd);
            $vg_cmds->{$sample}->{$kmer} = $cmd;
        }
    }
}

sub get_processor_slots
{
    my $slots = [];
    for (my $i=0; $i<4; $i++) {
        push (@$slots, $i+4);
        push (@$slots, $i+4);
    }
    push (@$slots, 9) for (1..8);
    push (@$slots, 10) for (1..8);
    return $slots;
}

sub do_host_qsub
{
    my $host = shift;
    my $cmd = shift;
    my $hold_jid = shift;
    my $hold_param_str = ($hold_jid ? " -hold_jid " . $hold_jid . " " : "");
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N velvetg_" . $tr . 
            " -l h=" . $host . " -p -400 " . $hold_param_str . $options->{qsub_script} . " '" . $cmd . "'";
    my ($jobid, $qsub_str) = ('','');
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
    push (@holds, 0) for (1..$num_slots);
    
    my $i = 0;
    for my $sample (keys %$vg_cmds) {
        my $vg_sample_cmds = $vg_cmds->{$sample};
        for my $kmer (keys %$vg_sample_cmds) {
            my $j = $i % $num_slots;
            my $cmd = $vg_sample_cmds->{$kmer};
            my $host = "biocomp-0-" . $slots->[$j];
            my ($qsub_cmd, $jobid) = do_host_qsub($host, $cmd, $holds[$j]);
            $holds[$j] = $jobid;
            set_check_record($records, [$sample, "velvet", $tr, "kmer", $kmer], "velvetg_qsub_cmd", $qsub_cmd);
            set_check_record($records, [$sample, "velvet", $tr, "kmer", $kmer], "velvetg_qsub_jobid", $jobid);
            $i++;
            last if ($options->{submit_max} and $i >= $options->{submit_max});
        }
        last if ($options->{submit_max} and $i >= $options->{submit_max});
    }
}
        

gather_opts;
$records = LoadFile($options->{yaml_in});
my $sample_list = get_sample_list;
for my $sample (@$sample_list) {
    my $rec = $records->{$sample};
    get_vg_sample_cmds($rec);
}
my $slots = get_processor_slots;
qsub_cmds($slots);
DumpFile($options->{yaml_out}, $records);


