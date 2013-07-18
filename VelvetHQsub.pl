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

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/06_velvet_cmds.yml
        yaml_out yaml_files/07_velveth_qsub.yml
        qsub_script qsub_array.sh
        vh_batch_dir qsub_files/04_vh_cmds
        trim 1
        raw 0
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
    $options->{qsub_opts} = $options->{qsub_opts} . " -p -500 "; # . -N velvet_hg ";
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --sample_list <ID1,ID2,...,IDX (no spaces)>
                --qsub_array_script <filename>
                --memory
                --trim
                --raw
                --submit
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

sub write_batch_cmds
{
    my $rec = shift;
    my $aref = shift;
    my $vtype = shift; # should be 'vh' or 'vg';  
    my $bdir_key = $vtype . "_batch_dir";
    my @cmd_list = @$aref;
    my $batch_dir = $options->{$bdir_key};
    unless (-e $batch_dir) { mkpath $batch_dir; }
    my $batch_file = $batch_dir . "/" . $rec->{sample} . "_" . $tr . ".sh";
    open (FBATCH, '>', $batch_file) or die "Error: couldn't open file " . $batch_file . "\n";
    print FBATCH join ("\n", @cmd_list) . "\n";
    close (FBATCH);
    return $batch_file;
}

sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]+\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

# Get all the outstanding commands/kmers for the given sample
# together and submit as an array job.
sub get_vh_qsub
{
    my $rec = shift;
    my $aref = shift;
    my $sample = $rec->{sample};
    my @vh_cmd_list = @$aref;
    my $num_cmds = scalar @vh_cmd_list;
    if ($num_cmds > 0) {
        my $batch_filename = write_batch_cmds($rec, $aref, "vh");
        my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N velveth_" . $tr . " -t 1:" . $num_cmds . " " . $options->{qsub_script} . " " . $batch_filename;
        set_check_record($rec, ["velvet", $tr, "velveth"], "qsub_cmd", $qsub_cmd);
        set_check_record($rec, ["velvet", $tr, "velveth"], "cmd_file", $batch_filename);
        if ($options->{verbose}) {
            print "Qsubbing array of " . $num_cmds . " velveth commands for sample " . $rec->{sample} . "\n";
            print $qsub_cmd . "\n";
        }
        print $qsub_cmd . "\n";
        if ($options->{submit}) {
            my $qsub_str = `$qsub_cmd`;
            print $qsub_str;
            my $jobid = get_jobid($qsub_str);
            set_check_record($rec, ["velvet", $tr, "velveth"], "qsub_jobid", $jobid);
        }
    }
}

sub check_vh_sample_cmd
{
    my $rec = shift;
    my $krange = get_kmer_range($rec);
    my $vh_cmds = [];
    for my $kmer (@$krange) {
        my $cmd = check_vh_kmer_cmd($rec, $kmer);
        if ($cmd) {
            push (@$vh_cmds, $cmd);
        }
    }
    get_vh_qsub($rec, $vh_cmds);
}

gather_opts;
$records = LoadFile($options->{yaml_in});
for my $sample (keys %$records) {
    my $rec = $records->{$sample};
    check_vh_sample_cmd($rec);
}
DumpFile($options->{yaml_out}, $records);


