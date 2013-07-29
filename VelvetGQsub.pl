#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
# output files created by velveth/g
my $vh_outfiles = [qw(CnyUnifiedSeq CnyUnifiedSeq.names Log Roadmaps)];
my $vg_outfiles = [qw(Graph2 LastGraph PreGraph stats.txt)];
my @job_table_headers = ("Sample", "Trim/raw", "kmer", "Qsub Job ID", "Command");

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/10_velveth_qsub_host.yml
        yaml_out yaml_files/11_velvetg_qsub.yml
        job_table_in input_data/VGQsubJobIDs.tab
        job_table_out output_files/VGQsubJobIDsOut.tab
        qsub_script qsub_script.sh
        trim 1
        raw 1
        verbose 0
        submit 0
        submit_max 0
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    $options->{qsub_opts} = $options->{qsub_opts} . " -p -400 "; # . -N velvet_hg ";
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

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $trimraw = shift;
    my $kmin = get_check_record($rec, ["velvet", $trimraw, "min_kmer"]);
    my $kmax = get_check_record($rec, ["velvet", $trimraw, "max_kmer"]);
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

# Return the velvetg command to run, unless the three following conditions are met:
# 1. kmer dir exists
# 2. all output files have been created
# 3. No file size is 0.
sub check_vg_cmd
{
    my $rec = shift;
    my $trimraw = shift;
    my $kmer = shift;
    my $kdir = get_check_record($rec, ["velvet", $trimraw, "kmer", $kmer, "kmer_dir"]);
    my $cmd = get_check_record($rec, ["velvet", $trimraw, "kmer", $kmer, "velvetg_cmd"]);
    
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
        print_verbose "Not running for this kmer - found valid velvetg output files for sample " . $rec->{sample} . " trim/raw=" . $trimraw . " kmer=" . $kmer . "\n";
        $cmd = '';
    }
    return $cmd;
}

# Returns a list of sample/trimraw/kmer/velveth command combos where
# no valid output files were found - i.e. velveth needs to be run.
sub get_all_cmds
{
    my $records = shift;
    my $sample_list = shift;
    my $qsub_param_list = [];
    for my $sample (@$sample_list) {
        my $rec = $records->{$sample};
        for my $trimraw (qw(trim raw)) {
            if ($options->{$trimraw}) {
                my $krange = get_kmer_range($rec, $trimraw);
                for my $kmer (@$krange) {
                    my $cmd = check_vg_cmd($rec, $trimraw, $kmer);
                    if ($cmd) {
                        push (@$qsub_param_list, [$sample, $trimraw, $kmer, $cmd]);
                    }
                }
            }
        }
    }
    return $qsub_param_list;
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
    my $trimraw = shift;
    my $host = shift;
    my $cmd = shift;
    my $hold_jid = shift;
    my $hold_param_str = ($hold_jid ? " -hold_jid " . $hold_jid . " " : "");
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N velvetg_" . $trimraw . 
            " -l h=" . $host . " " . $hold_param_str . $options->{qsub_script} . " '" . $cmd . "'";
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
    my $records = shift;
    my $qsub_param_list = shift;
    my $jobid_list = shift;
    my $slots = get_processor_slots;
    my $num_slots = scalar @$slots;
    my @holds = ();
    push (@holds, 81213) for (1..$num_slots);
    my $num_to_submit = scalar @$qsub_param_list;
    if ($options->{submit_max} and $options->{submit_max} < $num_to_submit) {
        $num_to_submit = $options->{submit_max};
    }
    my $new_qsub_ids = [];
    for (my $i=0; $i<$num_to_submit; $i++) {
        my ($sample, $trimraw, $kmer, $cmd) = @{$qsub_param_list->[$i]};
        my $j = $i % $num_slots;
        my $host = "biocomp-0-" . $slots->[$j];
        my ($qsub_cmd, $jobid) = do_host_qsub($trimraw, $host, $cmd, $holds[$j]);
        $holds[$j] = $jobid;
        set_check_record($records, [$sample, "velvet", $trimraw, "kmer", $kmer], "velvetg_qsub_cmd", $qsub_cmd);
        set_check_record($records, [$sample, "velvet", $trimraw, "kmer", $kmer], "velvetg_qsub_jobid", $jobid);
        push (@$jobid_list, [$sample, $trimraw, $kmer, $jobid, $cmd]);
    }
    return $jobid_list;
}

sub add_old_jobids
{
    my $records = shift;
    my $fname = ($options->{job_table_in} ? $options->{job_table_in} : '');
    my $jobid_list = [];
    if ($fname and -e $fname) {
        open (JID_IN, '<', $fname) or die "Error: couldn't open input job table $fname.\n";
        while (my $line = <JID_IN>) {
            chomp $line;
            my @fields = split (/\s+/, $line);
            if (scalar @fields == 4) {
                my ($sample, $trimraw, $kmer, $old_jobid, $old_cmd) = @fields;
                set_check_record($records, [$sample, "velvet", $trimraw, "kmer", $kmer], "velvetg_qsub_cmd", $old_cmd);
                set_check_record($records, [$sample, "velvet", $trimraw, "kmer", $kmer], "velvetg_qsub_jobid", $old_jobid);
                push (@$jobid_list, [$sample, $trimraw, $kmer, $old_jobid, $old_cmd]);
            }
        }
        close (JID_IN);
    }
    return $jobid_list;
}
                    
sub write_all_jobids
{
    my $jobid_list = shift;
    my $fname = ($options->{job_table_out} ? $options->{job_table_iout} : '');
    if ($fname) {
        open (JID_OUT, '>', $fname) or die "Error: couldn't open output job table $fname.\n";
        print JID_OUT join("\t", @job_table_headers) . "\n";
        for my $arr (@$jobid_list) {
            print JID_OUT join("\t", @$arr) . "\n";
        }
        close (JID_OUT);
    }
}

sub run_all
{
    gather_opts;
    my $records = LoadFile($options->{yaml_in});
    my $jobid_list = add_old_jobids($records);
    my $sample_list = get_sample_list($records);
    my $qsub_param_list = get_all_cmds($records, $sample_list);
    qsub_cmds($records, $qsub_param_list, $jobid_list);
    write_all_jobids($jobid_list);
    DumpFile($options->{yaml_out}, $records);
}

run_all;
