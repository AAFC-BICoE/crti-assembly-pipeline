#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);

my $options = {};
my $records = {};
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
my @sample_list = ();


# output files created by velveth
my @vh_outputs = qw(CnyUnifiedSeq CnyUnifiedSeq.names Log Roadmaps);

sub set_default_opts;
{
    my %defaults = qw(
        vh_batch_dir qsub_files/04_vh_cmds
        trim 1
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
                verbose
                sample_list <ID1,ID2,...,IDX (no spaces)>
                qsub_array_script <filename>
                memory
                trim
                raw
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'qsub_array_script=s',
        'verbose',
        'sample_list=s',
        'memory_requirements|m=s',
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

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $kmin = get_check_record($rec, ["velvet", $tr, "min_kmer"]);
    my $kmax = get_check_record($rec, ["velvet", $tr, "max_kmer"]);
    my $krange = [];
    for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
        push (@$kref, $kmer);
    }
    return $kref;
}

# Check if:
# 1. kmer dir exists
# 2. all output files have been created
# 3. No file size is 0.
sub check_vh_kmer_cmd
{
    my $rec = shift;
    my $kmer = shift;
    my $kdir = get_check_record($rec, ["velvet", $tr, "kmer_dirs", $kmer]);
    my $cmd = get_check_record($rec, ["velveth", $tr, "cmd", $kmer]);
    my $outputs_exist = 1;
    for my $fname (@vh_outputs) {
        my $fpath = $kdir . "/" . $fname;
        if (!-e $fpath or !-s $fpath) { # file nonexistent or empty
            $outputs_exist = 0;
        }
    }
    if ($outputs_exist) {
        return $cmd;
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
    submit_vh_array($vh_cmds);
}

sub get_unfinished_vh_cmds
{
    for my $sample (@sample_list) {
        my $rec = $records->{$sample};
        my $cmd = check_vh_cmd($rec, $

}



sub get_vh_qsub
{
    my $rec = shift;
    my $aref = shift;
    my $sample = $rec->{sample};
    my @vh_cmd_list = @$aref;
    my $num_cmds = scalar @vh_cmd_list;
    my $batch_filename = write_batch_cmds($rec, $aref, "vh");
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N " . $tr . "_velveth -t 1:" . $num_cmds . " " . $options->{qsub_script} . " " . $batch_filename;
    push (@vh_qsub_list, $qsub_cmd);
    $rec->{velveth}->{$tr}->{qsub_cmd} = $qsub_cmd;
    $rec->{velveth}->{$tr}->{cmd_file} = $batch_filename;
}

sub get_vg_qsub
{
    my $rec = shift;
    my $aref = shift;
    my $sample = $rec->{sample};
    my @vg_cmd_list = @$aref;
    my $num_cmds = scalar @vg_cmd_list;
    my $batch_filename = write_batch_cmds($rec, $aref, "vg");
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N " . $tr . "_velvetg -t 1:" . $num_cmds . " " . $options->{qsub_script} . " " . $batch_filename;
    push (@vg_qsub_list, $qsub_cmd);
    $rec->{velvetg}->{$tr}->{qsub_cmd} = $qsub_cmd;
    $rec->{velvetg}->{$tr}->{cmd_file} = $batch_filename;
}

sub vh_write_qsub_batch
{
    my $outfile = $options->{vh_qsub_batch};
    open (FQS, '>', $outfile) or die "Error: couldn't open file $outfile.\n";
    print FQS join("\n", @vh_qsub_list) . "\n";
    close (FQS);
}

sub vg_write_qsub_batch
{
    my $outfile = $options->{vg_qsub_batch};
    open (FQS, '>', $outfile) or die "Error: couldn't open file $outfile.\n";
    print FQS join("\n", @vg_qsub_list) . "\n";
    close (FQS);    
}

# For each sample, look at the
# kmers that have successfully
# completed velvetg. Submit the rest.

my $records = LoadFile($options->{yaml_in});

vh_write_qsub_batch;
vg_write_qsub_batch;