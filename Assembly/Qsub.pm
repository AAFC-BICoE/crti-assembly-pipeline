#!/usr/bin/env perl
# Assembly::Velvet.pm
# Author: Jeff Cullis
# Date: August 11, 2013

package Assembly::Qsub;

use strict;
use warnings;
use Assembly::Utils;

sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
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
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N velveth_" . $trimraw . 
            " -l h=" . $host . "  " . $hold_param_str . $options->{qsub_script} . " '" . $cmd . "'";
    my ($jobid, $qsub_str) = ('','');
    print $qsub_cmd . "\n";
    if ($options->{submit}) {
        $qsub_str = `$qsub_cmd`;
        print_verbose $qsub_str . "\n";
        $jobid = get_jobid($qsub_str);
        print_verbose "Got jobid " . $jobid . "\n";
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
    push (@holds, 0) for (1..$num_slots);
    my $num_to_submit = scalar @$qsub_param_list;
    if ($options->{submit_max} and $options->{submit_max} < $num_to_submit) {
        $num_to_submit = $options->{submit_max};
    }
    for (my $i=0; $i<$num_to_submit; $i++) {
        my ($sample, $trimraw, $kmer, $cmd) = @{$qsub_param_list->[$i]};
        my $j = $i % $num_slots;
        my $host = "biocomp-0-" . $slots->[$j];
        my ($qsub_cmd, $jobid) = do_host_qsub($trimraw, $host, $cmd, $holds[$j]);
        $holds[$j] = $jobid;
        set_check_record($records, [$sample, "velvet", $trimraw, "kmer", $kmer], "velveth_qsub_cmd", $qsub_cmd);
        set_check_record($records, [$sample, "velvet", $trimraw, "kmer", $kmer], "velveth_qsub_jobid", $jobid);
        push (@$jobid_list, [$sample, $trimraw, $kmer, $jobid, $cmd]);
    }
}

return 1;