#!/usr/bin/env perl
# Assembly::Velvet.pm
# Author: Jeff Cullis
# Date: August 11, 2013

package Assembly::Qsub;

use strict;
use warnings;
use Math::Random;

sub build_slot_info
{
    my $self = shift;
    my $slot_info = [];
    my $host_numbers = [];
    #for (my $i=0; $i<8; $i++) {
    for my $i (3,5,6,7,8) {
        push (@$host_numbers, $i) for (1..2);
    }
    push (@$host_numbers, 9) for (1..8);
    push (@$host_numbers, 10) for (1..8);
    $self->{num_slots} = scalar @$host_numbers;
    my $job_ids = [];
    for (my $i=0; $i<$self->{num_slots}; $i++) {
        my $rec = {
                hostname => "biocomp-0-" . $host_numbers->[$i],
                hold_jid => 93975,
                };
        push (@$slot_info, $rec);
    }
    $self->{slot_info} = $slot_info;
}

sub new
{
    my $class = shift;
    my $qsub_script = shift;
    my $testing = shift;
    my $submit_max = shift;
    my $verbose = shift;
    my $self = {
        qsub_script => $qsub_script,
        testing => $testing,
        submit_max => $submit_max,
        verbose => $verbose,
        counter => 0,
        total_submitted => 0,
        qsub_bin => "/opt/gridengine/bin/lx26-amd64/qsub",
        };
    build_slot_info($self);
    bless $self, $class;
    return $self;
}

sub print_verbose
{
    my $self = shift;
    if ($self->{verbose}) {
        print (@_);
    }
}

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

sub inc_counter
{
    my $self = shift;
    $self->{counter} = ($self->{counter} + 1) % $self->{num_slots};
}

sub inc_total
{
    my $self = shift;
    $self->{total_submitted} = $self->{total_submitted} + 1;
}

sub submit_velvet
{
    my $self = shift;
    my $velvet_cmd = shift;
    my $velvet_type = shift; # type=velveth or velvetg
    my $count = $self->{counter};
    my $hostname = $self->{slot_info}->[$count]->{hostname};
    my $hold_jid = $self->{slot_info}->[$count]->{hold_jid};
    my $hold_param_str = ($hold_jid ? " -hold_jid " . $hold_jid . " " : "");
    my $qsub_cmd = $self->{qsub_bin} . " -N " . $velvet_type . " -l h=" . $hostname . 
            " " . $hold_param_str . $self->{qsub_script} . " '" . $velvet_cmd . "'";
    my ($jobid, $qsub_str) = ('',''); 
    #print $qsub_cmd . "\n";
    unless ($self->{testing}) {
        my $exceed_max = 0;
        if ($self->{submit_max} and $self->{total_submitted} > $self->{submit_max}) {
            $exceed_max = 1;
        }
        unless ($exceed_max) {
            $qsub_str = `$qsub_cmd`;
            print_verbose $self, $qsub_str . "\n";
            $jobid = get_jobid($self, $qsub_str);
            print_verbose $self, "Got jobid " . $jobid . "\n";
            inc_total($self);
        }
    } elsif ($self->{testing}) {
        # testing jobid
        $jobid = random_uniform_integer(1, 100000, 190000);
    }
    $self->{slot_info}->[$count]->{hold_jid} = $jobid; 
    my $cmd_id = [$qsub_cmd, $jobid];
    return $cmd_id;    
}

sub submit_vh
{
    my $self = shift;
    my $velveth_cmd = shift;
    my $trimraw = shift;
    my $velvet_type = "velveth_" . $trimraw;
    my $cmd_id = submit_velvet($self, $velveth_cmd, $velvet_type);
    inc_counter($self);
    return $cmd_id
}

sub submit_vg
{
    my $self = shift;
    my $velvetg_cmd = shift;
    my $trimraw = shift;
    my $velvet_type = "velvetg_" . $trimraw;
    my $cmd_id = submit_velvet($self, $velvetg_cmd, $velvet_type);
    inc_counter($self);
    return $cmd_id;
}

sub submit_vhg
{
    my $self = shift;
    my $velveth_cmd = shift;
    my $velvetg_cmd = shift;
    my $trimraw = shift;
    #print "Submitting cmds $velveth_cmd \n\n $velvetg_cmd \n\n";
    my $vh_type = "velveth_" . $trimraw;
    my $vg_type = "velvetg_" . $trimraw;
    my $vh_cmd_id = submit_velvet($self, $velveth_cmd, $vh_type);
    my $vg_cmd_id = submit_velvet($self, $velvetg_cmd, $vg_type);
    inc_counter($self);
    my $merged_cmd_id = [@$vh_cmd_id, @$vg_cmd_id];
    return $merged_cmd_id;
}

return 1;