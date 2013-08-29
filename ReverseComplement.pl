#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use Assembly::Utils;
use File::Basename;

my $options = {};
my $fx_rc_bin = "/opt/bio/fastx/bin/fastx_reverse_complement";
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/07_genome_lengths.yml
        yaml_out yaml_files/08_reverse_complement.yml
        verbose 1
        qsub_script qsub_script.sh
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
    #$options->{qsub_opts} = $options->{qsub_opts} . " "; # add some defaults here...
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --run
                --qsub_script
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'run',
        'verbose',
        'qsub_script=s',
        );
    set_default_opts;
    check_opts;
}

sub print_verbose
{
    if ($options->{verbose}) {
        print (@_);
    }
}

sub get_qsub_cmd
{
    my $cmd = shift;
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N revcomp" . 
            " " . $options->{qsub_script} . " '" . $cmd . "'";
    return $qsub_cmd;
}

sub do_revcomp
{
    my ($records, $species, $strain, $sample_type, $trimraw) = @_;
    if (Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type])) {
        my $trimdata = $trimraw . "data";
        my $r1data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R1", $trimdata]);
        my $r2data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R2", $trimdata]);
        if ($r1data and $r2data) {
            my $sample_dir = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "sample_dir"]);
            my $r1out = $sample_dir . "/rev_" . basename($r1data);
            my $r2out = $sample_dir . "/rev_" . basename($r2data);
            my $r1rev = "rev" . $trimdata;
            my $r2rev = "rev" . $trimdata;
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc", "R1"], $r1rev, $r1out);
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc", "R2"], $r2rev, $r2out);
            my $r1_fx_rc_cmd = $fx_rc_bin . " -Q 33 -i " . $r1data . " -o " . $r1out;
            my $r2_fx_rc_cmd = $fx_rc_bin . " -Q 33 -i " . $r2data . " -o " . $r2out;
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc", "R1"], "cmd", $r1_fx_rc_cmd);
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc", "R2"], "cmd", $r1_fx_rc_cmd);
            my $r1_qsub_cmd = get_qsub_cmd($r1_fx_rc_cmd);
            my $r2_qsub_cmd = get_qsub_cmd($r2_fx_rc_cmd);
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc", "R1"], "cmd", $r1_qsub_cmd);
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc", "R2"], "cmd", $r2_qsub_cmd);
            
            print_verbose "Running command:\n" . $r1_qsub_cmd . "\n";
            if ($options->{run} and not -e $r1out) {
                system($r1_qsub_cmd);
            }
            print_verbose "Running command:\n" . $r2_qsub_cmd . "\n";
            if ($options->{run} and not -e $r2out) {
                system($r2_qsub_cmd);
            }
        }
    }
}

sub all_revcomp
{
    my $records = shift;
    for my $species (keys %$records) {
        my $specref = $records->{$species}->{DNA};
        for my $strain (keys %$specref) {
            for my $sample_type (qw(MP MP3 MP8)) {
                for my $trimraw (qw(raw trim)) {
                    my ($r1_cmd, $r2_cmd) = do_revcomp($records, $species, $strain, $sample_type, $trimraw);
                }
            }
        }
    }
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
all_revcomp($records);
DumpFile($options->{yaml_out}, $records);
