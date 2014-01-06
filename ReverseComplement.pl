#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use Assembly::Utils;
use File::Basename;
use Cwd;

my $options = {};
my $fx_rc_bin = "/opt/bio/fastx/bin/fastx_reverse_complement";
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
my $reverse_complement_bin = getcwd . "/ReverseComplement.sh";

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
                --testing
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
        'testing',
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
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N revcomp " . 
            #" -l h=biocomp-0-5 " . 
            $options->{qsub_script} . " '" . $cmd . "'";
    return $qsub_cmd;
}

sub build_cmd_gz
{
    my $frec = shift;
    my $rval = shift;
    my $datafile = shift;
    my $sample_dir = shift;
    my $trimdata = shift;
    
    my $outfile = $sample_dir . "/rev_" . basename($datafile);
    my $rc_cmd = $reverse_complement_bin . " " . $datafile . " " . $outfile;
    my $rc_qsub_cmd = get_qsub_cmd ($rc_cmd);
    my $revdata = "rev" . $trimdata;
    Assembly::Utils::set_check_record($frec, [$rval], $revdata, $outfile);
    Assembly::Utils::set_check_record($frec, [$rval], "cmd", $rc_cmd);
    Assembly::Utils::set_check_record($frec, [$rval], "qsub_cmd", $rc_qsub_cmd);

    print_verbose "Running command:\n" . $rc_qsub_cmd . "\n";
    unless ($options->{testing} or -e $outfile) {
        system($rc_qsub_cmd);
    }

}

sub build_cmd_nogz
{
    my $frec = shift;
    my $rval = shift;
    my $datafile = shift;
    my $sample_dir = shift;
    my $trimdata = shift;
    
    my $outfile = $sample_dir . "/rev_" . basename($datafile) . ".gz";
    my $rc_cmd = $fx_rc_bin . " -Q 33 -i " . $datafile . " -z -o " . $outfile;
    my $rc_qsub_cmd = get_qsub_cmd ($rc_cmd);
    my $revdata = "rev" . $trimdata;
    Assembly::Utils::set_check_record($frec, [$rval], $revdata, $outfile);
    Assembly::Utils::set_check_record($frec, [$rval], "cmd", $rc_cmd);
    Assembly::Utils::set_check_record($frec, [$rval], "qsub_cmd", $rc_qsub_cmd);
    
    print_verbose "Running command:\n" . $rc_qsub_cmd . "\n";
    unless ($options->{testing} or -e $outfile) {
        system($rc_qsub_cmd);
    }
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
            my $frec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "fastx_rc"]);
            if ($r1data =~ /\.gz\s*$/) {
                build_cmd_gz ($frec, "R1", $r1data, $sample_dir, $trimdata);
            } else {
                build_cmd_nogz ($frec, "R1", $r1data, $sample_dir, $trimdata);
            }
            if ($r2data =~ /\.gz\s*$/) {
                build_cmd_gz ($frec, "R2", $r2data, $sample_dir, $trimdata);
            } else {
                build_cmd_nogz ($frec, "R2", $r2data, $sample_dir, $trimdata);
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
