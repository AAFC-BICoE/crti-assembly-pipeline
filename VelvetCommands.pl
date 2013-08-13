#!/usr/bin/env perl
use strict;
use warnings;
use YAML::XS qw(LoadFile DumpFile);
use File::Path;
use Getopt::Long;
use Assembly::Utils;

# From the records, pull just the genomic data
# Read in genome lengths from file
# Write out all the velveth/g commands per sample per kmer

# Ideally: have separate folders within velvet that include
# e.g. the different genome lengths, QC trims etc. that we try?

my $options = {};

# @ kbins is only used by get_kmer_bin function below.
# Assumes we have binaries of form
# velvetg_x and velveth_x for x>0 and x in @kbins.

sub set_default_opts
{
    my %defaults = qw(
            yaml_in yaml_files/08_velvetk.yml
            yaml_out yaml_files/09_velvet_cmds.yml
            min_kmer 21 
            max_kmer 95
            trim 1
            raw 1
            use_velvetk 1
            velvetk_radius 6
            );
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file>
            Optional:
                --trim
                --raw
                --submit
                --verbose
                --min_kmer <value (default 21)>
                --max_kmer <value (default 95)>
                --use_velvetk
                --velvetk_radius <nearby kmers max>
                ";
    }
}


sub gather_opts
{
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'trim',
            'raw',
            'submit',
            'verbose',
            'min_kmer',
            'max_kmer',
            'use_velvetk',
            'velvetk_radius=s',
            );
    set_default_opts;
    check_opts;
}

sub add_velveth_idx
{
	my $sample_ref = shift;
	my $sample_type = shift;
	my $trimraw = shift;
	my $idx = shift;
	
	my $trdata = $trimraw . "data";
	if ($sample_type =~ /MP/) {
		$trdata = "rev" . $trdata;
	}
	my $r1file = Assembly::Utils::get_check_record($sample_ref, [$sample_type, "R1", $trdata]);
	my $r2file = Assembly::Utils::get_check_record($sample_ref, [$sample_type, "R2", $trdata]);

	my $cmd_str = "-shortPaired";
    if ($i) { 
        $cmd_str .= $i+1;
    }
    $cmd_str .= " " . $r1file . " " . $r2file;
	
	return $cmd_str;
}

sub add_velvetg_idx
{
    my $sample_ref = shift;
    my $sample_type = shift;
    my $trimraw = shift;
    my $idx = shift;
    
    my $cmd_str = "-ins_length";
    if ($i) {
        $cmd_str .= $i+1;
    }
    if ($sample_type =~ /PE/) {
        $cmd_str .= " 300";
    } elsif ($sample_type =~ /MP3/) {
        $cmd_str .= " 3000";
    } elsif ($sample_type =~ /MP8/) {
        $cmd_str .= " 8000";
    } else {
        print "Warning: unknown insert length for sample_type " . $sample_type . " using 300\n";
        $cmd_str = " 300";
    }
    if ($sample_type =~ /MP/) {
        $cmd_str .= " -shortMatePaired";
        if ($i) {
            $cmd_str .= ($i+1);
        }
    }
    return $cmd_str;
}

sub base_velveth_cmd
{
    my $kmer_bin = shift;
    my $outdir = shift;
    my $velveth_bin = $velvet_bin_dir . "/velveth_" . $kmer_bin;

    my $velveth_cmd = $velveth_bin . " " . $outdir . " " . $kmer . 
        "  -fastq -create_binary ";
    return $velveth_cmd;
}

sub base_velvetg_cmd
{
    my $kmer_bin = shift;
    my $outdir = shift;
    my $exp_cov = shift;
    my $min_contig_opt = " ";
    my $scaffolding_opt = " -scaffolding yes ";    
    my $velvetg_bin = $velvet_bin_dir . "/velvetg_" . $kmer_bin;
    my $velvetg_cmd = $velvetg_bin . " " . $outdir . " -exp_cov " . $exp_cov .
        " " . $min_contig_opt . $scaffolding_opt . " -amos_file no -cov_cutoff auto";
    return $velvetg_cmd;
}

sub get_velvet_cmds
{
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $kmer = shift;
    my $cov_vars = shift;
    my $assembly_outdir = shift;
    
    my $trdata = $trimraw . "data";
    my $kmer_bin = get_kmer_bin($kmer);
    my $exp_cov = calc_exp_cov($cov_vars, $kmer);
    
    # Create the kmer directory
    my $kmer_dir = $assembly_outdir . "/velvet/" . $trimraw . "/assem_kmer-" . $kmer . "_exp_cov-" . $exp_cov . "_covcutoff-auto";
    unless (-e $kmer_dir) {
        mkpath $kmer_dir;
    }
    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "kmer_dir", $kmer_dir);
    
    my $velveth_cmd = base_velveth_cmd($kmer_bin, $outdir, $kmer);
    my $velvetg_cmd = base_velvetg_cmd($kmer_bin, $outdir);
    my $i = 0;
    for my $sample_type (qw(PE PER MP MP3 MP8)) {
        my $sample_ref = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type]);
        if ($sample_ref) {
            $velveth_cmd .= add_velveth_idx($sample_rec, $sample_type, $trimraw, $i);
            $velvetg_cmd .= add_velvetg_idx($sample_rec, $sample_type, $trimraw, $i);
            $i++;
        }
    }
    if ($i) {
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "velveth_cmd");
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "velvetg_cmd");
    } else {
        $velveth_cmd = '';
        $velvetg_cmd = '';
    }
    return ($velveth_cmd, $velvetg_cmd);
}

sub build_assembly_cmds
{
    my $records = shift;
    for my $species (keys %$records) {
        my $spec_ref = $records->{$species}->{DNA};
        for my $strain (keys %$spec_ref) {
			for my $trimraw (qw(trim raw)) {
				if ($options->{$trimraw}) {
                    my $assembly_dir = Assembly::Velvet::get_assembly_outdir($records, $species, $strain, $trimraw);
                    my $cov_vars = get_coverage_vars($records, $species, $strain, $trimraw);
                    my $kmer_list = Assembly::Velvet::create_kmer_range($records, $species, $strain, $trimraw);
                    for my $kmer (@$kmer_list) {
					    get_velvet_cmds($species, $strain, $trimraw, $kmer, $cov_vars, $assembly_dir);
				}
			}
		}
	}
}

gather_opts;
$records = LoadFile($options->{yaml_in});
build_assembly_cmds;
DumpFile($options->{yaml_out}, $records);
