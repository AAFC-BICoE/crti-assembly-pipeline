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

sub get_coverage_vars
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $sample = shift;
    my $trimraw = shift;
    my $trdata = $trimraw . "data";
    print "Getting cov vars for " . $sample . "\n";
    my $var = {};
    $var->{R1_nreads} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R1", $trdata, "num_reads"]);
    $var->{R1_readlen} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R1", $trdata, "read_length"]);
    $var->{R2_nreads} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R2", $trdata, "num_reads"]);
    $var->{R2_readlen} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "samples", $sample, "data_stats", "R2", $trdata, "read_length"]);
    $var->{genome_length} = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "related_genome_length", "RG_Est_Genome_Length"]);
    my $pass = 1;
    for my $key (keys %$var) {
        if ($var->{$key} !~ /^\s*\d+\s*$/) {
            print "Got bad var " . $var->{$key} . " at key " . $key . "\n";
            $pass = 0;
        }
    }
    if ($pass) {
        $var->{total_coverage} = ($var->{R1_nreads} * $var->{R1_readlen} + $var->{R2_nreads} * $var->{R2_readlen}) / $var->{genome_length};
        $var->{avg_readlen} = ($var->{R1_readlen} + $var->{R2_readlen}) / 2;
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "total_coverage", $var->{total_coverage});
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "average_read_length", $var->{avg_readlen});
        return $var;
    } else {
        return '';
    }
}

sub calc_exp_cov
{
    my $var = shift;
    my $kmer = shift;
    my $exp_cov_float = $var->{total_coverage} * ($var->{avg_readlen} - $kmer + 1) / $var->{avg_readlen};
    my $exp_cov_rounded_int = int($exp_cov_float + 0.5);
    return $exp_cov_rounded_int;
} 




sub get_velveth_cmd
{
    my $rec = shift;
    my $trimraw = shift;
    my $kmer = shift;
    my $kmer_bin = shift;
    my $outdir = shift;
    my $trdata = $trimraw . "data";
    my $velveth_bin = $velvet_bin_dir . "/velveth_" . $kmer_bin;
    my $r1_file = $rec->{R1}->{$trdata};
    my $r2_file = $rec->{R2}->{$trdata};
    my $velveth_cmd = $velveth_bin . " " . $outdir . " " . $kmer . 
        "  -fastq -shortPaired -create_binary -separate " . $r1_file . " " . $r2_file;
    Assembly::Utils::set_check_record($rec, ["velvet", $trimraw, "kmer", $kmer], "velveth_cmd", $velveth_cmd);
    return $velveth_cmd;
}

sub get_velvetg_cmd
{
    my $rec = shift;
    my $trimraw = shift;
    my $exp_cov = shift;
    my $kmer = shift;
    my $kmer_bin = shift;
    my $working_dir = shift;
    my $min_contig_opt = " ";
    my $scaffolding_opt = " -scaffolding yes ";
    my $velvetg_bin = $velvet_bin_dir . "/velvetg_" . $kmer_bin;
    my $velvetg_cmd = $velvetg_bin . " " . $working_dir . " " .
        "-ins_length 300 -exp_cov " . $exp_cov . " " . $min_contig_opt . 
        $scaffolding_opt . " -amos_file no -cov_cutoff auto";
    Assembly::Utils::set_check_record($rec, ["velvet", $trimraw, "kmer", $kmer], "velvetg_cmd", $velvetg_cmd);
    return $velvetg_cmd;
}

sub get_velvet_cmds
{
    my $rec = shift;
    my $trimraw = shift;
    my $cov_vars = shift;
    my $assembly_outdir = shift;
    my $kmer = shift;
    my $kmer_bin = get_kmer_bin($kmer);
    my $exp_cov = calc_exp_cov($cov_vars, $kmer);
    my $kmer_dir = get_velvet_kdir($rec, $trimraw, $assembly_outdir, $kmer, $exp_cov);
    my $vh_cmd = get_velveth_cmd($rec, $trimraw, $kmer, $kmer_bin, $kmer_dir);
    my $vg_cmd = get_velvetg_cmd($rec, $trimraw, $exp_cov, $kmer, $kmer_bin, $kmer_dir);
    return ($vh_cmd, $vg_cmd);
}

sub add_velveth_files
{
	my $sample_ref = shift;
	my $sample_type = shift;
	my $trimraw = shift;
	my $trdata = $trimraw . "data";
	if ($sample_type =~ /MP/) {
		$trdata = "rev" . $trdata;
	}
	my $r1file = Assembly::Utils::get_check_record($sample_ref, [$sample_type, "R1", $trdata]);
	my $r2file = Assembly::Utils::get_check_record($sample_ref, [$sample_type, "R2", $trdata]);
	return $r1file . " " . $r2file;
}

sub build_strain_cmds
{
	my $records = shift;
	my $species = shift;
	my $strain = shift;
	my $trimraw = shift;
    my $kmer_range = Assembly::Velvet::get_kmer_range($records, $species, $strain, $trimraw);
    for my $kmer (@$kmer_range) {
        my $velveth_cmd = $velveth_bin . " -fastq ";
        my $i = 0;
        for my $sample_type (qw(PE PER MP MP3 MP8)) {
            my $sample_ref = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type]);
            if ($sample_ref) {
                $velveth_cmd .= "-shortPaired";
                if ($i) { $velveth_cmd .= $i+1; }
                $i++;
                $velveth_cmd .= " -separate ";
                $velveth_cmd .= add_velveth_files($sample_ref, $sample_type, $trimraw);
            }
		}
        if ($i) {
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "velveth_cmd");
            push (vh_cmds
        }
	}
	 
}
	

sub build_assembly_cmds
{
    my $records = shift;
    for my $species (keys %$records) {
        my $spec_ref = $records->{$species}->{DNA};
        for my $strain (keys %$spec_ref) {
			for my $trimraw (qw(trim raw)) {
				if ($options->{$trimraw}) {
					build_strain_cmds($records, $species, $strain, $trimraw);
				}
			}
		}
	}
}
            my $samp_ref = $spec_ref->{$strain}->{samples};
            my @sample_list = keys %$samp_ref;
            if (scalar @sample_list == 1)  {
                for my $trimraw (qw(trim raw)) {
                    if ($options->{$trimraw}) {
                        #my $rec = ($records->{$sample} ? $records->{$sample} : '');
                        #unless (defined($rec)) { die "Couldn't get a record from yaml file for sample $sample.\n"; }
                        my $cov_vars = get_coverage_vars($records, $species, $strain, $sample, $trimraw);
                        if ($cov_vars) {
                            my $assembly_outdir = get_assembly_outdir($rec, $trimraw);
                            my $kmer_range = get_kmer_range($rec, $trimraw);
                            for my $kmer (@$kmer_range) {
                                my ($vh_cmd, $vg_cmd) = get_velvet_cmds($rec, $trimraw, $cov_vars, $assembly_outdir, $kmer);
                            }
                        } else {
                            print "Warning: Couldn't parse coverage vars from input sample $sample.\n";
                        }
                    }
                }
            }
        }
    }
}

gather_opts;
$records = LoadFile($options->{yaml_in});
build_assembly_cmds;
DumpFile($options->{yaml_out}, $records);
