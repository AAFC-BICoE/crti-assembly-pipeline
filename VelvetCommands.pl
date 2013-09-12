#!/usr/bin/env perl
use strict;
use warnings;
use YAML::XS qw(LoadFile DumpFile);
use File::Path;
use Getopt::Long;
use Assembly::Utils;
use Assembly::Velvet;
use Assembly::Qsub;

# From the records, pull just the genomic data
# Read in genome lengths from file
# Write out all the velveth/g commands per sample per kmer

# Ideally: have separate folders within velvet that include
# e.g. the different genome lengths, QC trims etc. that we try?

my $options = {};
my $velvet_bin_dir = "/opt/bio/velvet";
#my $vh_outfiles = [qw(CnyUnifiedSeq CnyUnifiedSeq.names Log Roadmaps)];
#my $vg_outfiles = [qw(Graph2 LastGraph PreGraph stats.txt contigs.fa)];
# Modify the above - we plan on deleting all but Log, Roadmaps, stats.txt, contigs.fa
# in order to save space. Don't want the script to re-run where we've deleted files.
my $vh_outfiles = [qw(Log Roadmaps)];
my $vg_outfiles = [qw(stats.txt contigs.fa)];


# @ kbins is only used by get_kmer_bin function below.
# Assumes we have binaries of form
# velvetg_x and velveth_x for x>0 and x in @kbins.

sub set_default_opts
{
    my %defaults = qw(
            yaml_in yaml_files/10_velvetk.yml
            yaml_out yaml_files/11_velvet_cmds.yml
            qsub_script qsub_script.sh
            submit 0
            submit_max 0
            min_kmer 75 
            max_kmer 75
            trim 1
            raw 1
            use_velvetk 1
            velvetk_radius 6
            specimen_dir ../../processing_test2
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
                --qsub_script <filename>
                --submit
                --submit_max <#>
                --verbose
                --min_kmer <value (default 21)>
                --max_kmer <value (default 95)>
                --use_velvetk
                --velvetk_radius <nearby kmers max>
                --specimen_dir <dirname>
                ";
    }
}


sub gather_opts
{
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'qsub_script=s',
            'trim',
            'raw',
            'submit',
            'submit_max=s',
            'verbose',
            'min_kmer',
            'max_kmer',
            'use_velvetk',
            'velvetk_radius=s',
            'specimen_dir=s',
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

sub add_velveth_idx
{
	my $sample_ref = shift;
	my $sample_type = shift;
	my $trimraw = shift;
	my $idx = shift;
	
	my $trdata = $trimraw . "data";
	my ($r1file, $r2file) = ('','');
	if ($sample_type =~ /MP/) {
		my $revtrdata = "rev" . $trdata;
	    $r1file = Assembly::Utils::get_check_record($sample_ref, ["fastx_rc", "R1", $revtrdata]);
	    $r2file = Assembly::Utils::get_check_record($sample_ref, ["fastx_rc", "R2", $revtrdata]);
	} else {
	    $r1file = Assembly::Utils::get_check_record($sample_ref, ["R1", $trdata]);
	    $r2file = Assembly::Utils::get_check_record($sample_ref, ["R2", $trdata]);
	}
    #unless ($r1file and $r2file) { print "Missing data for r1file and r2file sample type is $sample_type\n"; }
	my $cmd_str = '';
	if ($r1file and $r2file) {
        $cmd_str .= "-shortPaired";
        if ($idx) { 
            $cmd_str .= $idx+1;
        }
        $cmd_str .= " -separate -fastq " . $r1file . " " . $r2file. " ";
    }
	
	return $cmd_str;
}

sub add_velvetg_idx
{
    my $sample_ref = shift;
    my $sample_type = shift;
    my $trimraw = shift;
    my $idx = shift;
    
    my $cmd_str = " -ins_length";
    if ($idx) {
        $cmd_str .= $idx+1;
    }
    if ($sample_type =~ /PE/) {
        $cmd_str .= " 300 ";
    } elsif ($sample_type =~ /MP3/) {
        $cmd_str .= " 3000 ";
    } elsif ($sample_type =~ /MP8/) {
        $cmd_str .= " 8000 ";
    } else {
        print "Warning: unknown insert length for sample_type " . $sample_type . " using 300\n";
        $cmd_str = " 300 ";
    }
    # Assume we don't need the below code since default is no and 
    # mate-paired reads are not contaminated with paired-end reads.
    #if ($sample_type =~ /MP/) {
    #   $cmd_str .= " -shortMatePaired";
    #   if ($idx) {
    #       $cmd_str .= ($idx+1);
    #   }
    #   $cmd_str .= " no ";
    #}
    return $cmd_str;
}

sub base_velveth_cmd
{
    my $kmer = shift;
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
        " " . $min_contig_opt . $scaffolding_opt . " -amos_file no -cov_cutoff auto ";
    return $velvetg_cmd;
}

sub get_velvet_cmds
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $kmer = shift;
    my $total_coverage = shift;
    my $avg_readlen = shift;
    my $assembly_outdir = shift;
    
    my $trdata = $trimraw . "data";
    my $kmer_bin = Assembly::Velvet::get_kmer_bin($kmer);
    my $exp_cov = Assembly::Velvet::calc_exp_cov($total_coverage, $avg_readlen, $kmer);
    
    # Create the kmer directory
    my $kmer_dir = $assembly_outdir . "/" . $trimraw . "/assem_kmer-" . $kmer . "_exp_cov-" . $exp_cov . "_covcutoff-auto";
    unless (-e $kmer_dir) {
        mkpath $kmer_dir;
    }
    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "kmer_dir", $kmer_dir);
    
    my $velveth_cmd = base_velveth_cmd($kmer, $kmer_bin, $kmer_dir, $kmer);
    my $velvetg_cmd = base_velvetg_cmd($kmer_bin, $kmer_dir, $exp_cov);
    my $i = 0;
    for my $sample_type (qw(PE PER MP MP3 MP8)) {
        my $sample_ref = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type]);
        if ($sample_ref) {
            $velveth_cmd .= add_velveth_idx($sample_ref, $sample_type, $trimraw, $i);
            $velvetg_cmd .= add_velvetg_idx($sample_ref, $sample_type, $trimraw, $i);
            $i++;
        }
    }
    if ($i) {
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "velveth_cmd", $velveth_cmd);
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer], "velvetg_cmd", $velvetg_cmd);
    } else {
        $velveth_cmd = '';
        $velvetg_cmd = '';
    }
    return ($velveth_cmd, $velvetg_cmd);
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

sub submit_cmds
{
    my $vqs = shift;
    my $rec = shift;
    my $trimraw = shift;
    my $kdir = Assembly::Utils::get_check_record($rec, ["kmer_dir"]);
    my $vh_cmd = Assembly::Utils::get_check_record($rec, ["velveth_cmd"]);
    my $vg_cmd = Assembly::Utils::get_check_record($rec, ["velvetg_cmd"]);
    
    my $vh_files_exist = outfiles_exist($kdir, $vh_outfiles);
    my $vg_files_exist = outfiles_exist($kdir, $vg_outfiles);

    if ($vh_files_exist and !$vg_files_exist) {
        my $sub_cmds = $vqs->submit_vg($vh_cmd, $trimraw);
        my ($vg_qsub_cmd, $vg_jobid) = $sub_cmds;
        Assembly::Utils::set_check_record($rec, [], "velvetg_qsub_cmd", $vg_qsub_cmd);
        Assembly::Utils::set_check_record($rec, [], "velvetg_qsub_jobid", $vg_jobid);
    } elsif (!$vh_files_exist) {
        my $sub_cmds = $vqs->submit_vhg($vh_cmd, $vg_cmd, $trimraw);
        my ($vh_qsub_cmd, $vh_jobid, $vg_qsub_cmd, $vg_jobid) = @$sub_cmds;
        Assembly::Utils::set_check_record($rec, [], "velveth_qsub_cmd", $vh_qsub_cmd);
        Assembly::Utils::set_check_record($rec, [], "velveth_qsub_jobid", $vh_jobid);
        Assembly::Utils::set_check_record($rec, [], "velvetg_qsub_cmd", $vg_qsub_cmd);
        Assembly::Utils::set_check_record($rec, [], "velvetg_qsub_jobid", $vg_jobid);
    } elsif ($vh_files_exist and $vg_files_exist) {
        print_verbose "Not running for this kmer - found valid velvetg output files\n";
    }
}
        

sub build_assembly_cmds
{
    my $records = shift;
    my $vqs = new Assembly::Qsub($options->{qsub_script}, $options->{submit}, $options->{submit_max}, $options->{verbose});
    for my $species (keys %$records) {
        my $spec_ref = $records->{$species}->{DNA};
        for my $strain (keys %$spec_ref) {
			for my $trimraw (qw(trim raw)) {
				if ($options->{$trimraw}) {
                    my $assembly_dir = Assembly::Velvet::get_assembly_outdir($records, $species, $strain, $trimraw);
                    if ($assembly_dir) {
                        my ($total_coverage, $avg_readlen) = Assembly::Velvet::get_coverage_vars($records, $species, $strain, $trimraw);
                        my $kmer_list = Assembly::Velvet::create_kmer_range($records, $species, $strain, $trimraw);
                        for my $kmer (@$kmer_list) {
                            my ($velveth_cmd, $velvetg_cmd) = get_velvet_cmds($records, $species, $strain, $trimraw, $kmer, $total_coverage, $avg_readlen, $assembly_dir);
                            my $rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer]);
                            if ($rec) {
                                if ($velveth_cmd !~ /shortPaired2/ and $velvetg_cmd !~ /ins_length2/) {
                                    submit_cmds($vqs, $rec, $trimraw);
                                }
                            }
                        }
					}
				}
			}
		}
	}
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
build_assembly_cmds($records);;
DumpFile($options->{yaml_out}, $records);
