#!/usr/bin/env perl
use strict;
use warnings;
use YAML::XS qw(LoadFile DumpFile);
use File::Path;
use Getopt::Long;

# From the records, pull just the genomic data
# Read in genome lengths from file
# Write out table of velvet assembly params
# Write out qsubs using these commands.

# Ideally: have separate folders within velvet that include
# e.g. the different genome lengths, QC trims etc. that we try?

$options = {};
$records = {};
my @genomic_samples = ();
my @glen_col_headers = qw(IL_Species IL_Genotype IL_Biomaterial IL_Biomaterial_type 
    IL_Sample_ID IL_Kingdom IL_Resident_Expert RG_Species RG_Strain RG_Est_Genome_Length
     RG_Source RG_Citation RG_URL RG_Notes);

my $tr;
my $trdata;
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
my $velvet_bin_dir = "/opt/bio/velvet";
my @kbins = (0, 31, 63, 127, 145); 
# @ kbins is only used by get_kmer_bin function below.
# Assumes we have binaries of form
# velvetg_x and velveth_x for x>0 and x in @kbins.
                                
my @vg_qsub_list = ();
my @vg_qsub_list = ();

sub set_default_opts
{
    my %defaults = qw(
            min_kmer 21 
            max_kmer 95
            vh_batch_dir qsub_files/04_vh_cmds
            vg_batch_dir qsub_files/04_vg_cmds
            qsub_script qsub_array.sh
            vh_qsub_batch qsub_files/04_vh_qsubs.sh
            vg_qsub_batch qsub_files/04_vg_qsubs.sh
            trim 1
            genome_length_file input_data/GenomeLengthEst.tab
            );
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
    if ($options->{raw}) {
        $tr = "raw";
        $trdata = "rawdata";
    } else {
        $tr = "trim";
        $trdata = "trimdata";
    }
    $options->{qsub_opts} = '' unless ($options->{qsub_opts});
}

sub check_opts
{
    unless ($opts->{yaml_in} and $opts->{yaml_out}) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file>
            Optional:
                --genome_length_file <filename>
                --trim
                --raw
                --qsub_opts <qsub options string in quotes>
                --qsub_script <script name>
                --vh_qsub_batch <output batch filename>
                --vg_qsub_batch <output batch filename>
                --submit
                --verbose
                --min_kmer <value (default 21)>
                --max_kmer <value (default 95)>
                --vh_batch_file <filename>
                --vg_batch_file <filename>
                ";
    }
}


sub gather_opts
{
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'genome_length_file|g=s',
            'trim',
            'raw',
            'qsub_opts=s',
            'qsub_script=s',
            'qsub_batch=s',
            'submit',
            'verbose',
            'min_kmer',
            'max_kmer',
            'vh_batch_file=s',
            'vg_batch_file=s',
            );
    check_opts;
    set_default_opts;
}

sub parse_genome_lengths
{
    my $fname = ($options->{genome_length_file} ? $options->{genome_length_file} : '');
    if ($fname) {
        open (FGLEN, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FGLEN>; <FGLEN>; # Skip first two col header lines.
        while (my $line = <FGLEN>) {
            chomp $line;
            my @fields = split(/\t/, $line);
            my %fh = map { $glen_col_headers[$i] => ($fields[$i] ? $fields[$i] : '') } (0..$#glen_col_headers);
            my $sample = $fh->{IL_Sample_ID};
            my $rec = ($records->{$sample} ? $records->{$sample} : '');
            if ($rec) {
                $rec->{related_genome_length} = {};
                my $rgl_ref = $rec->{related_genome_length};
                for my $ch (@glen_col_headers) {
                    if ($ch =~ /^RG/) {
                        $rgl_ref->{$ch} = $fh->{$ch};
                    }
                }
            } else {
                print "Warning: in file $fname, sample $sample was not found in input yaml file.\n";
            }
        }
        close (FGLEN);
    }
}

sub get_genomic_records
{
    for my $sample (keys %$records) {
        if ($records->{$species}->{bio_type} =~ /DNA/) {
            push (@genomic_samples, $sample);
        }
    }
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

sub get_coverage_vars
{
    my $rec = shift;
    my $var = {};
    my $var->{R1_nreads} = get_check_record($rec, ["data_stats", "R1", $trdata, "num_reads"]);
    my $var->{R1_readlen} = get_check_record($rec, ["data_stats", "R1", $trdata, "read_length"]);
    my $var->{R2_nreads} = get_check_record($rec, ["data_stats", "R2", $trdata, "num_reads"]);
    my $var->{R2_readlen} = get_check_record($rec, ["data_stats", "R2", $trdata, "read_length"]);
    my $var->{genome_length} = get_check_record($rec, ["related_genome_length", "RG_Est_Genome_Length"]);
    my $var->{total_coverage} = ($var->{R1_nreads} * $var->{R1_readlen} + $var->{R2_nreads} * $var->{R2_readlen}) / $var->{genome_length};
    my $var->{avg_readlen} = ($var->{R1_readlen} + $var->{R2_readlen}) / 2;
}

sub calc_exp_cov
{
    my $var = shift;
    my $kmer = shift;
    my $exp_cov_float = $var->{total_coverage} * ($var->{avg_readlen} - $kmer + 1) / $var->{avg_readlen};
    my $exp_cov_rounded_int = int($exp_cov_float + 0.5);
    return $exp_cov_rounded_int;
} 

sub get_assembly_outdir
{
    unless (defined ($rec->{velvet})) {
        $rec->{velvet} = {};
    }
    unless (defined ($rec->{velvet}->{$tr})) {
        $rec->{velvet}->{$tr};
    }
    my $assembly_outdir = $rec->{sample_dir} . "/assemblies/velvet/$tr";
    mkpath $assembly_outdir;
    $rec->{velvet}->{$tr}->{assembly_outdir} = $assembly_outdir;
    return $assembly_outdir;
}

sub get_velvet_kdir
{
    my $assembly_outdir = shift;
    my $kmer = shift;
    my $exp_cov = shift;
    my $kdir = $assembly_outdir . "/assem_kmer-" . $kmer . "_exp-" . $exp_cov . "_covcutoff-auto";
    mkpath $kdir;
    return $kdir;
}

sub get_kmer_bin
{
    my $kmer = shift;
    my $bin = '';
    if ($kmer < $kbins[0]) { 
        print "Error: kmer value $kmer must be a positive integer!\n";
    } elsif ($kmer > $kbins[$#kbins]) {
        print "kmer value $kmer not supported. Recompile velvet for higher kmer.\n";
    } else {
        for (my $i = 0; $i < $#kbins; $i++) {
            if ($kbins[$i] < $kmer and $kmer <= $kbins[$i+1]) {
                $bin = $kbins[$i+1];
            }
        }
    }
    return $bin;
}
    

sub get_velveth_cmd
{
    my $rec = shift;
    my $kmer = shift;
    my $kmer_bin = shift;
    my $outdir = shift;
    my $velveth_bin = $velvet_bin_dir . "/velveth_" . $kmer_bin;
    my $r1_file = $rec->{R1}->{$trdata};
    my $r2_file = $rec->{R2}->{$trdata};
    my $velveth_cmd = $velveth_bin . " " . $outdir . " " . $kmer . 
        "  -fastq -shortPaired -create_binary -separate " . $r1_file . " " . $r2_file;
    return $velveth_cmd;
}

sub get_velvetg_cmd
{
    my $exp_cov = shift;
    my $kmer = shift;
    my $kmer_bin = shift;
    my $working_dir = shift;
    my $min_contig_opt = " ";
    my $scaffolding_opt = " -scaffolding yes ";
    my $velvetg_bin = $velvet_bin_dir . "/velvetg_" . $kmer_bin;
    my $velveth_cmd = $velvetg_bin . " " . $working_dir . " " .
        "-ins_length 300 -exp_cov " . $exp_cov . " " $min_contig_opt . 
        $scaffolding_opt . " -amos_file no -cov_cutoff auto";
    return $velvetg_cmd;
}

sub get_velvet_cmds
{
    my $rec = shift;
    my $cov_vars = shift;
    my $assembly_outdir = shift;
    my $kmer = shift;
    my $kmer_bin = get_kmer_bin($kmer);
    my $exp_cov = calc_exp_cov($cov_vars, $kmer);
    my $kdir = get_velvet_kdir($assembly_outdir, $kmer, $exp_cov);
    my $vh_cmd = get_velveth_cmd($rec, $cov_vars->{avg_readlen}, $kmer, $kmer_bin, $kdir);
    my $vg_cmd = get_velvetg_cmd($exp_cov, $kmer, $kmer_bin, $kdir);
    return ($vh_cmd, $vg_cmd);
}

sub vh_write_batch
{
    my $vh_batch_dir = $options->{vh_batch_dir};
    open (FHBATCH, '>', $vh_batch_file) or die "Error: couldn't open file " . $vh_batch_file . "\n";
    print FHBATCH join ("\n", @vh_list) . "\n";
    close (FHBATCH);
}

sub vg_write_batch
{
    my $sample = shift;
    my $vg_batch_dir = $options->{vg_batch_dir};
    my $sample_batch
    open (FGBATCH, '>', $vg_batch_file) or die "Error: couldn't open file " . $vg_batch_file . "\n";
    print FGBATCH join ("\n", @vg_list) . "\n";
    close (FGBATCH);
}

sub get_vh_qsub
{
    my $rec = shift;
    my $aref = shift;
    my $sample = $rec->{sample};
    my @vh_cmd_list = @$aref;
    my $num_cmds = scalar @vh_cmd_list;
    my $batch_filename = vh_write_batch($aref);
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N " . $tr . "_velveth -t 1:" . $num_cmds . " " . $options->{qsub_script} . " " . $batch_filename;
    push (@vh_qsub_list, $qsub_cmd);
}

sub get_vg_qsub
{
    my $rec = shift;
    my $aref = shift;
    my $sample = $rec->{sample};
    my @vh_cmd_list = @$aref;
    my $num_cmds = scalar @vg_cmd_list;
    my $batch_filename = vg_write_batch($aref);
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N " . $tr . "_velvetg -t 1:" . $num_cmds . " " . $options->{qsub_script} . " " . $batch_filename;
    push (@vg_qsub_list, $qsub_cmd);
}

sub build_assembly_cmds
{
    for my $sample (@genomic_samples) {
        my $rec = ($records->{$sample} ? $records->{$sample} : '');
        unless ($rec) die "Couldn't get a record from yaml file for sample $sample.\n";
        my $cov_vars = get_coverage_vars($rec);
        my $assembly_outdir = get_assembly_outdir;
        my @vh_cmd_list = ();
        my @vg_cmd_list = ();
        for my $kmer ($options->{min_kmer}..$options->{max_kmer}) {
            my ($vh_cmd, $vg_cmd) = get_velvet_cmds($rec, $cov_vars, $assembly_outdir, $kmer);
            push (@vh_cmd_list, $vh_cmd);
            push (@vg_cmd_list, $vg_cmd);
        }
        my $vh_qsub_cmd = get_vh_qsub($rec, \@vh_cmd_list);
        my $vg_qsub_cmd = get_vg_qsub($rec, \@vg_cmd_list);
        push (@vh_qsub_list, $vh_qsub_cmd);
        push (@vg_qsub_list, $vg_qsub_cmd);
    }
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

sub create_qsubs
{
    # Go through vh's 
    # submit all for a sample, then record jid
    # then add hold_jid to the vg's qsubs
}


$records = LoadFile($options->{yaml_in});
get_genomic_records;
parse_genome_lengths;
build_assembly_cmds;
vh_write_qsub_batch;
vg_write_qsub_batch;
DumpFile($options->{yaml_out}, $records);
