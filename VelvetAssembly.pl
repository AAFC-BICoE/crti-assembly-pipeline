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
my @genomic_keys = ();
my @glen_col_headers = qw(IL_Species IL_Genotype IL_Biomaterial IL_Biomaterial_type 
    IL_Sample_ID IL_Kingdom IL_Resident_Expert RG_Species RG_Strain RG_Est_Genome_Length
     RG_Source RG_Citation RG_URL RG_Notes);


sub check_opts
{
    unless ($opts->{yaml_in} and $opts->{yaml_out}) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file>
            Optional:
                --genome_length_file <filename>
                --no_qsub
                --qsub_opts <qsub options string in quotes>
                --qsub_script <script name>
                --qsub_batch <output batch filename>
                --verbose
                ";
    }
}


sub gather_opts
{
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'genome_length_file|g=s',
            'verbose',
            'qsub_opts=s',
            'qsub_script=s',
            'qsub_batch=s',
            );
    check_opts;
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
            push (@genomic_keys, $sample);
        }
    }
}

sub getExp { 
    # nucCoverage, seq length, kmer size 
    my $k = shift;
    my $C = $options->{num_reads} * $options->{read_length} / $options->{genome_length} 
    my $L = $options->{read_length};

    my $exp_raw = $C * ($L - $k + 1)/$L;
    my $exp_rnd = int($exp_raw + 0.5); # round expected coverage value.
    return $exp_rnd;
}

sub build_assembly_cmds
{
    for my $sample (@genomic_keys) {
        my $rec = $records->{$sample};
        my $r1_numreads = $rec->{data_stats}->{R1}->{trimdata}->{num_reads};
        my $r1_read_length = $rec->{data_stats}->{R1}->{trimdata}->{read_length};
        my $r2_numreads = $rec->{data_stats}->{R2}->{trimdata}->{num_reads};
        my $r2_read_length = $rec->{data_stats}->{R2}->{trimdata}->{read_length};
        my $genome_length = $rec->{related_genome_length}->{RG_Est_Genome_Length};
        my $assembly_outdir = $rec->{sample_dir} . "/assemblies/velvet";
        mkpath $assembly_outdir;


$records = LoadFile($options->{yaml_in});
get_genomic_records;
parse_genome_lengths;
build_assembly_cmds;
DumpFile($options->{yaml_out}, $records);
