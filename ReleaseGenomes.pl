#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use POSIX qw(strftime);
use Assembly::Utils;
use Cwd;
use File::Path;

my $options = {};
my @release_headers = qw (Species Strain Sequencing_Types Trim_Raw Release_Version Est_Genome_Size Best_Kmer Best_N50);
my @wiki_genome = ();
my @wiki_release = ();
my @wiki_combined = ();

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/13_contig_stats.yml
        yaml_out yaml_files/14_velvet_release.yml
        release_table input_data/Release.tab
        verbose 1
        testing 0
        copy_reads 1
        link_reads 0
        no_version 1
        wiki_release_table output_files/wiki_release_table.txt
        wiki_genome_table output_files/wiki_genome_table.txt
        wiki_combined_table output_files/wiki_combined_table.txt
    );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file> -r <release table>
            Optional:
                --verbose
                --testing (place output yaml files in /tmp)
                --wiki_release_table (print out table of releases for pasting to wiki)
                --wiki_genome_table (print table of genome stats for pasting to wiki)
                --wiki_combined_table (print combined release/genome stats for pasting to wiki)
            ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'release_table|r=s',
        'verbose',
        'testing',
        'copy_reads',
        'link_reads',
        'no_version',
        'wiki_release_table',
        'wiki_genome_table',
        'wiki_combined_table',
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

# Function to add commas to a number, from: 
# http://perldoc.perl.org/perlfaq5.html#How-can-I-output-my-numbers-with-commas-added%3F
sub commify {
    local $_ = shift;
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
}

sub get_software_versions
{
    my $versions = {};
    if ($options->{no_version}) {
        $versions->{fastqc} = "FastQC v0.10.1\n";
        $versions->{fastx_trimmer} = "Part of FASTX Toolkit 0.0.13.2 by A. Gordon (gordon\@cshl.edu)\n";
        $versions->{velvetk} = "2012\n";
        $versions->{velveth} = "Version 1.2.08\n";
        $versions->{velvetg} = "Version 1.2.08\n";
    } else {
        $versions->{fastqc} = `fastqc --version`;
        $versions->{fastx_trimmer} = `fastx_trimmer -h | grep FASTX`;
        $versions->{velvetk} = "2012\n";
        $versions->{velveth} = `velveth_31 | grep Version`;
        $versions->{velvetg} = `velvetg_31 | grep Version`;
    }
    return $versions;
}

sub parse_release_table
{
    my $fname = shift;
    my $base_dir = shift;
    my $release_candidates = {};
    unless ($fname) { return; }
    open (RTAB, '<', $fname) or die "Error: couldn't open release table file $fname\n";
    <RTAB>; # skip header
    while (my $line = <RTAB>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        my $num_headers = scalar @release_headers;
        my $num_fields = scalar @fields;
        my $rec = {};
        for (my $i=0; $i<$num_headers; $i++) { 
            my $key = $release_headers[$i];
            my $val = ($fields[$i] ? $fields[$i] : '');
            $rec->{$key} = $val; 
        }
        my $species = Assembly::Utils::format_species_key($rec->{Species});
        my $strain = Assembly::Utils::format_strain_key($rec->{Strain});
        $release_candidates->{$species} = {} unless $release_candidates->{$species};
        $release_candidates->{$species}->{$strain} = $rec;
        #print "Added release candidate record for species $species strain $strain\n";
    }
    close (RTAB);
    return $release_candidates;
}

sub get_release_prefix
{
    my $species = shift;
    my $strain = shift;
    my $cand_rec = shift;
    my $strain_rec = shift;
    
    my $species_abbr = Assembly::Utils::get_check_record($strain_rec, ["PE", "species_abbr"]);
    my $cand_release_version = $cand_rec->{Release_Version};
    my $release_prefix = $species_abbr . "_" . $strain . "_" . $cand_release_version;
    return $release_prefix;
}

sub release_exists
{
    my $release_dir = shift;
    my $release_prefix = shift;
    
    unless (-e $release_dir) {
        return 0;
    }
    my $release_yaml_file = $release_dir . "/" . $release_prefix . ".yml";
    unless (-e $release_yaml_file and -s $release_yaml_file) {
        # yml file doesn's exist or has size=0
        return 0;
    }
    # check for any other files here ... reads / contigs?
}   

sub get_release_stanza
{
    my $species = shift;
    my $strain = shift;
    my $inrec = shift;
    my $outrec = {};
    $outrec->{species} = Assembly::Utils::get_check_record($inrec, ["PE", "sequencing_metadata", "Organism"]);
    $outrec->{strain} = Assembly::Utils::get_check_record($inrec, ["PE", "sequencing_metadata", "Genotype"]);
    $outrec->{version} = Assembly::Utils::get_check_record($inrec, ["release", "version"]);
    $outrec->{date} = strftime "%m/%d/%Y", localtime;
    return $outrec;
}
    
sub get_genome_stanza
{
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $inrec = shift;
    my $release_kmer = shift;
    my $outrec = {};
    $outrec->{kingdom} = Assembly::Utils::get_check_record($inrec, ["related_genome_length", "IL_Kingdom"]);
    $outrec->{species_dir} = Assembly::Utils::get_check_record($inrec, ["PE", "species_dir"]);
    $outrec->{estimated_genome_length} = Assembly::Utils::get_check_record($inrec, ["related_genome_length", "RG_Est_Genome_Length"]);
    $outrec->{release_kmer} = $release_kmer;
    $outrec->{N50} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "N50"]);
    $outrec->{max_contig} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "max_contig"]);
    $outrec->{reads_used} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "reads_used"]);
    $outrec->{total_length} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "total_length"]);
    $outrec->{num_contigs} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "num_contigs"]);
    $outrec->{min_contig_len} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "min_contig_len"]);
    $outrec->{median_contig_len} = Assembly::Utils::get_check_record($inrec, ["velvet", $trimraw, "kmer", $release_kmer, "median_contig"]);
    my @types = ();
    for my $st (qw(PE PER MP MP3 MP8)) {
    #for my $st (qw(PE)) {
        if (Assembly::Utils::get_check_record($inrec, [$st])) {
            push (@types, $st);
        }
    }
    $outrec->{assembly_sample_types} = join (",", @types);
    return $outrec;
}

# Determine whether to copy or link reads files to release dir
# based on input options.
sub copy_link_reads
{
    my $file_in = shift;
    my $file_out = shift;
    unless (-e $file_out and (-s $file_in) == (-s $file_out)) {
        print_verbose "cp/ln -s $file_in $file_out\n";
        unless ($options->{testing}) {
            if ($options->{link_reads}) {
                symlink($file_in, $file_out);
            } elsif ($options->{copy_reads}) {
                system("cp $file_in $file_out");
            }
        }
    }
} 

sub get_sample_stanza
{
    my $strain_rec = shift;
    my $sample_type = shift;
    my $trimraw = shift;
    my $sample_rec = Assembly::Utils::get_check_record($strain_rec, [$sample_type]);
    my $versions = get_software_versions;
    my $outrec = {};
    $outrec->{sequencing_metadata} = Assembly::Utils::get_check_record($sample_rec, ["sequencing_metadata"]);
    # $outrec->{species_abbr} = Assembly::Utils::get_check_record($sample_rec, ["species_abbr"]);
    $outrec->{read_data} = {};
    $outrec->{read_data}->{sample_dir} = Assembly::Utils::get_check_record($sample_rec, ["sample_dir"]);
    $outrec->{read_data}->{R1} = {};
    $outrec->{read_data}->{R1}->{raw} = {};
    $outrec->{read_data}->{R1}->{raw}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R1", "rawdata"]);
    $outrec->{read_data}->{R1}->{raw}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "rawdata", "num_reads"]);
    $outrec->{read_data}->{R1}->{raw}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "rawdata", "read_length"]);
    $outrec->{read_data}->{R2} = {};
    $outrec->{read_data}->{R2}->{raw} = {};
    $outrec->{read_data}->{R2}->{raw}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R2", "rawdata"]);
    $outrec->{read_data}->{R2}->{raw}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "rawdata", "num_reads"]);
    $outrec->{read_data}->{R2}->{raw}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "rawdata", "read_length"]);

    # get info for output file path for release raw reads
    my $bio_type = Assembly::Utils::get_check_record($sample_rec, ["bio_type"]);
    my $sample_id = Assembly::Utils::get_check_record($sample_rec, ["sample"]);
    my $release_prefix = Assembly::Utils::get_check_record($strain_rec, ["release", "prefix"]);
    my $release_dir = Assembly::Utils::get_check_record($strain_rec, ["release", "release_dir"]);
    my $fpath_base = $release_dir . "/" . $release_prefix . "_" . $bio_type . "_" . $sample_type . 
            "_" . $sample_id . "_";

    my $r1_in = Assembly::Utils::get_check_record($sample_rec, ["R1", "rawdata"]);
    my $r1_out = $fpath_base . "R1.fq";
    if ($r1_in =~/\.gz\s*$/){ $r1_out .= ".gz"; }
    my $r2_in = Assembly::Utils::get_check_record($sample_rec, ["R2", "rawdata"]);
    my $r2_out = $fpath_base . "R2.fq";
    if ($r2_in =~/\.gz\s*$/){ $r2_out .= ".gz"; }

    $outrec->{release} = [];
    $outrec->{release}->[0] = {};
    $outrec->{release}->[0]->{input_file} = $r1_in;
    $outrec->{release}->[0]->{output_file} = $r1_out;
    $outrec->{release}->[1] = {};
    $outrec->{release}->[1]->{input_file} = $r2_in;
    $outrec->{release}->[1]->{output_file} = $r2_out;

    # copy the files here??
    copy_link_reads($r1_in, $r1_out);
    copy_link_reads($r2_in, $r2_out);

    if ($trimraw =~ /trim/) {
        $outrec->{read_data}->{R1}->{trim} = {};
        $outrec->{read_data}->{R1}->{trim}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R1", "trimdata"]);
        $outrec->{read_data}->{R1}->{trim}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "trimdata", "num_reads"]);
        $outrec->{read_data}->{R1}->{trim}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R1", "trimdata", "read_length"]);
        $outrec->{read_data}->{R2}->{trim} = {};
        $outrec->{read_data}->{R2}->{trim}->{file_path} = Assembly::Utils::get_check_record($sample_rec, ["R2", "trimdata"]);
        $outrec->{read_data}->{R2}->{trim}->{num_reads} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "trimdata", "num_reads"]);
        $outrec->{read_data}->{R2}->{trim}->{read_length} = Assembly::Utils::get_check_record($sample_rec, ["data_stats", "R2", "trimdata", "read_length"]);

        my $tr1_in = Assembly::Utils::get_check_record($sample_rec, ["R1", "trimdata"]);
        my $tr1_out = $fpath_base . "trim_R1.fq";
        if ($tr1_in =~/\.gz\s*$/){ $tr1_out .= ".gz"; }
        my $tr2_in = Assembly::Utils::get_check_record($sample_rec, ["R2", "trimdata"]);
        my $tr2_out = $fpath_base . "trim_R2.fq";
        if ($tr2_in =~/\.gz\s*$/){ $tr2_out .= ".gz"; }
        
        $outrec->{release}->[2] = {};
        $outrec->{release}->[2]->{input_file} = $tr1_in;
        $outrec->{release}->[2]->{output_file} = $tr1_out;
        $outrec->{release}->[3] = {};
        $outrec->{release}->[3]->{input_file} = $tr2_in;
        $outrec->{release}->[3]->{output_file} = $tr2_out;
        
        copy_link_reads ($tr1_in, $tr1_out);
        copy_link_reads ($tr2_in, $tr2_out);
    }
    
    # Now add the fastqc info to the pipeline
    my $qcr1 = {};
    $qcr1->{description} = "FastQC";
    $qcr1->{version} = $versions->{fastqc};
    chomp $qcr1->{version};
    $qcr1->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "raw_cmd"]);
    $qcr1->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "raw_qsub"]);
    $qcr1->{run_dir} = getcwd;
    
    my $qcr2 = {};
    $qcr2->{description} = "FastQC";
    $qcr2->{version} = $versions->{fastqc};
    chomp $qcr2->{version};
    $qcr2->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "raw_cmd"]);
    $qcr2->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "raw_qsub"]);
    $qcr2->{run_dir} = getcwd;
    
    my $qtr1 = {};
    $qtr1->{description} = "Fastx_trimmer";
    $qtr1->{version} = $versions->{fastx_trimmer};
    chomp $qtr1->{version};
    $qtr1->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R1", "trim_cmd"]);
    $qtr1->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R1", "qsub_cmd"]);
    $qtr1->{run_dir} = getcwd;
    
    my $qtr2 = {};
    $qtr2->{description} = "Fastx_trimmer";
    $qtr2->{version} = $versions->{fastx_trimmer};
    chomp $qtr2->{version};
    $qtr2->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R2", "trim_cmd"]);
    $qtr2->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastx_trimmer", "R2", "qsub_cmd"]);
    $qtr2->{run_dir} = getcwd;    
    
    # Now add the fastqc info to the pipeline
    my $qct1 = {};
    $qct1->{description} = "FastQC";
    $qct1->{version} = $versions->{fastqc};
    chomp $qct1->{version};
    $qct1->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "trim_cmd"]);
    $qct1->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R1", "trim_qsub"]);
    $qct1->{run_dir} = getcwd;
    #$qct1->{release} = [];
    #$qct1->{release}->[0] = {};
    #$qct1->{release}->[0]->{input_file} = Assembly::Utils::get_check_record($sample_rec, ["R1", "trimdata"]);
    #$qct1->{release}->[0]->{output_file} = 
    # this is handled in the location above... should it be??
    
    my $qct2 = {};
    $qct2->{description} = "FastQC";
    $qct2->{version} = $versions->{fastqc};
    chomp $qct2->{version};
    $qct2->{command} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "trim_cmd"]);
    $qct2->{qsub_cmd} = Assembly::Utils::get_check_record($sample_rec, ["fastqc", "R2", "trim_qsub"]);
    $qct2->{run_dir} = getcwd;
    
    my $qc =[];
    push (@$qc, $qcr1);
    push (@$qc, $qcr2);
    if ($trimraw =~ /trim/) {
        push (@$qc, $qtr1);
        push (@$qc, $qtr2);
        push (@$qc, $qct1);
        push (@$qc, $qct2);
    }
    return ($outrec, $qc);
}

sub get_pipeline_stanza
{
    my $strain_rec = shift;
    my $trimraw = shift;
    my $release_kmer = shift;
    my $release_dir = shift;
    my $release_prefix = shift;
    my $outrec = []; # note: it's an array this time
    my $krec = Assembly::Utils::get_check_record($strain_rec, ["velvet", $trimraw, "kmer", $release_kmer]);
    my $versions = get_software_versions;
    
    my $vh = {};
    $vh->{description} = "velveth";
    $vh->{run_dir} = getcwd;
    $vh->{version} = $versions->{velveth};
    chomp $vh->{version};
    $vh->{command} = Assembly::Utils::get_check_record($krec, ["velveth_cmd"]);
    $vh->{qsub_cmd} = Assembly::Utils::get_check_record($krec, ["velveth_qsub_cmd"]);
    
    my $vg = {};
    $vg->{description} = "velvetg";
    $vg->{run_dir} = getcwd;
    $vg->{version} = $versions->{velvetg};
    chomp $vg->{version};
    $vg->{command} = Assembly::Utils::get_check_record($krec, ["velvetg_cmd"]);
    $vg->{qsub_cmd} = Assembly::Utils::get_check_record($krec, ["velvetg_qsub_cmd"]);
    my $kmer_dir = Assembly::Utils::get_check_record($krec, ["kmer_dir"]);
    my $contigs_infile = $kmer_dir . "/contigs.fa";
    my $contigs_outfile = $release_dir . "/" . $release_prefix . "_contigs.fa";
    $vg->{release} = [];
    $vg->{release}->[0] = {};
    $vg->{release}->[0]->{input_file} = $contigs_infile;
    $vg->{release}->[0]->{output_file} = $contigs_outfile;
    
    unless (-e $contigs_outfile and -s $contigs_outfile) {
        print_verbose "cp $contigs_infile $contigs_outfile\n";
        unless ($options->{testing}) {
            system ("cp $contigs_infile $contigs_outfile");
        }
    }
    
    push (@$outrec, $vh);
    push (@$outrec, $vg);
    return $outrec;
}   

sub create_release_yaml
{
    my $species = shift;
    my $strain = shift;
    my $cand_rec = shift;
    my $strain_rec = shift;
    #my $release_yaml_rec = {};
    my $release_yaml_rec = [];
    my $trimraw = $cand_rec->{Trim_Raw};
    my $release_kmer = $cand_rec->{Best_Kmer};
    
    $release_yaml_rec->[0]->{release} = get_release_stanza($species, $strain, $strain_rec);
    $release_yaml_rec->[1]->{genome_assembly} = get_genome_stanza($species, $strain, $trimraw, $strain_rec, $release_kmer);
    my $pipeline = [];
    $release_yaml_rec->[2]->{samples} = [];
    for my $sample_type (qw(PE PER MP MP3 MP8)) {
        my $sample_rec = Assembly::Utils::get_check_record($strain_rec, [$sample_type]);
        if ($sample_rec) {
            my ($sample_stanza, $qc_cmds) = get_sample_stanza($strain_rec, $sample_type, $trimraw);
            #$release_yaml_rec->{$sample_type} = $sample_stanza;
            push (@{$release_yaml_rec->[2]->{samples}}, $sample_stanza);
            for my $rec (@$qc_cmds) {
                push (@$pipeline, $rec);
            }
        }
    }
    # Now get the rest of the pipeline
    my $release_prefix = Assembly::Utils::get_check_record($strain_rec, ["release", "prefix"]);
    my $release_dir = Assembly::Utils::get_check_record($strain_rec, ["release", "release_dir"]);
    my $pipeline_remainder = get_pipeline_stanza($strain_rec, $trimraw, $release_kmer, $release_dir, $release_prefix);
    for my $rec (@$pipeline_remainder) {
        push (@$pipeline, $rec);
    }
    $release_yaml_rec->[3]->{pipeline} = $pipeline;
    return $release_yaml_rec;
}

sub get_wiki_release_line
{
    my $release = shift;
    my $release_dir = shift;
    my $species = $release->[0]->{release}->{species};
    my $strain = $release->[0]->{release}->{strain};
    my $version = $release->[0]->{release}->{version};
    my $type = "Genome";
    my $link_target = "http://biocluster/project_data/CRTI-09S-462RD/specimen/";
    if ($release_dir =~ /specimen\/(.*)/) {
        $link_target = $link_target . $1;
    } else {
        $link_target = $link_target . $release_dir;
    }
    my $link_text = "Release"; # could also be the release prefix here.
    my $link = "[[" . $link_target . "][" . $link_text . "]]";
    my $line = "|" . join("|", ($species, $strain, $version, $type, $link)) . "|";
    push (@wiki_release, $line);
}

sub get_wiki_genome_line
{
    my $release = shift;
    my $species = $release->[0]->{release}->{species};
    my $strain = $release->[0]->{release}->{strain};
    my $N50 = $release->[1]->{genome_assembly}->{N50};
    my $max_contig = $release->[1]->{genome_assembly}->{max_contig};
    my $reads_used = $release->[1]->{genome_assembly}->{reads_used};
    my $total_length = $release->[1]->{genome_assembly}->{total_length};
    my $est_genome_length = $release->[1]->{genome_assembly}->{estimated_genome_length};
    my @fields = ($species, $strain, $est_genome_length, $N50, $max_contig, $reads_used,
            $total_length);
    my $line = "|" . join("|", @fields) . "|";
    push (@wiki_genome, $line);
}

sub get_wiki_combined_line
{
    my $release = shift;
    my $release_dir = shift;
    my @fields = ();
    push (@fields, $release->[1]->{genome_assembly}->{kingdom});
    push (@fields, $release->[0]->{release}->{species});
    push (@fields, $release->[0]->{release}->{strain});
    push (@fields, $release->[0]->{release}->{version});
    push (@fields, "Genome"); #Type

    push (@fields, "yes"); # paired-end = true
    push (@fields, "no"); # MP 3kb = false
    push (@fields, "no"); # MP 8kb = false
    push (@fields, $release->[1]->{genome_assembly}->{release_kmer}); # Kmer used for release

    my $link_target = "http://biocluster/project_data/CRTI-09S-462RD/specimen/";
    if ($release_dir =~ /specimen\/(.*)/) {
        $link_target = $link_target . $1;
    } elsif ($release_dir =~ /\.\.\/\.\.\/processing_test2\/(.*)/) {
        $link_target = $link_target . $1;
    } else {
        $link_target = $link_target . $release_dir;
    }
    my $link_text = "Release"; # could also be the release prefix here.
    my $link = "[[" . $link_target . "][" . $link_text . "]]";
    push (@fields, $link);
    
    for my $key (qw(total_length estimated_genome_length num_contigs min_contig_len median_contig_len
                max_contig N50)) {
        my $var = ($release->[1]->{genome_assembly}->{$key} ? $release->[1]->{genome_assembly}->{$key} : '');
        push (@fields, commify($var));
    }
    
    # Add commas to both numerator and denominator of reads used.
    my $reads_used = ($release->[1]->{genome_assembly}->{reads_used} ? $release->[1]->{genome_assembly}->{reads_used} : '');
    if ($reads_used =~ /^(\d+)\/(\d+)/) {
        my $numer = commify($1);
        my $denom = commify($2);
        my $expr = $numer . "/" . $denom;
        push (@fields, $expr);
    }
    
    my $line = "|" . join("|", @fields) . "|";
    push (@wiki_combined, $line);
}  


sub create_release
{
    my $species = shift;
    my $strain = shift;
    my $cand_rec = shift;
    my $yaml_rec = shift;
    
    my $release_prefix = get_release_prefix($species, $strain, $cand_rec, $yaml_rec);
    my $species_dir = Assembly::Utils::get_check_record($yaml_rec, ["PE", "species_dir"]);
    my $release_dir = $species_dir . "/release/" . $release_prefix;
    if ($options->{run} and not -e $release_dir) {
        mkpath $release_dir;
    }
    Assembly::Utils::set_check_record($yaml_rec, ["release"], "prefix", $release_prefix);
    Assembly::Utils::set_check_record($yaml_rec, ["release"], "release_dir", $release_dir);
    Assembly::Utils::set_check_record($yaml_rec, ["release"], "version", $cand_rec->{Release_Version});
    unless (release_exists($release_dir, $release_prefix) or $release_dir =~ /^\/release/) {
        mkpath $release_dir;
        my $release_yaml_fname = '';
        if ($options->{testing}) {
            $release_yaml_fname = "/tmp/" . $release_prefix . "_metadata.yml";
        } else {
            $release_yaml_fname = $release_dir . "/" . $release_prefix . "_metadata.yml";
        }
        
        Assembly::Utils::set_check_record($yaml_rec, ["release"], "yaml_file", $release_yaml_fname);
        my $trimraw = $cand_rec->{Trim_Raw};
        my $release = create_release_yaml($species, $strain, $cand_rec, $yaml_rec);
        DumpFile($release_yaml_fname, $release);
        # Optionally print out table of release info to paste into wiki data release
        # section and genome assembly section.
        if ($options->{wiki_release_table}) {
            get_wiki_release_line($release, $release_dir);
        }
        if ($options->{wiki_genome_table}) {
            get_wiki_genome_line($release, $release_dir);
        }
        if ($options->{wiki_combined_table}) {
            get_wiki_combined_line($release, $release_dir);
        }
        # add yml info to release yml record
        # write the release yml record
        # copy over the required files
    }
}

sub create_all_releases
{
    my $release_candidates = shift;
    my $records = shift;
    for my $species (keys %$records) {
        for my $strain (keys %{$records->{$species}->{DNA}}) {
            print "Working on species $species strain $strain\n";
            #my $cand_rec = Assembly::Utils::get_check_record($release_candidates, [$species, $strain]);
            my $cand_rec = $release_candidates->{$species};
            if ($cand_rec) {
                print "Matched species $species\n";
                print "strain is '$strain'\n";
                $cand_rec = $release_candidates->{$species}->{$strain};
                #print "possible strain values: '" . join ("'\t'", (keys %{$release_candidates->{$species}})) . "'\n";
            }
            if ($cand_rec) {
                print "Got release candidate for species $species strain $strain\n";
                my $yaml_rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain]);
                create_release ($species, $strain, $cand_rec, $yaml_rec);
            } else {
                #print "No candidate record for species $species strain $strain\n";
            }
        }
    }
}
        
sub print_wiki_tables
{
    if ($options->{wiki_release_table}) {
        my @release_headers = ("Species", "Strain", "Release Version", "Type", "Link");
        @wiki_release = sort @wiki_release; # I agree, this isn't pretty
        my $fname = $options->{wiki_release_table};
        open (FREL, '>', $fname) or die "Error: could not open output table file $fname\n";
        print FREL "|" . join("|", @release_headers) . "|\n";
        print FREL join("\n", @wiki_release) . "\n";
        close (FREL);
    }
    if ($options->{wiki_genome_table}) {
        my @genome_headers = ("Species", "Strain", "Estimated Genome Length", "N50", "Max Contig",
            "Reads Used", "Total Length");
        @wiki_genome = sort @wiki_genome;
        my $fname = $options->{wiki_genome_table};
        open (FGEN, '>', $fname) or die "Error: could not open output table file $fname\n";
        print FGEN "|" . join("|", @genome_headers) . "|\n";
        print FGEN join("\n", @wiki_genome) . "\n";
        close (FGEN);
    }
    if ($options->{wiki_combined_table}) {
        my @combined_headers = ("Kingdom", "Species", "Strain", "Release", "Type", "PE", "MP-3kb", "MP-8kb", "K-mer", "Link", 
            "Total Length", "Estimated Genome Length", "Number of Contigs", "Min Contig", 
            "Median Contig", "Max Contig", "N50", "Reads Used");
        @wiki_combined = sort @wiki_combined;
        my $fname = $options->{wiki_combined_table};
        open (FCMB, '>', $fname) or die "Error: could not open output table file $fname\n";
        print FCMB "|" . join("|", @combined_headers) . "|\n";
        print FCMB join("\n", @wiki_combined) . "\n";
        close (FCMB);
    }
}

gather_opts;
my $release_candidates = parse_release_table($options->{release_table});
my $records = LoadFile($options->{yaml_in});
create_all_releases($release_candidates, $records);
print_wiki_tables;
DumpFile($options->{yaml_out}, $records);

