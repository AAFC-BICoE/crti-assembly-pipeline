#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use File::Path;
use File::Basename;
use Assembly::Utils;

# Given an input spreadsheet with entries for at least the
# following fields:
# Species; strain; sample ID; genome assembly release metadata file; 
# Create the transcriptome assembly within the given base folder
# and generate a release file for the transcript assembly.

# One-liner to start from scratch with all assemblies
# cd /isilon/biodiversity/projects/CRTI/specimen/processing_test2
# find . -maxdepth 6 -path "*RNA/assemblies" -exec find {} -name "trim" \; | xargs -I {} rm -rf {}/

my $options = {};
my @colheaders = ("Species", "Strain", "Sample_ID", "Reference_Strain", "Reference_Metafile");
my $bowtie_path = "/opt/bio/bowtie2/bowtie2-build";
my $tophat_path = "/opt/bio/tophat/bin/tophat";
my $cufflinks_path = "/opt/bio/cufflinks/bin/cufflinks";
my $qsub_path = "/opt/gridengine/bin/lx26-amd64/qsub";
my $qsub_script = "./qsub_script.sh";
my $qsub_nprocs = 10;

my $output_filenames = {
    "bowtie" => [qw(.1.bt2 .2.bt2 .3.bt2 .4.bt2 .fa)],
    "tophat" => [qw(accepted_hits.bam deletions.bed insertions.bed junctions.bed unmapped.bam)],
    "cufflinks" => [qw(genes.fpkm_tracking isoforms.fpkm_tracking skipped.gtf transcripts.gtf)],
    };

sub set_default_opts
{
    my %defaults = qw(
        rna_assembly_table input_data/RnaAssemblyTable.tab
        specimen_dir ../../processing_test2
        yaml_in yaml_files/06b_assembly_setup.yml
        yaml_out yaml_files/15_rna_setup.yml
        create_links 1
        );
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
}

sub check_options
{
        unless ($options->{rna_assembly_table} and $options->{yaml_in} ) {
                die "Usage: $0 -a <RNA assembly table> -i <yaml input file> -o <yaml output file>
                        Optional:
                                --testing
                                --verbose
                                --create_links
                                --specimen_dir <path to specimen/ dir>
                                ";
        }
}

sub gather_options
{
        GetOptions($options,
                'rna_assembly_table|a=s',
                'testing|t',
                'verbose|v',
                'create_links|l',
                'specimen_dir|p=s',
                'yaml_in|i=s',
                'yaml_out|o=s',
                );
        set_default_opts;
        check_options;
}               

sub print_verbose
{
    if ($options->{verbose}) {
        print (@_);
    }
}

# Parse the input table 
sub parse_assembly_table
{
    my $table_recs = {};
    my $fname = ($options->{rna_assembly_table} ? $options->{rna_assembly_table} : '');
    if ($fname and -s $fname) {
        open (FRNA, '<', $fname) or die "Error: couldn't open file $fname\n";
        while (my $line = <FRNA>) {
            my @fields = split (/\t/, $line);
            if (scalar (@fields) == scalar (@colheaders)) {
                my $line_record = {};
                foreach my $i (0..$#colheaders) {
                    # Add hash entry key=column header, value=parsed field value
                    $line_record->{$colheaders[$i]} = $fields[$i];
                }
                my $species_key = Assembly::Utils::format_species_key($line_record->{Species});
                $line_record->{Species} = $species_key; # Reset value to standard access format.
                my $strain_key = Assembly::Utils::format_strain_key($line_record->{Strain});
                $line_record->{Strain} = $strain_key; # Reset value to standard access format.
                my $reference_strain_key = Assembly::Utils::format_strain_key($line_record->{Reference_Strain});
                $line_record->{Reference_Strain} = $reference_strain_key;
                my $sample = $line_record->{"Sample_ID"};
                $table_recs->{$sample} = $line_record;
            } else {
                print_verbose ("Error on line $. of rna assembly table file - incorrect number of cols.\n");
            }
        }
    }
    return $table_recs;
}

# Function below is from Assembly::Qsub.pm. Need to make a static version
# so don't have to create Qsub object to use this.
sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = 1;
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

sub submit_job
{
    my $qsub_cmd = shift;
    my $qsub_jobid = '';
    if ($options->{testing}) {
        print_verbose ("Testing - no submission. Test qsub command is:\n$qsub_cmd\n");
        $qsub_jobid = int(rand(10)) + 1; # Use > 1 because it won't stop anyhting holding on 1 from running
    } else {
        my $submit_str = `$qsub_cmd`;
        $qsub_jobid = get_jobid ($submit_str);
        print_verbose ("Submitted job with jobid $qsub_jobid:\n$qsub_cmd\n");
    }
    return $qsub_jobid;
}

# Note that trailing slash is required on input filename.
# This is because we are using a file prefix in one case.
sub output_files_exist
{
    my $cmd_name = shift;
    my $base_dir = shift;
    my $exist = 1;
    my $file_list = ($output_filenames->{$cmd_name} ? $output_filenames->{$cmd_name} : []);
    foreach my $fname (@$file_list) {
        unless (-e $base_dir . $fname) {
            $exist = 0;
        }
    }
    return $exist;
}

# Check if output files exist; if not, submit job and return jobid.
sub check_submit_cmd
{
    my $cmd_name = shift;
    my $base_dir = shift;
    my $qsub_cmd = shift;
    my $qsub_jobid = 1;
    if (output_files_exist ($cmd_name, $base_dir)) {
        print_verbose ("Output $cmd_name files already exist - not issuing the following command:\n$qsub_cmd\n");
    } else {
        $qsub_jobid = submit_job ($qsub_cmd);
    }
    return $qsub_jobid;
}

# Two funcs below could be a mini-package to handle linking of files/creating of dirs?
sub create_symlink
{
    my $infile = shift;
    my $outfile = shift;
    if ($options->{create_links}) {
        unless ($options->{testing}) {
            if (-e $infile and -s $infile) {
                my $dirname = dirname ($outfile);
                if (-d $dirname) {
                    unless (-e $outfile) {
                        symlink ($infile, $outfile);
                    }
                }
            }
        }
    }
}

sub create_dir
{
    my $dirname = shift;
    unless ($options->{testing}) {
        unless (-e $dirname) {
            print_verbose ("Creating directory $dirname\n");
            mkpath $dirname;
        }
    }
}


# Given release metadata file, link the original contigs.fa (not release copy)
# to the rna assembly working directory
sub get_genome_assembly
{
    my $metadata_file = shift;
    
    my $genome_prefix = '';
    my $metadata_base = basename ($metadata_file);
    if ($metadata_base =~ /(.*)_metadata\.yml/) {
        $genome_prefix = $1;
    }
    
    my $meta_record = LoadFile ($metadata_file);
    my @p = @$meta_record;
    my $pi = $p[$#p];
    my $pip = Assembly::Utils::get_check_record ($pi, ["pipeline"]);

    my $last_cmd_idx = scalar @$pip - 1;
    my $cmd_rec = $pip->[$last_cmd_idx];
    my $genome_file = $cmd_rec->{"release"}->[0]->{"input_file"};
    return ($genome_file, $genome_prefix);
}

# Create rna assembly dirname; add info to rna assembly yaml rec; create dir
sub add_assembly_dir
{
    my ($yrec, $asm_rec, $strain, $sample, $trimraw) = @_;
    # my $species_abbr = Assembly::Utils::get_check_record($yrec, ["species_abbr"]);
    my $species_dir = Assembly::Utils::get_check_record($yrec, ["species_dir"]);
    print "strain $strain Got species dir $species_dir\n";
    my $dirname = $species_dir . "/RNA/assemblies/" . $strain . "/" . $sample . "/" . $trimraw . "/";
    Assembly::Utils::set_check_record($asm_rec, [], "rna_assembly_dir", $dirname);
    #create_dir ($dirname);
    mkpath $dirname;
}

# Find raw reads; add info to rna assembly yaml rec; create links in rna assembly dir
sub add_read_files
{
    my $yrec = shift;
    my $asm_rec = shift;
    my $trimraw = shift;
    my $trdata = $trimraw . "data";
    my $r1data = Assembly::Utils::get_check_record($yrec, ["R1", $trdata]);
    my $r2data = Assembly::Utils::get_check_record($yrec, ["R2", $trdata]);
    
    # The first two steps below shouldn't be necessary in future ... all input yaml 
    # raw files will have the .gz extension ... eventually
    if ($r1data !~ /\.gz\s*$/ and $r1data !~ /trim/) { $r1data .= ".gz"; }
    if ($r2data !~ /\.gz\s*$/ and $r2data !~ /trim/) { $r2data .= ".gz"; }
    $r1data = `readlink -e $r1data`; 
    chomp $r1data;
    $r2data = `readlink -e $r2data`; 
    chomp $r2data;
    
    Assembly::Utils::set_check_record($asm_rec, [], "r1data", $r1data);
    Assembly::Utils::set_check_record($asm_rec, [], "r2data", $r2data);
    
    # Link in the files
    my $dirname = Assembly::Utils::get_check_record($asm_rec, ["rna_assembly_dir"]);
    my $r1outfile = $dirname . "/" . basename($r1data);
    my $r2outfile = $dirname . "/" . basename($r2data);
    create_symlink ($r1data, $r1outfile);
    create_symlink ($r2data, $r2outfile);
}

# Find existing genome assembly that was released; add info to rna assembly yaml rec; link in non-release copy.
sub add_genome_file
{
    my $asm_rec = shift;
    my $genome_infile = shift;
    my $genome_prefix = shift;
    $genome_infile = `readlink -e $genome_infile`; 
    chomp $genome_infile;
    my $dirname = Assembly::Utils::get_check_record($asm_rec, ["rna_assembly_dir"]);
    $genome_prefix = $dirname . "/" . $genome_prefix;
    my $genome_outfile = $genome_prefix . ".fa"; # "_" . basename($genome_infile)
    Assembly::Utils::set_check_record($asm_rec, [], "genome_infile", $genome_infile);
    Assembly::Utils::set_check_record($asm_rec, [], "genome_prefix", $genome_prefix);;
    create_symlink ($genome_infile, $genome_outfile);
}    

# Add full file info for genome file, output dir, and r1, r2 read files
# to yaml record for future access. Create rna assembly dirs and link in required files.
sub add_file_info
{
    my ($yaml_recs, $species, $strain, $sample, $reference_strain, $trimraw, $genome_file, $genome_prefix) = @_;
    my $yrec = Assembly::Utils::get_check_record($yaml_recs, [$species, "RNA", $strain, $sample]);
    unless ($yrec) { print "ERROR: species '$species', strain '$strain', sample '$sample'\n"; }
    Assembly::Utils::set_check_record($yaml_recs, [$species, "RNA", $strain, $sample, "RNA_assembly"], $trimraw, {});
    my $asm_rec = Assembly::Utils::get_check_record($yaml_recs, [$species, "RNA", $strain, $sample, "RNA_assembly", $trimraw]);
    add_assembly_dir ($yrec, $asm_rec, $strain, $sample, $trimraw);
    add_genome_file($asm_rec, $genome_file, $genome_prefix);
    add_read_files ($yrec, $asm_rec, $trimraw);
    return $asm_rec;
}

sub run_bowtie
{
    my $asm_rec = shift;
    my $genome_file = Assembly::Utils::get_check_record ($asm_rec, ["genome_infile"]);
    my $genome_prefix = Assembly::Utils::get_check_record ($asm_rec, ["genome_prefix"]);
    my $assembly_dir = Assembly::Utils::get_check_record ($asm_rec, ["rna_assembly_dir"]);
    my $bowtie_cmd = "$bowtie_path $genome_file $genome_prefix";
    my $bowtie_qsub_cmd = $qsub_path . " -N bowtie2-build " . $qsub_script . " '" . $bowtie_cmd . "'";
    my $bowtie_qsub_jobid = check_submit_cmd ("bowtie", $genome_prefix, $bowtie_qsub_cmd);
    Assembly::Utils::set_check_record($asm_rec, [], "bowtie_cmd", $bowtie_cmd);
    Assembly::Utils::set_check_record($asm_rec, [], "bowtie_qsub_cmd", $bowtie_qsub_cmd);
    Assembly::Utils::set_check_record($asm_rec, [], "bowtie_qsub_jobid", $bowtie_qsub_jobid);
    # Line below should be changed.
    my $bowtie_version = `ssh biocomp-0-0 '$bowtie_path --version'`;
    Assembly::Utils::set_check_record($asm_rec, [], "bowtie_version", $bowtie_version);
    return $bowtie_qsub_jobid;
}

sub run_tophat
{
    my $asm_rec = shift;
    my $bowtie_jobid = shift;
    my $sample = shift;
    my $assembly_dir = Assembly::Utils::get_check_record($asm_rec, ["rna_assembly_dir"]);
    my $tophat_dir = $assembly_dir . "/" . $sample . "_tophat";
    create_dir ($tophat_dir);
    my $genome_prefix = Assembly::Utils::get_check_record($asm_rec, ["genome_prefix"]);
    my $r1data = Assembly::Utils::get_check_record($asm_rec, ["r1data"]);
    my $r2data = Assembly::Utils::get_check_record($asm_rec, ["r2data"]);
    my $tophat_cmd = $tophat_path . " -p " . $qsub_nprocs . " -o " . $tophat_dir . " " . 
        $genome_prefix . " " . $r1data . " " . $r2data;
    my $tophat_qsub_cmd = $qsub_path . " -N tophat -hold_jid " . $bowtie_jobid . 
        " -pe smp " . $qsub_nprocs . " " . $qsub_script . " '" . $tophat_cmd . "'";
    my $tophat_qsub_jobid = check_submit_cmd("tophat", $tophat_dir . "/", $tophat_qsub_cmd);
    Assembly::Utils::set_check_record($asm_rec, [], "tophat_dir", $tophat_dir);
    Assembly::Utils::set_check_record($asm_rec, [], "tophat_cmd", $tophat_cmd);
    Assembly::Utils::set_check_record($asm_rec, [], "tophat_qsub_cmd", $tophat_qsub_cmd);
    Assembly::Utils::set_check_record($asm_rec, [], "tophat_qsub_jobid", $tophat_qsub_jobid);
    my $tophat_version = `ssh biocomp-0-0 '$tophat_path --version'`;
    Assembly::Utils::set_check_record($asm_rec, [], "tophat_version", $tophat_version);
    return $tophat_qsub_jobid;
}

sub run_cufflinks
{
    my $asm_rec = shift;
    my $tophat_jobid = shift;
    my $sample = shift;
    my $assembly_dir = Assembly::Utils::get_check_record($asm_rec, ["rna_assembly_dir"]);
    my $cufflinks_dir = $assembly_dir . "/" . $sample . "_cufflinks";
    create_dir ($cufflinks_dir);
    my $tophat_dir = Assembly::Utils::get_check_record($asm_rec, ["tophat_dir"]);
    my $accepted_hits_file = $tophat_dir . "/accepted_hits.bam";
    my $cufflinks_cmd = $cufflinks_path . " -p " . $qsub_nprocs . " -o " . 
        $cufflinks_dir . " " . $accepted_hits_file;
    my $cufflinks_qsub_cmd = $qsub_path . " -N cufflinks -hold_jid " . $tophat_jobid . 
        " -pe smp " . $qsub_nprocs . " " . $qsub_script . " '" . $cufflinks_cmd . "'";
    my $cufflinks_qsub_jobid = check_submit_cmd ("cufflinks", $cufflinks_dir . "/", $cufflinks_qsub_cmd);
    Assembly::Utils::set_check_record ($asm_rec, [], "cufflinks_dir", $cufflinks_dir);
    Assembly::Utils::set_check_record ($asm_rec, [], "cufflinks_cmd", $cufflinks_cmd);
    Assembly::Utils::set_check_record ($asm_rec, [], "cufflinks_qsub_cmd", $cufflinks_qsub_cmd);
    Assembly::Utils::set_check_record ($asm_rec, [], "cufflinks_qsub_jobid", $cufflinks_qsub_jobid);
    my $cufflinks_version = `ssh biocomp-0-0 '$cufflinks_path 2>&1 | egrep "^cufflinks"'`;
    Assembly::Utils::set_check_record($asm_rec, [], "cufflinks_version", $cufflinks_version);
    return $cufflinks_qsub_jobid;
} 

sub create_assembly
{
    my ($yaml_recs, $species, $strain, $sample, $reference_strain, $trimraw, $genome_file, $genome_prefix) = @_;
    my $asm_rec = add_file_info ($yaml_recs, $species, $strain, $sample, $reference_strain, $trimraw, $genome_file, $genome_prefix);
    
    my $bowtie_jobid = run_bowtie ($asm_rec);
    my $tophat_jobid = run_tophat ($asm_rec, $bowtie_jobid, $sample);
    my $cufflinks_jobid = run_cufflinks ($asm_rec, $tophat_jobid, $sample);
}

sub run_rna_assembly
{
    my $yaml_recs = shift;
    my $table_recs = shift;
    for my $sample (keys %$table_recs) {
    #for my $sample ("S001278") {
        # Get the parameters of interest.
        my $species = $table_recs->{$sample}->{"Species"};
        my $strain = $table_recs->{$sample}->{"Strain"};
        my $reference_strain = $table_recs->{$sample}->{"Reference_Strain"};
        my ($genome_file, $genome_prefix) = get_genome_assembly ($table_recs->{$sample}->{"Reference_Metafile"});
        next unless (-s $genome_file);
        for my $trimraw (qw(trim)) {
            create_assembly ($yaml_recs, $species, $strain, $sample, $reference_strain, $trimraw, $genome_file, $genome_prefix);
        }
    }
}

gather_options;
my $yaml_recs = LoadFile($options->{yaml_in});
my $table_recs = parse_assembly_table;
run_rna_assembly ($yaml_recs, $table_recs);
DumpFile ($options->{yaml_out}, $yaml_recs);
