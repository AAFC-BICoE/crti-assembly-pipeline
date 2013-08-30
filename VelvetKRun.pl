#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw (LoadFile DumpFile);
use Assembly::Utils;
use Cwd;

my $options = {};
my $cdir = getcwd;
my $velvetk_bin = $cdir . "/velvetk.pl";
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/09_velvet_advisor.yml
        yaml_out yaml_files/10_velvetk.yml
        trim 1
        raw 0
        verbose 0
        run 0
        velvetk_infile input_data/AssemblySetup.tab
        velvetk_outfile output_files/VelvetKBestOut.tab
        qsub_script qsub_script.sh
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --verbose
                --sample_list <ID1,ID2,...,IDX (no spaces)>
                --trim
                --raw
                --run
                --velvetk_infile <filename>
                --velvetk_outfile <filename>
                --qsub_script <filename>
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'verbose',
        'sample_list=s',
        'trim',
        'raw',
        'velvetk_infile=s',
        'velvetk_outfile=s',
        'run',
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
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " -N velvetk " . 
            " " . $options->{qsub_script} . " '" . $cmd . "'";
    return $qsub_cmd;
}

sub get_jobid
{
    my $qsub_str = shift;
    my $hold_jobid = '';
    if ($qsub_str =~ /Your job[^\s]*\s(\d+)[\.\s]/) {
        $hold_jobid = $1;
    }
    return $hold_jobid;
}

sub run_velvetk
{
    my $records = shift;
    my @newbest = ();
    my $fname = ($options->{velvetk_infile} ? $options->{velvetk_infile} : '');
    if ($fname) {
        open (FIN, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FIN>;
        my @headers = qw (Species Strain PE PER MP MP3 MP8 trim_raw velvetk_kmer velveth_done velvetg_done best_kmer best_n50);
        while (my $line = <FIN>) {
            chomp $line;
            my @fields = split (/\t/, $line);
            my $num_headers = scalar @headers;
            my $num_fields = scalar @fields;
            my %rec;
            for (my $i=0; $i<$num_headers; $i++) { 
                my $key = $headers[$i];
                my $val = ($fields[$i] ? $fields[$i] : '');
                $rec{$key} = $val; 
            }
            my $species = Assembly::Utils::format_species_key($rec{Species});
            my $strain = Assembly::Utils::format_strain_key($rec{Strain});
            my $trimraw = $rec{trim_raw};
            my $trimdata = $rec{trim_raw} . "data";
            my $genome_size = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "related_genome_length", "RG_Est_Genome_Length"]);
            unless ($genome_size) {
                print "Couldn't get genome size for species $species strain $strain\n";
            }
            
            if ($genome_size ) { 
                
                my $vk_cmd = $velvetk_bin . " --size " . $genome_size . " --best ";
                for my $sample_type (qw(PE PER MP MP3 MP8)) {
                    #my $hr = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type]);
                    my $r1data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R1", $trimdata]);
                    my $r2data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R2", $trimdata]);
                    
                    $vk_cmd .= $r1data . " " . $r2data . " ";
                }
                print "got velvetk_cmd: $vk_cmd\n";
                my $vk_qsub_cmd = get_qsub_cmd($vk_cmd);
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_cmd", $vk_cmd);
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_qsub_cmd", $vk_qsub_cmd);
                print "got qsub_cmd: $vk_qsub_cmd\n";
                
                if ($options->{run} and $vk_cmd) {
                    my $best = `$vk_cmd`;
                    chomp $best;
                    print "velvetk.pl found best kmer: " . $best . "\n";
                    Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer", $best);
                    push (@newbest, [$species, "DNA", $strain, $best]);
                    
                    # qsub version:
                    # have to pull the jobid, then submit another job that holds on this one
                    # and then goes through the qsub output files for this jobid and pulls out
                    # the best kmer from that.
                    # my $qsub_str = `$vk_qsub_cmd`;
                    # my $jobid = get_jobid($qsub_str);
                    # get_best_from_outfile_in_cwd(hold on $jobid)
                }
            } elsif ($rec{velvetk_kmer}) {
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer", $rec{velvetk_kmer});
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_cmd", "");
            }
        }
    }
    close (FIN);
    foreach my $arr (@newbest) {
        print join("\t", @$arr) . "\n";
    }
}

sub write_new_kmers
{
    my $records = shift;
    my $sample = '';
    my $fname = ($options->{velvetk_outfile} ? $options->{velvetk_outfile} : '');
    if ($fname) {
        open (FTAB, '>', $fname) or die "Error: couldn't open file $fname\n";
        print FTAB join("\t", qw(Species Strain Trim/Raw VK_Best_Kmer VK_Command)) . "\n";
        for my $species (keys %$records) {
            my $ss = $records->{$species}->{DNA};
            for my $strain (keys %$ss) {
                for my $trimraw (qw(trim raw)) {
                    my $vk_cmd = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_cmd"]);
                    my $vk_best = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_best_kmer"]);
                    print FTAB join("\t", ($species, $strain,$trimraw, $vk_best, $vk_cmd)) . "\n";
                }
            }
        }
    }
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
run_velvetk($records);
DumpFile($options->{yaml_out}, $records);









