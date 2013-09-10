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
        raw 1
        verbose 0
        run 0
        velvetk_infile input_data/VelvetKBest.tab
        velvetk_outfile output_files/VelvetKBestOut.tab
        velvet_stats_file output_files/VelvetStats.tab
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
                --velvet_stats_file <filename>
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
        'velvet_stats_file=s',
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

# Step 1 - go through the existing VelvetK.tab file and pull out kmers to records
sub read_velvetk_table
{
    my $records = shift;
    my $fname = ($options->{velvetk_infile} ? $options->{velvetk_infile} : '');
    unless ($fname) { return; }
    my @headers = qw (Species Strain PE PER MP MP3 MP8 trim_raw velvetk_kmer);
    open (FIN, '<', $fname) or die "Error: couldn't open file $fname\n";
    <FIN>; # skip headers
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
        if ($rec{velvetk_kmer}) {
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer", $rec{velvetk_kmer});
        }
    }
    close (FIN);
}
        
# Step 2 go through the VelvetStats.tab file and pull out kmers to records
sub read_prev_best_table
{
    my $records = shift;
    my $fname = ($options->{velvet_stats_file} ? $options->{velvet_stats_file} : '');
    my @headers = qw (Species Strain trim_raw velvetk_kmer prev_best_kmer prev_best_n50 missing_kmers);
    open (FIN, '<', $fname) or die "Error: couldn't open file $fname\n";
    <FIN>; # skip headers
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
        if ($rec{prev_best_kmer}) {
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "prev_best_n50_kmer", $rec{prev_best_kmer});
        }
    }
    close (FIN);
}

sub submit_velvetk_cmd
{
    my $records = shift;
    my $species = shift;
    my $strain = shift;
    my $trimraw = shift;
    my $trimdata = $trimraw . "data";

    my $genome_size = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "related_genome_length", "RG_Est_Genome_Length"]);
    if ($genome_size) { 
        my $vk_cmd = $velvetk_bin . " --size " . $genome_size . " --best ";
        for my $sample_type (qw(PE PER MP MP3 MP8)) {
            my $r1data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R1", $trimdata]);
            my $r2data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R2", $trimdata]);
            $vk_cmd .= $r1data . " " . $r2data . " ";
        }
        my $vk_qsub_cmd = get_qsub_cmd($vk_cmd);
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_cmd", $vk_cmd);
        Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_qsub_cmd", $vk_qsub_cmd);
        
        if ($options->{run} and $vk_cmd) {
            print "Running command $vk_cmd\n";
            my $best = `$vk_cmd`;
            chomp $best;
            Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer", $best);
            #push (@newbest, [$species, "DNA", $strain, $best]);
            
            # qsub version:
            # have to pull the jobid, then submit another job that holds on this one
            # and then goes through the qsub output files for this jobid and pulls out
            # the best kmer from that.
            # my $qsub_str = `$vk_qsub_cmd`;
            # my $jobid = get_jobid($qsub_str);
            # get_best_from_outfile_in_cwd(hold on $jobid)
        }
    } else {
        print "Couldn't get genome size for species $species strain $strain\n";
    }
}

sub check_velvetk_records
{
    my $records = shift;
    for my $species (keys %$records) {
        for my $strain (keys %{$records->{$species}->{DNA}}) {
            for my $trimraw (qw(trim raw)) {
                my $velvetk_best = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_best_kmer");
                # my $prev_best_kmer = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "prev_best_n50_kmer");
                unless ($velvetk_best) {
                    submit_velvetk_cmd($records, $species, $strain, $trimraw);
                }
            }
        }
    }
}    


# The old function
sub run_velvetk
{
    my $records = shift;
    my @newbest = ();
    my $fname = ($options->{velvetk_infile} ? $options->{velvetk_infile} : '');
    if ($fname) {
        open (FIN, '<', $fname) or die "Error: couldn't open file $fname\n";
        <FIN>;
        my @headers = qw (Species Strain PE PER MP MP3 MP8 trim_raw velvetk_kmer);
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
            if ($genome_size and not $rec{velvetk_kmer}) { 
                
                my $vk_cmd = $velvetk_bin . " --size " . $genome_size . " --best ";
                for my $sample_type (qw(PE PER MP MP3 MP8)) {
                    #my $hr = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type]);
                    my $r1data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R1", $trimdata]);
                    my $r2data = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "R2", $trimdata]);
                    
                    $vk_cmd .= $r1data . " " . $r2data . " ";
                }
                #print "got velvetk_cmd: $vk_cmd\n";
                my $vk_qsub_cmd = get_qsub_cmd($vk_cmd);
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_cmd", $vk_cmd);
                Assembly::Utils::set_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw], "velvetk_qsub_cmd", $vk_qsub_cmd);
                #print "got qsub_cmd: $vk_qsub_cmd\n";
                
                if ($options->{run} and $vk_cmd) {
                    print "Running command $vk_cmd\n";
                    my $best = `$vk_cmd`;
                    chomp $best;
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

sub write_vk_kmers
{
    my $records = shift;
    my $sample = '';
    my $fname = ($options->{velvetk_outfile} ? $options->{velvetk_outfile} : '');
    if ($fname) {
        open (FTAB, '>', $fname) or die "Error: couldn't open file $fname\n";
        print FTAB join("\t", qw(Species Strain PE PER MP MP3 MP8 Trim/Raw VK_Best_Kmer)) . "\n";
        for my $species (keys %$records) {
            my $ss = $records->{$species}->{DNA};
            for my $strain (keys %$ss) {
                my $out_str = $species . "\t" . $strain . "\t";
                for my $sample_type (qw(PE PER MP MP3 MP8)) {
                    my $sample_id = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, $sample_type, "sample"]);
                    $out_str .= $sample_id . "\t";
                }
                for my $trimraw (qw(trim raw)) {
                    my $vk_best = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "velvetk_best_kmer"]);
                    print FTAB $out_str . $trimraw . "\t" . $vk_best . "\n";
                }
            }
        }
        close (FTAB);
    }
}

gather_opts;
my $records = LoadFile($options->{yaml_in});
#run_velvetk($records);
read_velvetk_table ($records);
read_prev_best_table ($records);
check_velvetk_records ($records);
write_vk_kmers($records);
DumpFile($options->{yaml_out}, $records);









