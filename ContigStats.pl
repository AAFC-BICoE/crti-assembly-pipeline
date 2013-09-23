#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);
use Bio::SeqIO;
use Assembly::Utils;

my $options = {};
my @contig_col_headers = qw(Species Strain Trim_Raw Kmer Num_Contigs Min_Contig_Len Median_Contig_Len);

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/12_velvet_stats.yml
        yaml_out yaml_files/13_contig_stats.yml
        contig_stats_in input_data/ContigStats.tab
        contig_stats_out output_files/ContigStatsOut.tab
        release_table input_data/Release.tab
        verbose 1
        );
    for my $kdef (keys %defaults) {
        $options->{$kdef} = $defaults{$kdef} unless $options->{$kdef};
    }
}

sub check_opts
{
    unless ($options->{yaml_in}) {
        die "Usage: $0 -i <input yaml file> -o <output yaml file>
            Optional:
                --contig_stats_in <input filename>
                --contig_stats_out <output filename>
                --verbose
                --stats_outfile <filename>
                ";
    }
}

sub gather_opts
{
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'contig_stats_in=s',
        'contig_stats_out=s',
        'verbose',
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

# Get an array of all the kmer values
sub get_kmer_range
{
    my $rec = shift;
    my $tr = shift;
    my $kmin = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "min_kmer"]);
    my $kmax = Assembly::Utils::get_check_record($rec, ["velvet", $tr, "max_kmer"]);
    my $krange = [];
    if ($kmin =~ /^\d+$/ and $kmax =~ /^\d+$/) {
        for (my $kmer = $kmin; $kmer <= $kmax; $kmer = $kmer + 2) {
            push (@$krange, $kmer);
        }
    }
    return $krange;
}


# Open current contig stats table and update input assembly records with relevant stats.
sub parse_input_stats
{
    my $records = shift;
    my $fname = ($options->{contig_stats_in} ? $options->{contig_stats_in} : '');
    if ($fname) {
        open (CSIN, '<', $fname) or die "Error: couldn't open input contig stats file $fname\n";
        while (my $line = <CSIN>) {
            chomp $line;
            my @fields = split (/\t/, $line);
            my $num_headers = scalar @contig_col_headers;
            my $num_fields = scalar @fields;
            if ($num_headers == $num_fields) {
                my %rec = ();
                my $i = 0;
                foreach my $header (@contig_col_headers) {
                    $rec{$header} = $fields[$i];
                    $i++;
                }
                my $species = Assembly::Utils::format_species_key($rec{Species});
                my $strain = Assembly::Utils::format_strain_key($rec{Strain});
                my $trimraw = $rec{Trim_Raw};
                my $kmer_rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $rec{Kmer}]);
                if ($kmer_rec) {
                    Assembly::Utils::set_check_record($kmer_rec, [], "num_contigs", $rec{Num_Contigs});
                    Assembly::Utils::set_check_record($kmer_rec, [], "min_contig_len", $rec{Min_Contig_Len});
                    Assembly::Utils::set_check_record($kmer_rec, [], "median_contig_len", $rec{Median_Contig_Len});
                } else {
                    print_verbose "No record found for contig stats line:\n$line\n";
                }
            } else {
                print_verbose "Line $. of contig stats file has incorrect number of fields.\n";
            }
        }
        close (CSIN);
    }
}

sub calc_contig_stats
{
    my $contigs_filename = shift;
    my $contigs = Bio::SeqIO->new(-file => $contigs_filename, -format => "Fasta");

    my @seqlen = ();
    while (my $seq = $contigs->next_seq()) {
        push(@seqlen, length($seq->seq));
    }

    @seqlen = sort { $b <=> $a } @seqlen;
    my $num_contigs = scalar @seqlen;

    my $min_contig_len = $seqlen[$num_contigs-1];

    my $mid_idx = int ($num_contigs - 0.5) / 2; # 0-base index
    my $median_contig_len = '';
    if ($num_contigs % 2 == 0) {
        $median_contig_len = ($seqlen[$mid_idx] + $seqlen[$mid_idx+1])/2;
    } else {
        $median_contig_len = $seqlen [$mid_idx];
    }
    
    return ($num_contigs, $min_contig_len, $median_contig_len);
}

sub get_contig_recs
{
    my $kmer_rec = shift;
    my $num_contigs = Assembly::Utils::get_check_record($kmer_rec, ["num_contigs"]);
    my $min_contig_len = Assembly::Utils::get_check_record($kmer_rec, ["min_contig_len"]);
    my $median_contig_len = Assembly::Utils::get_check_record($kmer_rec, ["median_contig_len"]);
    return ($num_contigs, $min_contig_len, $median_contig_len);
}

sub set_contig_recs
{
    my $kmer_rec = shift;
    my $num_contigs = shift;
    my $min_contig_len = shift;
    my $median_contig_len = shift;
    if ($num_contigs and $min_contig_len and $median_contig_len) {
        Assembly::Utils::set_check_record($kmer_rec, [], "num_contigs", $num_contigs);
        Assembly::Utils::set_check_record($kmer_rec, [], "min_contig_len", $min_contig_len);
        Assembly::Utils::set_check_record($kmer_rec, [], "median_contig", $median_contig_len);
    } else {
        print_verbose "Warning: contigs stats were not properly set - discarding\n";
    }
}   

sub get_stats
{
    my $kmer_rec = shift;
    my ($num_contigs, $min_contig_len, $median_contig_len) = get_contig_recs($kmer_rec);
    unless ($num_contigs and $min_contig_len and $median_contig_len) {
        print "Couldnt get rec !\n";
        my $kmer_dir = Assembly::Utils::get_check_record($kmer_rec, ["kmer_dir"]);
        my $contigs_file = $kmer_dir . "/contigs.fa";
        if (-e $contigs_file and -s $contigs_file) { # file exists and has nonzero size
            ($num_contigs, $min_contig_len, $median_contig_len) = calc_contig_stats ($contigs_file);
            print "Got num contigs $num_contigs min contig $min_contig_len, median contig $median_contig_len\n";
            set_contig_recs($kmer_rec, $num_contigs, $min_contig_len, $median_contig_len);
        }
    }
    return ($num_contigs, $min_contig_len, $median_contig_len);
}

sub all_contig_stats
{
    my $records = shift;
    my $output_stats = [];
    my $species = "Globodera_pallida";
    my $strain = "A2_Cyprus_pathotype_Pa23";
    my $trimraw = "trim";
    my $kmer = 43;
    for my $species (keys %$records) {
        for my $strain (keys %{$records->{$species}->{DNA}}) {
            my $strain_rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain]);
            for my $trimraw (qw(trim raw)) {
                my $kmer_range = get_kmer_range($strain_rec, $trimraw);
                for my $kmer (@$kmer_range) {
                    my $kmer_rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $trimraw, "kmer", $kmer]);
                    if ($kmer_rec) {
                        print_verbose "Working on contigs file for species $species, strain $strain, trim/raw $trimraw, kmer $kmer\n";
                        my ($num_contigs, $min_contig, $median_contig) = get_stats($kmer_rec);
                        if ($options->{contig_stats_out}) {
                            push (@$output_stats, [$species, $strain, $trimraw, $kmer, $num_contigs, $min_contig, $median_contig]);
                        }
                    }
                }
            }
        }
    }
    return $output_stats; 
}  

sub release_contig_stats
{
    my $records = shift;
    my $fname = ($options->{release_table} ? $options->{release_table} : '');
    my $output_stats = [];
    my @release_headers = qw (Species Strain Sequencing_Types Trim_Raw Release_Version Est_Genome_Size Best_Kmer Best_N50);
    unless ($fname and -e $fname) { return; }
    open (RTAB, '<', $fname) or die "Error: couldn't open release table file $fname\n";
    <RTAB>; # skip header
    while (my $line = <RTAB>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        my $num_headers = scalar @release_headers;
        my $num_fields = scalar @fields;
        my %rec = ();
        for (my $i=0; $i<$num_headers; $i++) { 
            my $key = $release_headers[$i];
            my $val = ($fields[$i] ? $fields[$i] : '');
            $rec{$key} = $val; 
        }
        my $species = Assembly::Utils::format_species_key($rec{Species});
        my $strain = Assembly::Utils::format_strain_key($rec{Strain});
        my $kmer_rec = Assembly::Utils::get_check_record($records, [$species, "DNA", $strain, "velvet", $rec{Trim_Raw}, "kmer", $rec{Best_Kmer}]);
        if ($kmer_rec) {
            print_verbose "Working on contigs file for species $species, strain $strain, trim/raw " . $rec{Trim_Raw} . ", kmer " . $rec{Best_Kmer} . "\n";
            my ($num_contigs, $min_contig, $median_contig) = get_stats($kmer_rec);
            if ($options->{contig_stats_out}) {
                push (@$output_stats, [$species, $strain, $rec{Trim_Raw}, $rec{Best_Kmer}, $num_contigs, $min_contig, $median_contig]);
            }
        }
    }
    close (RTAB);
    return $output_stats;
}

sub write_contig_stats
{
    my $output_stats = shift;
    my $fname = ($options->{contig_stats_out} ? $options->{contig_stats_out} : '');
    if ($fname) {
        open (CSOUT, '>', $fname) or die "Error: could not open output contig stats file $fname\n";
        foreach my $stats_arr (@$output_stats) {
            print CSOUT join("\t", (@$stats_arr)) . "\n";
        }
        close (CSOUT);
    }
}
        
gather_opts;
my $records = LoadFile($options->{yaml_in});
parse_input_stats($records);
#my $output_stats = all_contig_stats($records);
my $output_stats = release_contig_stats($records);
write_contig_stats($output_stats);
DumpFile($options->{yaml_out}, $records);














