#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $options = {};
$options->{stats_outfile} = "assembly_stats.tab";
$options->{stats_outfile_transpose} = "assembly_stats_transpose.tab";
GetOptions( $options,
    'velvet_log|v=s',
    'quast_report_tsv|q=s',
    'contig_file|c=s',
    'print_headers|p',
    'header_only|h',
    'stats_outfile|s=s',
    'stats_outfile_transpose|t=s',
    'kmer|k=s',
    'name|n=s',
    );

my @col_headers = qw(
    assembly_name
    contig_count
    contig_max
    contig_median
    contig_min
    contig_sum
    quast_l50
    quast_largest_contig
    quast_n50
    quast_ns_100kbp
    quast_num_contigs
    quast_total_length
    velvet_cov_depth
    velvet_max_len
    velvet_n50
    velvet_reads_used
    velvet_total_len
    );

# assembly_stats.pl --contig_file "../contigs.fa" --velvet_log "../Log" --quast_tab "quast_results/latest/transposed_report.tsv" --stats_out "assembly_stats.tab" --no_header

sub parse_velvet_log
{
    my $rec = shift;
    my ($n50, $max_len, $total_len, $reads_used, $cov_depth) = ('','','','','');
    my $velvet_log = $options->{velvet_log};
    if (-e $velvet_log) {
        open (FIN, '<', $velvet_log) or die "Error: couldn't open file ${velvet_log}\n";
        while (<FIN>) { 
            if (/n50 of (\d+), max (\d+), total (\d+), using ([0-9\/]+) reads/) {
                ($n50, $max_len, $total_len, $reads_used) = ($1, $2, $3, $4);
                $rec->{velvet_n50} = $n50;
                $rec->{velvet_max_len} = $max_len;
                $rec->{velvet_total_len} = $total_len;
                $rec->{velvet_reads_used} = $reads_used;
            }
            if (/Median coverage depth = ([0-9\.]+)/) {
                $cov_depth = $1;
                $rec->{velvet_cov_depth} = $cov_depth;
            }
        }
        close (FIN);
    }
}

sub parse_quast
{
    my $rec = shift;
    my $quast_report = $options->{quast_report_tsv}; # This should be the report.tsv file.
    if (-s $quast_report) {
        open (QIN, '<', $quast_report) or die "Error: couldn't open file $quast_report.\n";
        while (<QIN>) {
            if (/^# contigs\s+([0-9]+)/) {
                $rec->{quast_num_contigs} = $1;
            }
            if (/^Total length\s+([0-9]+)/) {
                $rec->{quast_total_length} = $1;
            }
            if (/^Largest contig\s+([0-9]+)/) {
                $rec->{quast_largest_contig} = $1;
            }
            if (/^N50\s+([0-9]+)/) {
                $rec->{quast_n50} = $1;
            }
            if (/^L50\s+([0-9]+)/) {
                $rec->{quast_l50} = $1;
            }
            if (/^# N\'s per 100 kbp\s+([0-9\.]+)/) {
                $rec->{quast_ns_100kbp} = $1;
            }
        }
        close (QIN);
    }
}

sub calc_contig_stats
{
    my $rec = shift;
    my $contigs = Bio::SeqIO->new(-file => $options->{contig_file}, -format => "Fasta");

    my @seqlen = ();
    while (my $seq = $contigs->next_seq()) {
        push(@seqlen, length($seq->seq));
    }

    @seqlen = sort { $b <=> $a } @seqlen;
    my $num_contigs = scalar @seqlen;

    my $max_contig_len = $seqlen[0];
    my $min_contig_len = $seqlen[$num_contigs-1];

    my $mid_idx = int ($num_contigs - 0.5) / 2; # 0-base index
    my $median_contig_len = '';
    if ($num_contigs % 2 == 0) {
        $median_contig_len = ($seqlen[$mid_idx] + $seqlen[$mid_idx+1])/2;
    } else {
        $median_contig_len = $seqlen [$mid_idx];
    }
    my $assembly_size = 0;
    for my $len (@seqlen) {
        $assembly_size += $len;
    }
    $rec->{contig_max} = $max_contig_len;
    $rec->{contig_min} = $min_contig_len;
    $rec->{contig_median} = $median_contig_len;
    $rec->{contig_count} = $num_contigs;
    $rec->{contig_sum} = $assembly_size;
}

sub print_stats
{
    my $rec = shift;
    my $outfile = $options->{stats_outfile};
    open (FSTATS, '>', $outfile) or die "Error: couldn't open stats output file $outfile.\n";
    for my $key (sort keys %$rec) {
        my $val = $rec->{$key};
        print FSTATS $key . "\t" . $val . "\n";
    }
    close (FSTATS);
}

sub print_stats_transpose
{
    my $rec = shift;
    my $outfile = $options->{stats_outfile_transpose};
    open (FSTATS_T, '>', $outfile) or die "Error: couldn't open stats output file $outfile.\n";
    my @col_order = sort keys %$rec;
    if ($options->{print_headers}) {
        print FSTATS_T join("\t", (@col_order)) . "\n";
    }
    my @values = map { $rec->{$_} } @col_order;
    print FSTATS_T join("\t", (@values)) . "\n";
    close (FSTATS_T);
}
         
sub report_stats
{
    my $rec = {};
    $rec->{assembly_kmer} = $options->{kmer} if $options->{kmer};
    $rec->{assembly_name} = $options->{name} if $options->{name};
    if ($options->{velvet_log}) {
        parse_velvet_log($rec);
    }
    if ($options->{quast_report_tsv}) {
        parse_quast($rec);
    }
    if ($options->{contig_file}) {
        calc_contig_stats($rec);
    }
    if ($options->{stats_outfile}) {
        print_stats($rec);
    }
    if ($options->{stats_outfile_transpose}) {
        print_stats_transpose($rec);
    }
    if ($options->{header_only}) {
        print join("\t", @col_headers), "\n";
    }
}

report_stats;




           