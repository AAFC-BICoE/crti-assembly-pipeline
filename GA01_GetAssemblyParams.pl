#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(LoadFile DumpFile);

# Gather and calculate as many stats as possible 
# that will be potentially useful (even if only as 
# error checks) for genome assembly with velvet.

my $options = {};
my $records = {};

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file> -s <illumina summary table file>
            Optional:
                --sample <sample ID>
                --no_raw_stats
                --input_read_table <input table filename>
                --output_read_table <output table filename>
                --verbose
                ";
    }         
}

sub gather_opts
{
    GetOptions($options,
        'yaml_in|i=s',
        'yaml_out|o=s',
        'illumina_summary_table|s=s',
        'sample=s',
        'no_raw_stats',
        'input_read_table',
        'output_read_table',
        'verbose',
        
        );
    check_opts;
}

sub get_illumina_stats
{
    my $fname = $options->{illumina_summary_table};
    open (FIS, '<', $fname) or die "Error: couldn't open file $fname\n";
    while (my $line = <FIS>) {
        chomp $line;
        my @fields = split(/\t/, $line); 
        my $sample = ($fields[13] ? $fields[13] : '');
        my $numreads = ($fields[17] ? $fields[17] : '');
        $sample =~ s/^(\s+)|(\s+$)//g;
        $numreads =~ s/^(\s+)|(\s+$)//g;
        if ($sample and $sample =~ /[A-Z0-9]+/ and $numreads and $numreads =~ /^\d+/) {
            if (defined $records->{$sample}) {
                my $rec = $records->{$sample};
                $rec->{illumina_numreads} = $numreads;
            }
        }
    }
    close (FIS);
}

sub record_read_info
{
    my $fref = shift;
    if (scalar @$fref == 5) {
        my ($sample, $trdata, $direction, $read_length, $num_reads) = @$fref;
        if (defined ($records->{$sample})) {
            my $ref = $records->{$sample};
            for my $key ("data_stats", $direction, $trdata) {
                $ref->{$key} = {} unless (defined ($ref->{$key}));
                $ref = $ref->{$key};
            }
            $ref->{read_length} = $read_length;
            $ref->{num_reads} = $num_reads;
        } elsif ($options->{verbose}) {
            print "No yaml records for sample $sample - referenced in the input read info table.\n";
        }
    } elsif ($options->{verbose}) {
        print "Invalid number of entries in input read table line $.\n";
    }
}
    

sub read_table_stats
{
    my $intable = ($options->{input_read_table} ? $options->{input_read_table} : '');
    if ($intable) {
        open (FINTAB, '<', $intable) or die "Error: couldn't open input read table $intable\n";
        <FINTAB>; # Skip first line, containing column headers.
        while (my $line = <FINTAB>) {
            chomp $line;
            my @fields = split(/\t/, $line); 
            record_read_info(\@fields);
        }
        close (FINTAB);
    }
}
        
# Given a raw or trimmed .fq file, return
# the number of reads and the read length.
sub get_read_info
{
    my $fq_infile = shift;
    open (FQIN, '<', $fq_infile);
    my ($num_reads, $read_length) = (0, 0);
    if ($options->{verbose}) { print "Determining read length/num reads from file $fq_infile\n"; }
    while (my $line = <FQIN>) {
        if ($. == 2) {
            chomp $line;
            $read_length = length($line);
        }
        if ($line =~ /^\@HWI/) {
            $num_reads++;
        }
    }
    close (FQIN);
    print join("\t", ($fq_infile, $read_length, $num_reads)) . "\n";
    return ($read_length, $num_reads);
}

sub calc_stats
{
    my $rec = shift;
    $rec->{data_stats} = {} unless (defined ($rec->{data_stats}));
    for my $direction (qw(R1 R2)) {
        $rec->{data_stats}->{$direction} = {} unless (defined ($rec->{data_stats}->{$direction}));
        for my $tr (qw(trim raw)) {
            my $trdata = $tr . "data";
            $rec->{data_stats}->{$direction}->{$trdata} = {} unless (defined ($rec->{data_stats}->{$direction}->{$trdata}));
            my $data_file = $rec->{$direction}->{$trdata};
            if ($data_file and (-e $data_file)) {
                if (defined ($rec->{data_stats}->{$direction}->{$trdata}->{read_length}) and 
                    defined ($rec->{data_stats}->{$direction}->{$trdata}->{num_reads})) {
                    if ($options->{verbose}) {
                        print "Already calculated read info for file " . $rec->{$direction}->{$trdata} . "\n";
                    }
                } else {
                    my ($read_length, $num_reads) = get_read_info($data_file);
                    $rec->{data_stats}->{$direction}->{$trdata}->{read_length} = $read_length;
                    $rec->{data_stats}->{$direction}->{$trdata}->{num_reads} = $num_reads;
                }
            } else {
                if ($options->{verbose}) {
                    print "Data file " . $data_file . " not found for sample " . $rec->{sample} . " direction $direction, type $tr\n";
                }
            }
        }
    }
}

sub calc_data_stats
{
    unless ($options->{no_raw_stats}) {
        my $rec;
        if ($options->{sample} and defined $records->{$options->{sample}}) {
            my $rec = $records->{$options->{sample}};
            calc_stats($rec);
        } else {
            for my $sample (keys %{$records}) {
                my $rec = $records->{$sample};
                calc_stats($rec);
            }
        }
    }            
}  

sub get_reads
{
    my $kref = shift;
    my $ref = $records;
    for my $key (@$kref) {
        if (defined $ref->{$key}) {
            $ref = $ref->{$key};
        } else {
            return ('','');
        }
    }
    return ($ref->{read_length}, $ref->{num_reads});
}

sub write_table_stats
{
    my $outtable = ($options->{output_read_table} ? $options->{output_read_table} : '');
    if ($outtable) {
        open (FOUTTAB, '>', $outtable) or die "Error: couldn't open output file $outtable\n";
        print FOUTTAB join("\t", qw(Sample Trim_or_Raw Direction Read_Length Num_Reads)) . "\n";
        for my $sample (keys %$records) {
            for my $direction (qw(R1 R2)) {
                for my $tr (qw(trim raw)) {
                    my ($read_length, $num_reads) = get_reads([$sample, "data_stats", $direction, $tr]);
                    print FOUTTAB join("\t", ($sample, $tr . "data", $direction, $read_length, $num_reads)) . "\n";
                }
            }
        }
    }
}

gather_opts;
$records = LoadFile($options->{yaml_in});
get_illumina_stats;
read_table_stats;
calc_data_stats;
write_table_stats;
DumpFile($options->{yaml_out}, $records);


