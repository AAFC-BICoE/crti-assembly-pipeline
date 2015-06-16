#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;

#my $bowtie2_build_path = "/opt/bio/bowtie2/bowtie2-build";
#my $bowtie2_path = "/opt/bio/bowtie2/bowtie2";
#my $samtools_path = "/opt/bio/samtools/samtools";
my $bowtie2_build_path = "bowtie2-build";
my $bowtie2_path = "bowtie2";
my $samtools_path = "samtools";

my $options = {};
GetOptions($options,
        'keep_unmapped_reads',
        'keep_mapped_reads',
        'genome|g=s',
        'reads_in_R1=s',
        'reads_in_R2=s',
        'reads_out_R1=s',
        'reads_out_R2=s',
        'sam_out=s',
        'delete_sam',);
unless ($options->{keep_mapped_reads} xor $options->{keep_unmapped_reads}) {
    die "Error: exactly one of --keep_mapped_reads or --keep_unmapped_reads must be set\n";
}
for my $opt (qw(genome reads_in_R1 reads_in_R2 reads_out_R1 reads_out_R2)) {
    unless ($options->{$opt}) {
        die "Error: option $opt must be set\n";
    }
}

sub run_bowtie2_build
{
    my $index_exists = 1;
    for my $ext (qw(1 2 3 4 rev.1 rev.2)) {
        my $idx_file = getcwd . "/" . $options->{genome} . "." . $ext . ".bt2";
        unless (-e $idx_file) {
            $index_exists = 0;
        }
    }
    unless ($index_exists) {
        my $output_prefix = getcwd . "/" . basename ($options->{genome});
        my $bowtie2_build_cmd = $bowtie2_build_path . " " . $options->{genome} . 
                " " . $output_prefix;
        print $bowtie2_build_cmd . "\n";
        my $out = `$bowtie2_build_cmd`;
        print $out;
    }
}
        

sub run_bowtie2
{
    my $samfile = $options->{sam_out};
    my $ri1 = $options->{reads_in_R1};
    my $ri2 = $options->{reads_in_R2};
    my $bowtie2_cmd = $bowtie2_path . " -x " . $options->{genome} . " -q -1 " . $ri1 . " -2 " .
        $ri2 . " -S " . $samfile;
    print $bowtie2_cmd . "\n";
    my $out = `$bowtie2_cmd`;
    print $out;
    return $samfile;
}

# See http://seqanswers.com/forums/showthread.php?t=16743
# for more info on the following two funcs
# Position 4 defines unmapped reads. Only
# difference in following two funcitons is -f4 vs -F4

sub get_unmapped_reads
{
    my $sam_infile = shift;
    my $sam_outfile = $sam_infile;
    $sam_outfile =~ s/\.sam/_unmapped.sam/;
    # samtools view -S -f4 Se2012_to_Pf.sam
    my $sam_cmd = $samtools_path . " view -S -f4 " . $sam_infile . " > " . $sam_outfile . "\n";
    print $sam_cmd;
    my $out = `$sam_cmd`;
    return $sam_outfile;
}

sub get_mapped_reads
{
    my $sam_infile = shift;
    my $sam_outfile = $sam_infile;
    $sam_outfile =~ s/\.sam/_mapped.sam/;    
    my $sam_cmd = $samtools_path . " view -S -F4 " . $sam_infile . " > " . $sam_outfile . "\n";
    print $sam_cmd;
    my $out = `$sam_cmd`;
    return $sam_outfile;
}

sub sam2fastq
{
    my $ump_samfile = shift;
    
    my $ro1 = $options->{reads_out_R1};
    my $ro2 = $options->{reads_out_R2};
    open (FSAM, '<', $ump_samfile) or die "Error opening file $ump_samfile.\n";
    open (FFQR1, '>', $ro1) or die "Error opening file $ro1.\n";
    open (FFQR2, '>', $ro2) or die "Error opening file $ro2.\n";
    while (my $line = <FSAM>) {
        chomp $line;
        if ($line !~ /^@/) {
            my @fields = split (/\s+/, $line);
            if (scalar @fields >= 13) {
                if (scalar @fields % 2 == 1) {
                    print FFQR1 "@" . $fields[0] . "\n";
                    print FFQR1 $fields[9] . "\n";
                    print FFQR1 "+\n";
                    print FFQR1 $fields[10] . "\n";
                } elsif (scalar @fields %2 == 0) {
                    print FFQR2 "@" . $fields[0] . "\n";
                    print FFQR2 $fields[9] . "\n";
                    print FFQR2 "+\n";
                    print FFQR2 $fields[10] . "\n";
                }
            }
        }
    }
    close (FSAM);
    close (FFQR1);
    close (FFQR2);
}

sub run_alignment
{
    run_bowtie2_build;
    my $bwt_samfile = run_bowtie2;
    
    my $ump_samfile = '';
    if ($options->{keep_unmapped_reads}) {
        print "Getting unmapped reads\n";
        $ump_samfile = get_unmapped_reads($bwt_samfile);
    } elsif ($options->{keep_mapped_reads}) {
        print "Getting mapped reads\n";
        $ump_samfile = get_mapped_reads($bwt_samfile);
    }
    if ($ump_samfile and -e $ump_samfile) {
        sam2fastq ($ump_samfile);
    }
    # cleanup
    if ($options->{delete_sam} and not $options->{testing}) {
        if (-e $bwt_samfile) {
             system("rm $bwt_samfile"); 
        }
        if (-e $ump_samfile) {
            system("rm $ump_samfile");
        }
    }
}

run_alignment;
        
