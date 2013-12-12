#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use YAML::XS qw(DumpFile LoadFile);
use Cwd;

my $options = {};
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";
my $fastx_trimmer_bin = "/opt/bio/fastx/bin/fastx_trimmer";
my $gunzip_bin = "/bin/gunzip";
my $fastx_gunzip_script = getcwd . "/FastxTrim.sh";
my @qsub_cmd_list = ();
my @html_lines = ();
my $records = {};

sub set_default_opts
{
    my %defaults = qw(
        yaml_in yaml_files/02_rawqc.yml
        yaml_out yaml_files/03_trim.yml
        trim_params_table input_data/ManualTrimParams.tab
        qsub_script qsub_script.sh
        qsub_batch_file qsub_files/02_fastx_trim.sh
        report_notrim output_files/untrimmed_reports.html
        );
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
    $options->{qsub_opts} = " -N fastx_trim " . $options->{qsub_opts};
}

sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out}) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file> 
            Optional args:
                --trim_params_table <file name>
                --report_notrim <output html filename>
            	--sample <sample id>
                --verbose
                --testing
                --qsub_opts <qsub options>
                --qsub_script <qsub script name>
                --qsub_batch_file
                ";
    }
}

sub gather_opts
{
    $options->{qsub_opts} = '';
    GetOptions($options,
            'yaml_in|i=s',
            'yaml_out|o=s',
            'trim_params_table|t=s',
            'report_notrim|r=s',
            'verbose|v',
            'testing',
            'qsub_opts=s',
            'qsub_script=s',
            'qsub_batch_file=s',
            );
    set_default_opts;
    check_opts;
}

sub get_subrecord
{
    my $rec = shift;
    my $kref = shift;
    my @keylist = @$kref;
    my $ref = $rec;
    for my $keyname (@keylist) {
        if (defined $ref->{$keyname}) {
            $ref = $ref->{$keyname};
        } else {
            return '';
        }
    }
    return $ref;
}

sub parse_trim_params
{
    my $ft = ($options->{trim_params_table} ? $options->{trim_params_table} : '');
    if (-e $ft) {
        open (FTRIM, '<', $ft) or die "Error: could not open trim params file $ft\n";
        <FTRIM>; # Skip header 1st line
        while (my $line = <FTRIM>) {
            chomp $line;
            my @fields = split(/\s+/, $line);
            if (scalar @fields >= 5) {
                my ($sample, $type, $direction, $Qval, $fval) = @fields[0..4];
                if ($sample =~ /[A-Z][a-z]_(S00[A-Z0-9]+)/) {
                    $sample = $1;
                }
                if (defined $records->{$sample}) {
                    my $rec = $records->{$sample};
                    unless (defined $rec->{fastx_trimmer}) { $rec->{fastx_trimmer} = {}; }
                    $rec->{fastx_trimmer}->{$direction} = {};
                    $rec->{fastx_trimmer}->{$direction}->{fval} = $fval;
                    $rec->{fastx_trimmer}->{$direction}->{Qval} = $Qval;
                }
            } else {
                print "Warning: too few fields in file $ft line $.\n";
            }
        }
    } else {
        die "Error: couldn't parse trim params file $ft\n";
    }
}

sub get_trim_cmd
{
    my $rec = shift;
    my $direction = shift;
    my $rawfile = $rec->{$direction}->{rawdata_symlink};
    my $outfile = $rec->{$direction}->{trimdata};
    my $Qval = $rec->{fastx_trimmer}->{$direction}->{Qval};
    my $fval = $rec->{fastx_trimmer}->{$direction}->{fval};
    my $cmd = '';
    if ($rawfile =~ /q\.gz\s*$/) {
        $cmd = join (" ", ($fastx_gunzip_script, $rawfile, $Qval, $fval, $outfile));
    } else {
        $cmd = $fastx_trimmer_bin . " -Q $Qval -f $fval -i $rawfile -z -o $outfile";
    }
    $rec->{fastx_trimmer}->{$direction}->{trim_cmd} = $cmd;
    return $cmd;
}

sub get_qsub_cmd
{
    my $rec = shift;
    my $direction = shift;
    my $trim_cmd = $rec->{fastx_trimmer}->{$direction}->{trim_cmd};
    my $qsub_cmd = $qsub_bin;
    if ($options->{qsub_opts}) {
        $qsub_cmd .= " " . $options->{qsub_opts};
    }
    if ($options->{qsub_script}) {
        $qsub_cmd .= " " . $options->{qsub_script};
    }
    $qsub_cmd .= " '" . $trim_cmd . "'";
    $rec->{fastx_trimmer}->{$direction}->{qsub_cmd} = $qsub_cmd;
    return $qsub_cmd;
}

# For each record, generate the trimmed output filename from params
# If output filename doesn't exist, generate the fastx_trimmer command
sub apply_trim_params
{
    foreach my $sample (keys %{$records}) {
        my $rec = $records->{$sample};
        if (defined $rec->{fastx_trimmer}) {
            foreach my $direction (qw(R1 R2)) {
                my $fval = get_subrecord($rec, ['fastx_trimmer', $direction, 'fval']);
                my $Qval = get_subrecord($rec, ['fastx_trimmer', $direction, 'Qval']);
                if ($fval and $Qval) {
                    #my $fval = $rec->{fastx_trimmer}->{$direction}->{fval};
                    #my $Qval = $rec->{fastx_trimmer}->{$direction}->{Qval};
                    my $rawfile = $rec->{$direction}->{rawdata_symlink};
                    my $trimfile = '';
                    if ($rawfile =~ /(.*)_(R[12].fq)/) {
                        my $trimval = $fval - 1; # The number of bases actually trimmed is f-1.
                        $trimfile = $1 . "_trim_" . $trimval . "-0-0_" . $2 . ".gz";
                        $rec->{$direction}->{trimdata} = $trimfile;
                    } else {
                        die "Error: raw file $rawfile cannot be used to get a trim file. Aborting";
                    }
                    my $trim_cmd = get_trim_cmd($rec, $direction);
                    my $qsub_cmd = get_qsub_cmd($rec, $direction);
                    if (-e $trimfile) {
                        if ($options->{verbose}) {
                            print "Trimmed file $trimfile already exists\n";
                        }
                    } else {
                        push (@qsub_cmd_list, $qsub_cmd);
                        unless ($options->{testing}) {
                            system($qsub_cmd);
                        }
                    }
                } else {
                    print "Error: found fastx_trimmer record for $sample $direction but fval and qval are undefined!\n";
                }
            }
        } elsif ($options->{report_notrim}) {
            for my $direction (qw(R1 R2)) {
                #my $html_report = $rec->{fastqc}->{$direction}->{raw_report_html};
                my $html_report = get_subrecord($rec, ['fastqc', $direction, 'raw_report_html']);
                $html_report = "../" . $html_report; # since the report will be put in the output_files/ subdir.
                if ($html_report) {
                    push (@html_lines, "<tr><td>$sample</td><td>$direction</td><td><a href='${html_report}'>report</a></td></tr>");
                } else {
                    print "Warning: no raw fastqc_report.html file specified for sample $sample, direction $direction\n";
                }
            }
        }
    }           
}

# Gather reports into output html file for
# any raw files we find that have no trimmed counterpart.
sub write_notrim_report
{
    my $outfile = $options->{report_notrim};
    open (FREP, '>', $outfile) or die "Error: could not open file $outfile\n";
    print FREP "<html><body><table>\n";
    print FREP join("\n", @html_lines) . "\n";
    print FREP "</table></body></html>";
    close (FREP);
}
                    
sub write_batch
{
    my $fb = $options->{qsub_batch_file};
    open (FBATCH, '>', $fb) or die "Error: couldn't open file $fb\n";
    print FBATCH join("\n", @qsub_cmd_list) . "\n";
    close (FBATCH);
}


gather_opts;
$records = LoadFile($options->{yaml_in});
if ($options->{trim_params_table}) {
    parse_trim_params;
}
apply_trim_params;
if ($options->{report_notrim}) {
    write_notrim_report;
}
if ($options->{qsub_batch_file}) {
    write_batch;
}
DumpFile($options->{yaml_out}, $records);














