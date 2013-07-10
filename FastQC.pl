#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use 5.010;
use YAML::XS qw(LoadFile DumpFile);
use File::Path;
use File::Basename;

my $fastqc_bin = "/opt/bio/FastQC/fastqc";
my $qsub_bin = "/opt/gridengine/bin/lx26-amd64/qsub";

my $options = {};
my $records = {};
my @qsub_cmd_list;


sub check_opts
{
    unless ($options->{yaml_in} and $options->{yaml_out} and ($options->{trim} or $options->{raw})) {
        die "Usage: $0 -i <yaml input file> -o <yaml output file> [--trim and/or --raw]
            Optional args:
            	--raw
            	--trim
            	--sample <sample id>
                --verbose
                --testing
                --qsub_opts <qsub options>
                --qsub_script <qsub script name>
                --qsub_batch_file
                --submit
                ";
    }
}

sub gather_opts
{
	GetOptions($options,
			'yaml_in|i=s',
			'yaml_out|o=s',
			'raw',
			'trim',
			'sample=s',
			'verbose|v',
			'testing',
			'qsub_opts=s',
			'qsub_script=s',
			'qsub_batch_file=s',
			'submit',
			);
	check_opts;
}

sub check_record_fields
{
    my $rec = shift;
    my $rval = shift;
    my $tr = shift;
    my $trdata = shift;
    unless (defined $rec->{sample}) {
        print "No sample for input record! Skipping.\n";
        return 0;
    }
    unless (defined $rec->{sample_dir}) {
        print "No sample_dir specified for sample " . $rec->{sample} . " - skipping\n";
        return 0;
    }
    unless (defined $rec->{$rval}) {
        print "Skipping sample " . $rec->{sample} . " - no data found for direction $rval.\n";
        return 0;
    }
    unless (defined $rec->{$rval}->{$trdata} and (-e $rec->{$rval}->{$trdata})) {
        print "Skipping sample " . $rec->{sample} . " - required $tr data file is missing or nonexistent.\n";
        return 0;
    }
    return 1;
}

sub get_fastqc_subdir
{
    my $rec = shift;
    my $rval = shift;
    my $tr = shift;
    my $tr_filename = shift;
    
    my $fastqc_subdir = $rec->{sample_dir} . "/data/fastqc_out/" . $tr_filename . "_fastqc";
    my $key = $tr . "_subdir";
    $rec->{fastqc}->{$rval}->{$key} = $fastqc_subdir;
    return $fastqc_subdir;
}

sub get_fastqc_cmd
{
    my $rec = shift;
    my $rval = shift;
    my $tr = shift;
    my $trdata = shift;
    my $fastqc_out = shift;
    my $fastqc_cmd = $fastqc_bin . " -o " . $fastqc_out . " " . $rec->{$rval}->{$trdata};
    #my $fastqc_cmd_key = $trdata . "_fastqc_cmd";
    #$rec->{$rval}->{$fastqc_cmd_key} = $fastqc_cmd;
    my $key = $tr . "_cmd";
    $rec->{fastqc}->{$rval}->{$key} = $fastqc_cmd;
    return $fastqc_cmd;
}

sub get_qsub_cmd
{
    my $rec = shift;
    my $rval = shift;
    my $tr = shift;
    my $fastqc_cmd = shift;
    $options->{qsub_opts} = ($options->{qsub_opts} ? $options->{qsub_opts} : '');
    $options->{qsub_opts} .= " -N " . $tr . "_FastQC ";
    $options->{qsub_script} = ($options->{qsub_script} ? $options->{qsub_script} : '');
    my $qsub_cmd = $qsub_bin . " " . $options->{qsub_opts} . " " . $options->{qsub_script} . " '" . $fastqc_cmd . "'";
    my $key = $tr . "_qsub";
    $rec->{fastqc}->{$rval}->{$key} = $qsub_cmd;
    return $qsub_cmd;
}

sub submit_qsub
{
    my $rec = shift;
    my $filename = shift;
    my $fastqc_subdir = shift;
    my $qsub_cmd = shift;
    #print join("\n", ("submit qsub:", $filename, $fastqc_subdir, $qsub_cmd)) . "\n";
    if (-e $fastqc_subdir) {
        if ($options->{verbose}) {
            print "Won't run fastqc on file $filename - output folder already exists:\n";
            print $fastqc_subdir . "\n";
        }
    } else {
        if ($options->{verbose}) {
            print "Submitting qsub command:\n${qsub_cmd}\n";
            my $fqcout = ($rec->{fastqc}->{fastqc_out} ? $rec->{fastqc}->{fastqc_out} : '');
            print "didn't find the out dir for this one - contents of fastqc_out:\n";
            print `/bin/ls $fqcout` if $fqcout;
        }
        if ($options->{qsub_batch_file}) {
            push (@qsub_cmd_list, $qsub_cmd);
        }
        if ($options->{submit}) {
            system($qsub_cmd);
        }
    }
}

sub get_sample_record
{
	my $temp = LoadFile($options->{yaml_in});
	my $sample = $options->{sample};
	if (defined $temp->{$sample}) {
		$records->{$sample} = $temp->{$sample};
	} else {
		die "Error: could not find a record for specified sample $sample in file " . $options->{yaml_in} . "\n";
	}
}

sub get_reports
{
    my $rec = shift;
    my $rval = shift;
    my $tr = shift;
    my $fastqc_subdir = shift;
    my $html_path = $fastqc_subdir . "/fastqc_report.html";
    my $text_path = $fastqc_subdir . "/fastqc_data.txt";
    my $html_key = $tr . "_report_html";
    my $text_key = $tr . "_data_txt";
    #if (-e $fastqc_subdir) {
    #    if (-e $html_path and -e $txt_path) {
    #        $rec->{fastqc}->{$rval}->{$key} = $html_path;
    #        $rec->{fastqc}->{$rval}->{$key} = $text_path;
    #    } else {
    #        die "Error: have fastqc output dir but it contains no html and/or text report\n";
    #    }
    #} else {
        $rec->{fastqc}->{$rval}->{$html_key} = $html_path;
        $rec->{fastqc}->{$rval}->{$text_key} = $text_path;
    #}
}

gather_opts;

if ($options->{sample}) {
	get_sample_record;
} else {
	$records = LoadFile($options->{yaml_in});
}

foreach my $sample (keys %{$records})
{  
    my $rec = $records->{$sample};
    foreach my $tr (qw(trim raw)) {
        if ($options->{$tr}) {
        	my $trdata = $tr . "data";
        	if ($tr =~ /raw/) {
        		$trdata = $trdata . "_symlink";
        	}
            $rec->{fastqc} = {};
            for my $rval (qw(R1 R2)) {
                if (check_record_fields($rec, $rval, $tr, $trdata)) {
                    my $fastqc_out = $rec->{sample_dir} . "/data/fastqc_out";
                    $rec->{fastqc}->{fastqc_out} = $fastqc_out;
                    mkpath $fastqc_out;
                    $rec->{fastqc}->{$rval} = {};
                    my $tr_filename = basename($rec->{$rval}->{$trdata});
                    my $fastqc_subdir = get_fastqc_subdir($rec, $rval, $tr, $tr_filename);
                    get_reports($rec, $rval, $tr, $fastqc_subdir);
                    my $fastqc_cmd = get_fastqc_cmd($rec, $rval, $tr, $trdata, $fastqc_out);
                    my $qsub_cmd = get_qsub_cmd($rec, $rval, $tr, $fastqc_cmd);
                    submit_qsub($rec, $rec->{$rval}->{$trdata}, $fastqc_subdir, $qsub_cmd);
                } else {
                    print "Failed record check for sample $sample direction $rval.\n";
                }
            }
        }
    }
}
        
if ($options->{qsub_batch_file}) {
    open (FQB, '>', $options->{qsub_batch_file}) or die "Error: could not open file " . $options->{qsub_batch_file} . "\n";
    print FQB join("\n", @qsub_cmd_list) . "\n";
    close (FQB);
}

if ($options->{yaml_out}) {
    DumpFile($options->{yaml_out}, $records);
}
    
    
    
    
    
    
       