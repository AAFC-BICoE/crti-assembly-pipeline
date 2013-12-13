#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Path;
use File::Basename;
use 5.010;
use YAML::XS qw(DumpFile);
use Assembly::Utils  qw(set_check_record get_check_record);
use Pod::Usage;

my $options = {};
my $records = {};
my %species_abbrs = ();


my @col_headers = qw (Shipping_Date Results_Received_NRC_PBI_LIMS_Request Project 
        Sequencing_Specifications Plate_Name Run_Name Lane Organism Genotype Growth_Media
        Growth_Condition_Notes Biomaterial Biomaterial_Type Sample_Name MID_Tag Barcode
        Total_Numreads_in_Lane NumReads Percent_of_Reads_in_Lane Download); 

=head1 NAME

DirSetup.pl

=head1 SYNOPSIS

DirSetup.pl --seq_sample_file input_data/Illumina_sample_summary.tab

=head1 OPTIONS

=over 8

=item B<--seq_sample_file>

Illumina summary table filename (or selected lines thereof) in tab-delimited text format.

=item B<--species_abbr_file>

File specifying two- or three-letter abbreviations for species names. Each abbreviation should be unique.

=item B<--specimen_dir>

The base folder for processing. All specimen folders and raw files will be created and linked at this location, and all further processing using the output metadata will occur within subfolders of this folder.

=item B<--yaml_out>

The output YAML file, containing metadata about file locations and steps performed. To be used as input in subsequent stages of pipeline processing.

=item B<--g_rosto_file>

Deprecated. Used to allow for nonstandard G rostochiensis file processing by pipeline. Deprecated as we now have our own (standard PBI) raw reads for G rostochiensis.

=item B<--all_samples>

Deprecated. Process all samples. Previous versions supported running specific samples only. Such can be performed now by selecting specific lines from the input sample sheet.

=item B<--verbose>

Print more detailed information on steps being performed.

=item B<--testing>

Step through the script but do not actually create any files or directories.

=item B<--defaults>

Print out all default settings for input option values.

=item B<--help>

Print out the help message.

=item B<--man>

Full documentation.

=back

=head1 DESCRIPTION

Initialize the folder structure in preparation for running the genome assembly pipeline.

=cut

sub set_default_opts
{
    my %defaults = qw(
        seq_sample_file input_data/Illumina_sample_summary.tab
        species_abbr_file input_data/SpeciesAbbreviations.tab
        specimen_dir ../../processing_test2
        yaml_out yaml_files/01_dirsetup.yml
        );
        # g_rosto_file input_data/G_rosto_sample_summary.tab
    for my $key (keys %defaults) {
        $options->{$key} = $defaults{$key} unless $options->{$key};
    }
    Assembly::Utils::print_defaults ($0, $options);
    #if ($options->{defaults}) {
    #    print "Defaults:\n";
    #    print Dumper($options);
    #    exit;
    #}
}

sub check_options
{
	if ($options->{help}) {
	    pod2usage(1);
	}
	if ($options->{man}) {
	    pod2usage(-exitval => 0, -verbose => 2);
	}
	unless ($options->{seq_sample_file} and $options->{species_abbr_file} ) {
		pod2usage(1);
	}
}

sub gather_options
{
	GetOptions($options,
		'seq_sample_file|s=s',
		'species_abbr_file|a=s',
		'all_samples',
		'testing|t',
		'verbose|v',
		'specimen_dir|p=s',
		'yaml_out|y=s',
		'g_rosto_file=s',
		'defaults|d',
		'help|h',
		'man|m',
		) or pod2usage(1);
	set_default_opts;
	check_options;
}		

sub parse_abbrevs
{
	open (FABBR, '<', $options->{species_abbr_file}) or die "Error: couldn't open file " . $options->{species_abbr_file} . "\n";
	while (my $line = <FABBR>) {
		chomp $line;
		my @fields = split(/\t/, $line);
		if (scalar @fields == 2) {
			my ($abbr, $species) = @fields;
			my $species_key = Assembly::Utils::format_species_key($species);
			$species_abbrs{$species_key} = $abbr;
		} else {
			die "Got wrong number of fields on line $. in abbrevs file " . $options->{species_abbr_file} . "\nline: $line\n";
		}
	}
	close (FABBR);
}

sub clean_field
{
	my $infield = shift;
	#strip any leading and trailing whitespace in field.
	$infield =~ s/^\s+|\s+$//g if $infield; 
	return $infield;
}

# Gather fields of interest from a line of the input .tab file.
sub parse_record
{
	my $line = shift;
	chomp $line;
	my @fields = split (/\t/, $line);
	my $rec = {};
	for (my $i=0; $i < scalar @col_headers; $i++) {
	    my $colname = $col_headers[$i];
	    my $value = ($fields[$i] ? $fields[$i] : ''); # i+1 because of blank first col in summary table.
	    $value = clean_field ($value);
	    Assembly::Utils::set_check_record($rec, ["sequencing_metadata"], $colname, $value);
	}
	$rec->{sample} = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Sample_Name"]);
	$rec->{bio_type} = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Biomaterial_Type"]);
	$rec->{plate} = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Plate_Name"]);
	$rec->{species} = Assembly::Utils::get_check_record($rec, ["sequencing_metadata", "Organism"]);
	$rec->{species} =~ s/\s\(Erwinia\)//g; # Not the best solution, but there it is.
	$rec->{species} =~ s/\s\-\sresequencing\s*$//; # as above
	return $rec;
}	

sub parse_line
{
    my $line = shift;
    chomp $line;
    my $rec = parse_record($line);
    if (defined $rec->{sample} and $rec->{sample} =~ /\S/) {
        $records->{$rec->{sample}} = $rec;
    }
}

# Go through input file and build records array matching the latter part
# (after last '_' char) of species name.
sub build_records
{
	open (FIN, '<', $options->{seq_sample_file}) or die "Error: could not open file " . $options->{seq_sample_file} . "\n";
	<FIN>; # Skip the parsing the first line (col headers).
	while (my $line = <FIN>) {
        parse_line($line);
	}
	close FIN;  
}

sub build_rosto_records
{
    open (FIN, '<', $options->{g_rosto_file}) or die "Error: couldn't open file " . $options->{g_rosto_records} . "\n";
    <FIN>; # Skip the parsing the first line (col headers).
	while (my $line = <FIN>) {
        #parse_line($line);
        chomp $line;
        my @fields = split(/\t/, $line);
        if (scalar @fields == 5) {
            my $rec = {};
            $rec->{plate} = $fields[0];
            $rec->{species} = $fields[1];
            $rec->{bio_type} = $fields[2];
            $rec->{sample} = $fields[3];
            $rec->{rawdata_dir} = $fields[4];
            $records->{$rec->{sample}} = $rec;
        }
	}
	close FIN;
}

sub species_to_dirname
{
	my $species_long = shift;
	if ($species_long) {
        my $dirname = '';
        if ($species_long =~ /^\s*([A-Z])\S+\s*([a-z]+)\s*/) {
            $dirname = $options->{specimen_dir} . "/" . $1 . "_" . $2;
        } elsif ($species_long =~ /([A-Z]_[a-z]+)/) {
            $dirname = $options->{specimen_dir} . "/" . $1;
        }
        return $dirname;
    }
    return '';
}

# Check whether driectory exists or just testing.
# Report status if verbose is set.
sub check_create_dir
{
    my $dir = shift;
    print "Attempting to create dir: $dir\n" if $options->{verbose};
    if (-e $dir) {
        if ($options->{verbose}) {
            print "Directory exists.\n";
        }
    } else {
        unless ($options->{testing}) {
            mkpath($dir);
        }
        if ($options->{verbose}) {
            print "Created directory\n";
        }
    }
}

sub find_rawdata_files
{
	my $rec = shift;
	my $plate = $rec->{plate};
	my $sample = $rec->{sample};
	my $species = $rec->{species};
	# Find any data files that go with the sample.
	my $rawdata_dir = "/isilon/biodiversity/data/raw/illumina/PBI/Project_" . $plate;
	if ($rec->{plate} =~ /G_ROSTO/) {
	    $rawdata_dir = $rec->{rawdata_dir};
	}
	print "Looking for directory " . $rawdata_dir . "\n" if $options->{verbose};
	my @rawfiles = ();
	if (-d $rawdata_dir) {
		opendir DIR, $rawdata_dir or die "Cannot open raw data directory: ${rawdata_dir}\n";
		my @files = readdir DIR;
		closedir DIR;
		@rawfiles = grep (/$sample/, @files);
		@rawfiles = map { $rawdata_dir . "/" . $_ } @rawfiles;
		print "Found rawdata files:\n" if $options->{verbose};
		print join ("\n", @rawfiles) . "\n" if $options->{verbose};
	} else {
		print "Warning: Raw data - no such directory: ${rawdata_dir}\n";
	}
	return \@rawfiles;
}

sub add_rawdata_record
{
	my $rec = shift;
	my $rawfile = shift;
	my $destfile = shift;
	if ($rawfile =~ /(R1|R2)\.f(ast)?q/) {
		my $direction = $1;
		$rec->{$direction} = {};
		$rec->{$direction}->{rawdata} = $rawfile;
		$rec->{$direction}->{rawdata_symlink} = $destfile;
	} else {
		print "Warning: could not parse direction R1 or R2 from raw filename: $rawfile\n";
	}
}

sub link_rawdata
{
	my $rec = shift;
	my $rawref = shift;
	my $datadir = shift;
	my @rawdata = @{$rawref};
	foreach my $rawfile (@rawdata) {
		my $destfile = $datadir . "/" . basename($rawfile);
		
		add_rawdata_record($rec, $rawfile, $destfile);
		if (-e $destfile) {
		    if ($options->{verbose}) {
			    print "Symlink: destination $destfile already exists\n";
		    }
		} elsif (-e $rawfile) {
            print "Symlinking $rawfile to $destfile\n" if $options->{verbose};
		    unless ($options->{testing}) {
			    symlink($rawfile, $destfile) or die "Error: could not create symlink.\n";
		    }
        }
	}
}

# Turn 'biomaterial_type' value into RNA or DNA.
sub get_type
{
	my $intype = shift;
	if ($intype =~ /(DNA|RNA)/) {
		$intype = $1;
	} elsif ($intype =~ /genomic/i) {
		$intype = "DNA";
	} else {
		print "Error: couldn't determine biomaterial type from field $intype.\n";
	}
	return $intype;
}

sub create_species_dirs
{
	my $rec = shift;
	my $specdir = $rec->{species_dir};
	check_create_dir($specdir);
	check_create_dir("$specdir/reference");
	check_create_dir("$specdir/releases");
	foreach my $subdir (qw(DNA RNA)) {
		check_create_dir("$specdir/$subdir");
		foreach my $ssd (qw(data assemblies annotations)) {
		    check_create_dir("$specdir/$subdir/$ssd");
		}
	}
}

sub create_sample_dirs
{
	my $rec = shift;
	my $specdir = $rec->{species_dir};
	my $seqtype = $rec->{bio_type};
	my $abbr = $rec->{species_abbr};
	
	my $sampledir = "$specdir/$seqtype/data/" . $abbr . "_" . $rec->{sample};
	$rec->{sample_dir} = $sampledir;
	check_create_dir($sampledir);
	return $sampledir;
}

sub create_directory_setup
{
	my $rec = shift;
	$rec->{species_dir} = species_to_dirname($rec->{species});
	$rec->{bio_type} = get_type($rec->{bio_type});
	my $species_key = Assembly::Utils::format_species_key($rec->{species});
	if (defined $species_abbrs{$species_key}) {
		$rec->{species_abbr} = $species_abbrs{$species_key};
		my $rawfiles_ref = find_rawdata_files($rec);
		my @rawfiles = @{$rawfiles_ref};
		if (scalar @rawfiles == 2) {
			create_species_dirs($rec);
			my $sampledir = create_sample_dirs($rec);
			link_rawdata($rec, $rawfiles_ref, $sampledir); # Also adds rawdata info to record.
		} elsif (scalar @rawfiles < 2) {
			print "Couldn't find raw data for species " . $rec->{species} . ", sample " . $rec->{sample} . ".\n";
		} else {
			#if (scalar @rawfiles > 2) {
			die "Error: too many raw data files with same sample ID! Aborting.\n";
		}
	} elsif ($options->{verbose}) {
		print "Could not determine species abbreviation - skipped creating dirs for species " . $rec->{species} . "\n";
	}
}

sub setup_all_dirs
{
	foreach my $sample (keys %$records) {
	    print "Got sample $sample\n";
	    if ($sample =~ /\S/) {
	        print "Working on sample $sample\n";
		    my $rec = $records->{$sample};
		    create_directory_setup($rec);
	    }
	}
}

gather_options;
parse_abbrevs;
build_records;
if ($options->{g_rosto_file}) {
    build_rosto_records;
}
setup_all_dirs;
DumpFile($options->{yaml_out}, $records);
