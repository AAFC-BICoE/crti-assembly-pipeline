#!/bin/bash

# Note that there are many command-line options available
# for each of the scripts below. However, defaults are set
# at the top of each file and generally don't need to be
# changed. Exceptions to this rule should be commented below.

# Setup:
# cd /isilon/biodiversity/projects/CRTI/specimen/processing
# svn co http://biodiversity/svn/source/AssemblyPipeline
# cd AssemblyPipeline
# mkdir output_files qsub_files yaml_files 

# Note: nearly all of the scripts have a 'set_default_opts' function
# at the top of the script. An easy thing to do before running any
# of the scripts below is head or less on the script file to
# quickly check that the defaults make sense for each run.

# Get the latest svn version of the Illumina_sample_summary.xlsx
# in tab-delimited text format in the input_files/ directory.
./IlluminaSamples2Tab.pl

# Match the Illumina sample summary metadata with rawdata and
# use this to build the processing directory structure.
# Note: set the --specimen-dir <dirname> flag to specify where
# the output dir hierarchy should be built. Default is ../../processing_test
./DirSetup.pl

# Submit FastQC jobs on all raw files where this hasn't yet been done.
./FastQC --raw --submit

# Add any new trim values to input_files/ManualTrimParams.tab, then run:
./FastxTrim.pl
# to submit the jobs.
# Note: for samples that are missing trim values - the output file
# output_files/untrimmed_reports.html can be opened and contains links
# to the FastQC html reports so trim values can be determined for these
# more easily. Then add the new trim values back into ManualTrimParams.tab and re-run.

# Submit FastQC jobs on all trimmed files
./FastQC --trim --submit

# Calculate all raw/trimmed read length and num reads.
# This can be a very slow step. Use the --no_raw_reads option
# if you only wish to incorporate already-known values into the 
# output yaml file.
./GetReadInfo.pl

# Add any new estimated genome sizes to the file
# input_data/GenomeLengthEst.tab
# Note that any genomic samples without estimated genome sizes and/or
# number of reads / read lengths defined cannot be assembled.
# Run the following to pull this information into the YAML data for 
# availability in the next steps of the pipeline.
./GenomeLengths.pl

# Pull best kmers found from velvet advisor
# http://dna.med.monash.edu.au/~torsten/velvet_advisor/
./VelvetAdvisor.pl

# Get predicted best k-mer sizes for all raw/trimmed data using velvetk.pl
./VelvetKRun.pl
# The output output_files/VelvetKBestOut.pl now contains all of these values.
# For the sake of future runs, do:
# svn commit
# cp output_files/VelvetKBestOut.tab input_files/VelvetKBest.tab
# svn commit

# Create all velveth/velvetg commands in the output yaml file, using 13 kmers centered 
# at the velvetk best for each.
./VelvetCommands.pl --use_velvetk

# Run velveth on all trimmed and raw data, 
./VelvetHQsub_host.pl --submit
# This may submit a large number of qsub jobs.
# Note: jobs are submitted to target hosts such that the max number running on any node
# is less than the expected memory usage of the jobs. Excess jobs hold using -hold_jid.
# Note: To keep track of jobids, should afterwards run:
# svn commit
# cp output_files/VHQsubJobIDsOut.tab input_files/VHQsubJobIDs.tab
# svn commit
# Note: jobs are only submitted if valid velveth outputs are not found in each k-mers output dir.
# Note: use --submit_max <value> to set the max number of qsub jobs to submit.

# Do the same as above for velvetg:
./VelvetGQsub.pl --submit
# Note: to run just for high-priority samples, use --sample_list SampleID1,SampleID2,...
# (no spaces, just commas to separate the list).

# Gather stats about the velvet assemblies
./VelvetStats.pl
# Stats output table: output_files/VelvetStats.tab

# Gather information about velvetg runtimes and memory usage (and possibly in future, failed jobs)
./VelvetGUsage.pl





