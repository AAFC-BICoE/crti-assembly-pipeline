

qsub_script.sh
takes one argument, which is the command to be run.
e.g.
qsub -N FastQC qsub_script.sh 'fastqc -o ...'


mkdir yaml_files
mkdir output_data
mkdir qsub_files

1.
./IlluminaSamples2Tab.pl

Downloads from 
http://biodiversity/svn/sequencing/Illumina_sample_summary.xlsx
Saves as 
./input_data/Illumina_sample_summary.xlsx
and
./input_data/Illimina_sample_summary.tab

by default. 

For custom download, specify -u <url> -x <output xlsx download path> -t <output tab file path>




2.

Assume we're in specimen/processing/AssemblyPipeline/ dir and running a test run in dir specimen/processing_test/

Example:
./DirSetup.pl

Mandatory:
    -s -a -yaml_out

All:
    --seq_sample_file <file> default: input_data/Illumina_sample_summary.tab
    --species_abbr_file <file> default: input_data/SpeciesAbbreviations.tab
    --specimen_dir <dirname> default: ../../processing_test/
    --all_samples default: set
    --testing
    --verbose
    --specimen_dir <path to specimen/ dir>

Error check: 

Run
cat yaml_files/01_dirsetup.yml | grep 'rawdata:' | wc -l

You should get back twice the number of samples you expect to have in your rawdata folder (sum of # of R1 and R2 files).

NOTE: You can run the pipeline on a subset of your samples by pulling the relevant records out of the illumina summary table, e.g.:

cat input_data/illumina_summary.tab | grep lanthierii | grep 1440 > input_data/illumina_summary_lanthierii.tab
./DirSetup.pl -s input_data/illumina_summary_lanthierii.tab

3. 

Run fastqc on the raw data.

./FastQC.pl --raw 

Optional args:
    --raw
    --trim
    --sample <sample id>
    --verbose
    --testing
    --qsub_opts <qsub options>
    --qsub_script <qsub script name>
    --qsub_batch_file

4. 

./FastxTrim.pl
firefox output_files/untrimmed_reports.html
record trim params in input_data/ManualTrimParams.tab
-Q 33
./FastxTrim.pl

$ ./FastxTrim.pl --verbose
Warning: too few fields in file input_data/ManualTrimParams.tab line 296
Your job 103501 ("fastx_trim") has been submitted
Your job 103502 ("fastx_trim") has been submitted




Parse trim parameters and generate trim commands/qsub scripts.

Read in any available trim parameters in table form. Apply these to the input yaml file. Construct the fastx_trimmer commands and qsubs. Create a batch script of qsub jobs.

Build some links to the fastqc reports for those files where we don't have manual trim parameters.

./FastxTrim.pl -i yaml_files/02_rawqc.yml -o yaml_files/03_trim.yml --trim_params_table input_data/ManualTrimParams.tab --qsub_script qsub_script.sh --qsub_batch_file qsub_files/02_fastx_trim.sh --report_notrim output_files/01_qc_reports_notrim.html


Optional args:
    --trim_params_table <file name>
    --sample <sample id>
    --verbose
    --testing
    --qsub_opts <qsub options>
    --qsub_script <qsub script name>
    --qsub_batch_file
    --run_fastx

Note that the --report_notrim command is not mandatory but is used here. This file generates links to the raw fastqc report html file for any of our records where we did not find any trim parameters. After running FastxTrim.pl, load the file to quickly perform fastqc on those records and come up with trim values to apply to the input_data/ManualTrimParams.tab table and run again.

Once trim commands/qsub commands have been generated for all samples, run:
cat qsub_files/02_fastx_trim.sh | bash to submit all the jobs.

5.

./FastQC.pl --trim

Run post-trimming fastqc

Simply use the FastQC.pl script again, this time using the --trim option instead of --raw.

./FastQC.pl -i yaml_files/03_trim.yml -o yaml_files/04_trimqc.yml --trim --qsub_script qsub_script.sh --qsub_batch_file qsub_files/03_trim_fastqc.sh

6.

screen
source ~/perl5/perlbrew/etc/bashrc
./GetReadLengths.pl --verbose

Go through raw/trimmed files to determine number of reads, read lengths.

Determining the number of reads is a very slow step. And is compounded by the fact that we do it for R1 and R2 for both raw and trimmed data for each sample.

Ideally, the raw/trimmed files should only actually be opened and parsed once; the outputs from the initial parsing should then be output to a table that can then be used as an input thereafter.

Also load 'number of reads' column from Illumina_sample_summary.tab file.

(Possibly compare the values obtained from each source in future?)




Take an input 'estimated genome lengths' table and add this information to the records.

7.

(Step 1 Genome Assembly part)


./AssemblySetup.pl --verbose

./GenomeLengths.pl --verbose

./ReverseComplement.pl --verbose

./VelvetAdvisor.pl --verbose

./VelvetKRun.pl --verbose

./VelvetCommands.pl --verbose

./VelvetStatsFull.pl --verbose

Next - release
Look at Velvet Stats file and add the best one to the input_data/Releases.tab file.

ContigStats.pl

Had to flip the comments so the last few lines in COntigStats.pl to look like:

my $output_stats = all_contig_stats($records);
#my $output_stats = release_contig_stats($records);

Then ran

ContigStats.pl --verbose

then look at output_files/ContigStatsOut.tab
Pull the line for the contig file to release and copy it to
input_data/ContigStats.tab

Then switch the comment for the lines above and run ContigStats.pl again?

From what I can see/remember, the output ContigStats.tab is not as important as just making sure the 
ContigStats.pl has been run so that the min/median/max contig lengths have been calculated and
added to the yaml file (yaml_files/13_contig_stats.yml) so that the info can be pulled by the release
script.

ReleaseGenomes.pl

First edited the input_data/Release.tab
added a line for the new species/strain (Arcobacter lanthierii AF_1440) using the stats from
output_files/VelvetStatsFull.pl for the best kmer (91).

Then had to do a
qlogin

because the ReleaseGenomes.pl script needs to be able to run e.g. 'fastqc --version' to add
the version number info
(Bug: This should be done just *before* we actually run each program, e.g. in VelvetCommands.pl).

Now we see:

[cullisj@biocluster AssemblyPipeline]$ ls ../../processing_lanthierii/A_lanthierii/release/Al_AF1440_R01/
Al_AF1440_R01_contigs.fa                Al_AF1440_R01_DNA_PE_S0017B0_R2.fq.gz*      Al_AF1440_R01_DNA_PE_S0017B0_trim_R2.fq.gz
Al_AF1440_R01_DNA_PE_S0017B0_R1.fq.gz*  Al_AF1440_R01_DNA_PE_S0017B0_trim_R1.fq.gz  Al_AF1440_R01_metadata.yml

The next step is to change the contig node names from the default values to the release prefix followed by an incremental integer.
Use the ReleaseContigs.pl script for this.

[cullisj@biocluster AssemblyPipeline]$ ./ReleaseContigs.pl --contigs_in ../../processing_lanthierii/A_lanthierii/release/Al_AF1440_R01/Al_AF1440_R01_contigs.fa --contigs_out ../../processing_lanthierii/A_lanthierii/release/Al_AF1440_R01/Al_AF1440_R01_contigs2.fa
(Potential bug: note that the release prefix that is used is actualyl the part before
_contigs.fa in the **input contigs filenam only**! This should be changed.)

Now check that the node names are correct in the output file
less !$

Now remove the original and replace with the new one

[cullisj@biocluster AssemblyPipeline]$ mv ../../processing_lanthierii/A_lanthierii/release/Al_AF1440_R01/Al_AF1440_R01_contigs.fa ~/.Trash/
[cullisj@biocluster AssemblyPipeline]$ mv ../../processing_lanthierii/A_lanthierii/release/Al_AF1440_R01/Al_AF1440_R01_contigs2.fa ../../processing_lanthierii/A_lanthierii/release/Al_AF1440_R01/Al_AF1440_R01_contigs.fa


Now the wiki table outputs are at
output_files/wiki_combined_table.txt

Move the folder to the correct folder within specimen/dir:
[cullisj@biocluster AssemblyPipeline]$ cd ../..
[cullisj@biocluster specimen]$ mv processing_lanthierii/A_lanthierii/release/Al_AF1440_R01 A_lanthierii/release/

Now modify the link in the wiki table above to reflect the new folder location, and add the 
new line to the wiki.

Done!

Now to assemble the transcriptome for the released genome...


