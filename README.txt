

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
./DirSetup.pl -s input_data/Illumina_sample_summary.tab -a input_data/SpeciesAbbreviations.tab --specimen_dir ../../processing_test/ --yaml_out yaml_files/Y01_dirsetup.yml --all_samples

Mandatory:
    -s -a -yaml_out

Optional:
    --species <species> 
    --sample <sample> 
    --sample_file <sample file>
    --all_samples
    --testing
    --verbose
    --record_file <output records table name>
    --specimen_dir <path to specimen/ dir>

Error check: 

Run
cat yaml_files/01_dirsetup.yml | grep 'rawdata:' | wc -l

You should get back twice the number of samples you expect to have in your rawdata folder (sum of # of R1 and R2 files).


3. 

Run fastqc on the raw data.

./FastQC.pl -i yaml_files/01_dirsetup.yml -o yaml_files/02_rawqc.yml --raw --qsub_script qsub_script.sh --qsub_batch_file qsub_files/01_raw_fastqc.sh

cat qsub_files/01_raw_fastqc.sh | bash

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

Run post-trimming fastqc

Simply use the FastQC.pl script again, this time using the --trim option instead of --raw.

./FastQC.pl -i yaml_files/03_trim.yml -o yaml_files/04_trimqc.yml --trim --qsub_script qsub_script.sh --qsub_batch_file qsub_files/03_trim_fastqc.sh

6.

Go through raw/trimmed files to determine number of reads, read lengths.

Determining the number of reads is a very slow step. And is compounded by the fact that we do it for R1 and R2 for both raw and trimmed data for each sample.

Ideally, the raw/trimmed files should only actually be opened and parsed once; the outputs from the initial parsing should then be output to a table that can then be used as an input thereafter.

Also load 'number of reads' column from Illumina_sample_summary.tab file.

(Possibly compare the values obtained from each source in future?)



Take an input 'estimated genome lengths' table and add this information to the records.

7.

(Step 1 Genome Assembly part)










