

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

Build some links to the fastqc reports for those files where we don't have manual trim parameters.

