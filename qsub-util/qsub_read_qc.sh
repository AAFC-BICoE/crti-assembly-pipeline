SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

# Note - need to be a bit careful about paired-end reads
# See 'Trimming paired-end reads' section of docs: https://cutadapt.readthedocs.org/en/latest/guide.html#basic-usage
# Forward and reverse adapter strings should be different.
run_cutadapt() {
    reads_file_in=$1
    reads_file_out=$2
    adapt_str=$3
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    cutadapt_cmd="/home/AAFC-AAC/cullisj/software/pyenv/versions/2.7.6/bin/cutadapt -a ${adapt_str} -o ${reads_file_out} ${reads_file_in}"
    jobname="cutadapt"
    jobid=`run_qsub 1 $qsub_holdid "$cutadapt_cmd" $jobname`
    echo $jobid
}

run_fastqc()
{
    reads=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    fastqc_cmd="/opt/bio/FastQC/fastqc $reads"
    jobname="fastqc"
    jobid=`run_qsub 1 $qsub_holdid "$fastqc_cmd" $jobname`
    echo $jobid
}

# Note: reads in must be unzipped reads
run_fastx_trimmer()
{
    reads_in=$1
    reads_out=$2
    first_base=$3
    last_base=$4
    qsub_holdid=1
    [ ! -z $5 ] && qsub_holdid=$5
    fastx_trimmer_cmd="/opt/bio/fastx/bin/fastx_trimmer -Q 33 -i $reads_in -o $reads_out -f $first_base -l $last_base"
    jobname="fastx_trimmer"
    jobid=`run_qsub 1 $qsub_holdid "$fastx_trimmer_cmd" $jobname`
    echo $jobid
}

# Note - reads must be gunzipped before running potrim.
run_potrim() 
{
    reads_r1=$1
    reads_r2=$2
    outname=$3
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    potrim_cmd="perl /opt/bio/popoolation/basic-pipeline/trim-fastq.pl --fastq-type sanger --input1 $reads_r1 --input2 $reads_r2 --quality-threshold 20 --min-length 50 --output $outname.trim"
    jobname="potrim"
    jobid=`run_qsub 1 $qsub_holdid "$potrim_cmd" $jobname`
    echo $jobid
}

run_read_subsample() 
{
    reads_r1=$1
    reads_r2=$2
    frac=$3 # fraction to keep - between 0 and 1.
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    reads_out1="${reads_r1}.s$frac.fq"
    reads_out2="${reads_r2}.s$frac.fq"
    qsub_holdid=1
    sample_cmd="read_sample.py $frac $reads_r1 $reads_r2 $reads_out1 $reads_out2"
    jobname="read_sample"
    jobid=`run_qsub 1 $qsub_holdid "$sample_cmd" $jobname`
    echo $jobid
}

# This function runs only one reads file, although an arbitrary list of reads
# files could be provided, possibly as a file listing instead so the number of args
# is always known.
run_read_stats() {
    reads=$1
    stats_file=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    stats_cmd="read_stats.py --reads $reads --stats_file $stats_file --stats_header"
    jobname="read_stats"
    jobid=`run_qsub 1 $qsub_holdid "$stats_cmd" $jobname`
    echo $jobid
}    

# join a few of the functs above.
run_fxtrim_fastqc_stats() {
    reads_in=$1
    reads_out=$2
    first_base=$3
    last_base=$4
    qsub_holdid=1
    [ ! -z $5 ] && qsub_holdid=$5
    fxjid=`run_fastx_trimmer $reads_in $reads_out $first_base $last_base $qsub_holdid`
    fqcjid=`run_fastqc $reads_out $fxjid`
    stats_file=$reads_out.stats
    fstatsjid=`run_read_stats $reads_out $stats_file $fqcjid`
    echo $fstatsjid
}

# Run quake
# Reverse-complement a read library
revcomp() {
    reads_in=$1
    # http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
    filename=$(basename "$reads_in")
    extension="${filename##*.}"
    get_qsub_script
    if [[ $extension = "gz" ]]; then
        reads_in_no_ext="${filename%.*}"
        reads_unlinked=`readlink -f $reads_in`
        qsub_gunzip_cmd="qsub -N gunzip -pe orte 1 qsub_script.sh \"gunzip_keep.sh $reads_unlinked ${reads_in_no_ext}\""
        echo $qsub_gunzip_cmd
        qsub_gunzip_out=`eval $qsub_gunzip_cmd`
        echo $qsub_gunzip_out
        qsub_gunzip_jobid=`echo $qsub_gunzip_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
        revcomp_cmd="/opt/bio/fastx/bin/fastx_reverse_complement -Q 33 -i ${reads_in_no_ext} -o rev_${reads_in_no_ext}"
        qsub_revcomp_cmd="qsub -N revcomp -pe orte 1 -hold_jid $qsub_gunzip_jobid qsub_script.sh \"$revcomp_cmd\""
        echo $qsub_revcomp_cmd
        qsub_revcomp_out=`eval $qsub_revcomp_cmd`
        echo $qsub_revcomp_out
        qsub_revcomp_jobid=`echo $qsub_revcomp_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
        qsub_gzip_cmd="qsub -N gzip -pe orte 1 -hold_jid $qsub_revcomp_jobid qsub_script.sh \"gzip rev_${reads_in_no_ext}\""
        echo $qsub_gzip_cmd
        qsub_gzip_out=`eval $qsub_gzip_cmd`
        echo $qsub_gzip_out
    else
        revcomp_cmd="/opt/bio/fastx/bin/fastx_reverse_complement -Q 33 -i $reads_in -o rev_${reads_in}"
        qsub_revcomp_cmd="qsub -N revcomp -pe orte 1 -hold_jid $qsub_gz_out qsub_script.sh \"$revcomp_cmd\""
        echo $qsub_revcomp_cmd
        qsub_revcomp_out=`eval $qsub_revcomp_cmd`
        echo $qsub_revcomp_out
        qsub_revcomp_jobid=`echo $qsub_revcomp_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
        qsub_gzip_cmd="qsub -N gzip -pe orte 1 -hold_jid $qsub_revcomp_jobid qsub_script.sh \"gzip rev_${reads_in}\""
        qsub_gzip_out=`eval $qsub_gzip_cmd`
        echo $qsub_gzip_out       
    fi
}
