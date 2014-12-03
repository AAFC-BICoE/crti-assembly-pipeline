SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

# Run the velvet assembly qsub script on a pair of reads files
# Must provide a template config file
run_velvet_assembly() {
    config_file=$1
    num_kmers=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    velvet_cmd="velvet_qs.sh -c $config_file -f all -n $num_kmers -h $qsub_holdid"
    >&2 echo $velvet_cmd
    jobid=`$velvet_cmd`
    echo $jobid
}

write_velvet_config()
{
    config_filename=$1
cat >$config_filename <<EOF
reads_R1_in=$2
reads_R2_in=$3
kmer_start=$4
kmer_end=$5
kmer_step=$6
velvet_dir=$7
insert_length=$8
readlen=$9
reads_prefix=${10}
est_genome_size=${11}
EOF
}
    

# Given just the *two* reads files, create a velvet config file using defaults (not always optimal)
# and submit the velvet job
run_auto_velvet_assembly() {
    reads_R1=`readlink -f $1` # get full path to specified reads file
    reads_R2=`readlink -f $2`
    config_prefix=$3
    insert_length=$4
    readlen=$5
    est_genome_size=$6
    qsub_holdid=1
    [ ! -z $7 ] && qsub_holdid=$7
    num_kmers=20
    reads_prefix="@M01696"
    
    write_velvet_config $config_prefix.cfg \
        $reads_R1 $reads_R2 \
        31 127 4 \
        $config_prefix $insert_length $readlen \
        $reads_prefix $est_genome_size
    
    velvet_cmd="velvet_qs.sh -c $config_prefix.cfg -f all -n $num_kmers -h $qsub_holdid"
    >&2 echo $velvet_cmd
    jobid=`$velvet_cmd`
    
    echo $jobid
}

# Create output reads name given input reads name and trim info.
# Note: will always create an output name in same dir as input reads.
get_fxtrim_name()
{
    reads_in=$1
    first_base=$2
    last_base=$3
    reads_out=${reads_in%.*}_fxtrim_${first_base}_${last_base}
    echo $reads_out
}

# Perform auto-assembly, but trimming the input reads first.
run_fxtrim_auto_velvet_assembly() {
    reads_r1=`readlink -f $1` # get full path to specified reads
    reads_r2=`readlink -f $2`
    config_prefix=$3
    insert_length=$4
    # readlen calc'd from trim params provided below.
    est_genome_size=$5
    first_trim_base=$6
    last_trim_base=$7
    qsub_holdid=1
    [ ! -z $8 ] && qsub_holdid=$8
    num_kmers=20
    
    reads_r1_trim=`get_fxtrim name $reads_r1 $first_trim_base $last_trim_base`
    reads_r2_trim=`get_fxtrim name $reads_r2 $first_trim_base $last_trim_base`
    readlen=$((last_trim_base-first_trim_base))
    
    fxjid1=`run_fxtrim_fastqc_stats $reads_r1 $reads_r1_trim $first_trim_base $last_trim_base $qsub_holdid`
    fxjid2=`run_fxtrim_fastqc_stats $reads_r2 $reads_r2_trim $first_trim_base $last_trim_base $qsub_holdid`
    fxjid12=`qsub_dummy_hold $fxjid1 $fxjid2`    
    
cat >>$config_prefix.cfg <<EOF
reads_R1_in=$reads_r1
reads_R2_in=$reads_r2
kmer_start=31
kmer_end=127
kmer_step=4
velvet_dir=$config_prefix
insert_length=$insert_length
readlen=$readlen
reads_prefix="@M01696"
est_genome_size=$est_genome_size
EOF
    velvet_cmd="velvet_qs.sh -c $config_prefix.cfg -f all -n $num_kmers -h $fxjid12"
    >&2 echo $velvet_cmd
    #jobid=`$velvet_cmd`
    echo $jobid
}
