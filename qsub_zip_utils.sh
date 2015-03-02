SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

# tar/zip a dir. first param is dirname.
# must pass zip program (e.g. 'gzip', 'bzip2'. 
# thrid param "del" will delete the original dir, "keep" to keep
# optional fourth param is a qsub jobid to hold on.
zip_dir () {
    dir="${1%/}"
    zipbin=$2
    delstr=""
    [ $3 = "del" ] && delstr=" --remove-files "
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    zipext="gz"
    [ $zipbin = "bzip2" ] && zipext="bz2"
    tarzip_cmd="tar -cf $dir.tar.$zipext --use-compress-prog=$zipbin $dir/ $delstr"
    jobname="tar$zipbin"
    jobid=`run_qsub 1 $qsub_holdid "$tarzip_cmd" $jobname`
    echo $jobid
}

run_gunzip() {
    gzipped_in=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    gunzip_cmd="gunzip $gzipped_in"
    jobname="gunzip"
    jobid=`run_qsub 1 $qsub_holdid "$gunzip_cmd" $jobname`
    echo $jobid
}

# Run gunzip but keep the input .gz file as is.
run_gunzip_keep() {
    gzipped_in=$1
    gunzipped_out=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    gunzip_cmd="/home/AAFC-AAC/cullisj/scripts/gunzip_keep.sh $gzipped_in $gunzipped_out"
    jobname="gunzip_keep"
    jobid=`run_qsub 1 $qsub_holdid "$gunzip_cmd" $jobname`
    echo $jobid
}

run_bunzip2() {
    bzip2_in=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$3
    bunzip2_cmd="bunzip2 $bzip2_in"
    jobname=bunzip2
    jobid=`run_qsub 1 $qsub_holdid "$bunzip2_cmd" $jobname`
    echo $jobid
}

# gzip/bzip2 a file
run_zip() {
    infile=$1
    zipbin=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    zip_cmd="$zipbin $infile"
    jobid=`run_qsub 1 $qsub_holdid "$gzip_cmd" $zipbin`
    echo $jobid
}

run_pbzip2() {
    infile=$1
    nprocs=$2 # Max procs should probably be less than 64
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    pbzip2_cmd="pbzip2 -p$nprocs $infile"
    jobname="pbzip2"
    jobid=`run_qsub $nprocs $qsub_holdid "$pbzip2_cmd"`
    echo $jobid
}

#run_gzip() {
#    infile=$1
#    qsub_holdid=1
#    [ ! -z $2 ] && qsub_holdid=$2
#    gzip_cmd="gzip $infile"
#    jobname="gzip"
#    jobid=`run_qsub 1 $qsub_holdid "$gzip_cmd" $jobname`
#    echo $jobid
#}
