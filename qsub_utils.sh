get_qsub_script() {
    [ -e qsub_script.sh ] || svn export http://svn.biodiversity.agr.gc.ca/repo/source/AssemblyPipeline/qsub_script.sh
    sed -i 's/#$ -pe .*//' qsub_script.sh
}

run_qsub() {
    nprocs=$1
    qsub_holdid=$2
    cmd=$3
    jobname="run_qsub"
    [ ! -z $4 ] && jobname=$4
    get_qsub_script
    qsub_cmd="qsub -N $jobname -pe smp $nprocs -hold_jid $qsub_holdid -q all.q qsub_script.sh \"$cmd\""
    >&2 echo ${qsub_cmd}
    qsub_out=`eval ${qsub_cmd}`
    >&2 echo ${qsub_out}
    qsub_jobid=`echo $qsub_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_jobid}
} 

# Given arbitrary number of input jobids, submit a single dummy job (run ls) that holds on them all.
qsub_dummy_hold() {
    hold_str=""
    for jobid in "$@"; do
        hold_str="${hold_str} -hold_jid $jobid"
    done
    get_qsub_script
    qsub_cmd="qsub -N qsub_dummy -pe smp 1 $hold_str qsub_script.sh \"ls\""
    >&2 echo ${qsub_cmd}
    qsub_out=`eval ${qsub_cmd}`
    >&2 echo ${qsub_out}
    qsub_jobid=`echo $qsub_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_jobid}     
} 
