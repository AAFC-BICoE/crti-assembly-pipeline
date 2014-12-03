usage() { echo "Usage: $0 -c <config file> -f <function> -n <num kmers in range>" 1>&2; exit 1; }

qsub_script="velvet_qsub.sh"
[ -s $qsub_script ] || svn export http://svn.biodiversity.agr.gc.ca/repo/source/AssemblyPipeline/${qsub_script}
assembly_script="velvet_assembly.sh"

assembly_script_args=
config_file=
num_kmers=
func=
holdid=1
velveth_host=" -l h=biocomp-0-9 " # Host to submit potentially memory-intensive velveth jobs to
while getopts "a:c:f:h:n:m:" opt; do
    case "${opt}" in
        c)
            config_file=${OPTARG}
            ;;
        f)
            func=${OPTARG}
            ;;
        h)
            holdid=${OPTARG} 
            ;;
        k)
            kmer=${OPTARG} # Not yet implemented. Idea would be to qsub velvetg on a single kmer w exp cov auto.
            ;;
        n)
            num_kmers=${OPTARG}
            ;;
    esac
done

# Pass a function name and a set of job ids for the function to hold on.
submit_job() 
{
    fn=$1
    min_count=0
    count=0
    hold_jid_str=
    for j in "$@"; do
        if [ $count -gt $min_count ]; then
            hold_jid_str="${hold_jid_str} -hold_jid $j"
        fi
        ((count++))
    done
    host=
    if [[ $fn = "velvetkh_raw" || $fn = "velvetkh_trim" ]]; then
	host=$velveth_host
    fi
    qout=`qsub -N "${fn}" -pe smp 1 $hold_jid_str $host $qsub_script "$assembly_script -c $config_file -f $fn"`
    echo $qout 1>&2
    rjid=`echo $qout | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo $rjid
}

# Submit all velvetg jobs in range, asking for 12 procs per job.
# Done mainly to ensure adequate memory, not because all procs are needed.
submit_velvetg_smp12()
{
    fn=$1
    prev_holdid=$2
    num_kmers=$3
    name=$4
    [ -z $name ] && name="vg"
    hold_jid_str=
    if [ -z $num_kmers ]; then
        echo "Error: Must pass the number of kmers to run velvetg on. Ending." 1>&2
        exit 1;
    fi    
    if [ ! -z $prev_holdid ]; then
        hold_jid_str="${hold_jid_str} -hold_jid ${prev_holdid}"
    fi
    vg_jids=
    for kmer_index in `seq 1 $num_kmers`; do
        qscmd="qsub -N ${name}_${kmer_index} -pe smp 12 $hold_jid_str $qsub_script \"${assembly_script} $assembly_script_args -c $config_file -f $fn -i $kmer_index\""
        echo $qscmd 1>&2
        qout=`eval $qscmd`
        echo $qout 1>&2
        rjid=`echo $qout | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
        echo $rjid 1>&2
        vg_jids="$vg_ids $rjid"
    done
    # Submit a final job that holds on all of the above, and return its job id for
    # subsequent qsubs to hold on.
    jid=`submit_job dummy_velvetg $vg_jids`
    echo $jid
}

submit_all()
{
    prev_holdid=$1
    setup_jid=`submit_job run_setup $prev_holdid`
    velvetkh_jid=`submit_job velvetkh $setup_jid`
    velvetg_jid=`submit_velvetg_smp12 run_velvetg $velvetkh_jid $num_kmers vg`
    stats_jid=`submit_job combine_assembly_stats $velvetg_jid`
    echo $stats_jid
}

if [[ "$func" = "all" ]]; then
    submit_all $holdid
elif [ "$func" = "run_velvetg" ]; then 
    submit_velvetg_smp12 $func $holdid $num_kmers "vg";
else
    submit_job $func $holdid
fi


