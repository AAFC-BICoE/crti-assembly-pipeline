usage() { echo "Usage: $0 [ -a <assembly script args> -f <function> ]" 1>&2; exit 1; }

qsub_script="qsub_script.sh"
qsub_script_smp12="qsub_script_smp12.sh"
[ -s $qsub_script ] || svn export http://biodiversity/svn/source/AssemblyPipeline/qsub_script.sh
[ -s $qsub_script_smp12 ] || svn export http://biodiversity/svn/source/AssemblyPipeline/qsub_script_smp12.sh
assembly_script="assembly.sh"

assembly_script_args=
config_file=
num_kmers_raw=
num_kmers_trim=
func=
holdid=
while getopts "a:c:f:h:n:m:" opt; do
    case "${opt}" in
        a)
            # pass any number of args to the assembly.sh script, enclosed in quotes.
            assembly_script_args=${OPTARG}
            ;;
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
            kmer_raw=${OPTARG} # Not yet implemented. Idea would be to qsub velvetg on a single kmer w exp cov auto.
            ;;
        l)
            kmer_trim=${OPTARG} # Not yet implemented. Idea would be to qsub velvetg on a single kmer w exp cov auto.
            ;;
        n)
            num_kmers_raw=${OPTARG}
            ;;
        m)
            num_kmers_trim=${OPTARG}
            ;;
        r)
            kmer_range_raw=${OPTARG} # Not yet implemented. Idea would be to submit a range of vg cmds w exp cov auto
            ;;
        s)
            kmer_range_trim=${OPTARG} # Not yet implemented. Idea would be to submit a range of vg cmds w exp cov auto
            ;;
        x)
            exp_cov_file_raw=${OPTARG} # Not impl. Submit vg cmds using kmer, exp cov specified in pre-existing tab'd file.
            ;;
        y)
            exp_cov_file_trim=${OPTARG} # Not impl. Submit vg cmds using kmer, exp cov specified in pre-existing tab'd file.
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
    qout=`qsub -N "${fn}" $hold_jid_str $qsub_script "$assembly_script $assembly_script_args -c $config_file -f $fn"`
    echo $qout 1>&2
    rjid=`echo $qout | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo $rjid
}

# Submit all velvetg jobs in range to separate hosts.
submit_velvetg_host()
{
    fn=$1
    prev_holdid=$2
    num_kmers=$3
    if [ -z $num_kmers ]; then
        echo "Error: Must pass the number of kmers to run velvetg on. Ending."
        exit 1;
    fi
    hold_jid_str=
    if [ ! -z $holdid ]; then
        hold_jid_str=" -hold_jid_str $holdid "
    fi
    if [ ! -z $prev_holdid ]; then
        hold_jid_str="${hold_jid_str} -hold_jid ${prev_holdid}"
    fi
    count=0
    hostnum=1
    subhost="biocomp-0-1"
    for kmer_index in {1..$num_kmers}; do
        #until [ ! -z `qhost | grep $subhost` ]; do
        #    count=$((count+1))
        #    hostnum=$((count%11+1))
        #    subhost="biocomp-0-$hostnum"
        #done
        qscmd="qsub -N vg${kmer_index} -l h=$subhost $hold_jid_str $qsub_script \"$assembly_script $assembly_script_args -c $config_file -f $fn -i $kmer_index\""
        echo $qscmd
        eval $qscmd
        count=$((count+1))
        hostnum=$((count%11+1)) # Add the one at the end so we can hit biocomp-0-11
        subhost="biocomp-0-$hostnum"        
    done
}

# Submit all velvetg jobs in range, asking for 12 procs per job.
# Done mainly to ensure adequate memory, not because all procs are needed.
submit_velvetg_smp12()
{
    fn=$1
    prev_holdid=$2
    num_kmers=$3
    hold_jid_str=
    if [ -z $num_kmers ]; then
        echo "Error: Must pass the number of kmers to run velvetg on. Ending."
        exit 1;
    fi    
    if [ ! -z $holdid ]; then
        hold_jid_str=" -hold_jid_str $holdid "
    fi
    if [ ! -z $prev_holdid ]; then
        hold_jid_str="${hold_jid_str} -hold_jid ${prev_holdid}"
    fi
    for kmer_index in {1..$num_kmers}; do
        qscmd="qsub -N vg${kmer_index} $hold_jid_str $qsub_scrip_smp12 \"${assembly_script} $assembly_script_args -c $config_file -f $fn -i $kmer_index\""
        echo $qscmd
        #eval $qscmd
    done
}

# Submit all tasks via qsub.
submit_all() 
{
    raw_fastqc_jid=`submit_job setup_raw_fastqc`
    velvetkh_raw_jid=`submit_job velvetkh_raw $raw_fastqc_jid`
    velvetg_raw_jid=`submit_velvetg_smp12 run_velvetg_raw $velvetkh_raw_jid $num_kmers_trim`
    
    trim_fastqc_jid=`submit_job trim_fastqc $raw_fastqc_jid`
    velvetkh_trim_jid=`submit_job velvetkh_trim $trim_fastqc_jid`
    velvetg_trim_jid=`submit_velvetg_smp12 run_velvetg_trim $velvetkh_trim_jid $num_kmers_raw`
}

if [[ "$func" = "all" ]]; then
    submit_all
elif [ "$func" = "run_velvetg_raw" ]; then 
    submit_velvetg_smp12 $func 1 $num_kmers_raw;
elif [ "$func" = "run_velvetg_trim" ]; then
    submit_velvetg_smp12 $func 1 $num_kmers_trim;
else
    submit_job $func $holdid
fi


