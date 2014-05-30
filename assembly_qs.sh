usage() { echo "Usage: $0 [ -a <assembly script args> -f <function> ]" 1>&2; exit 1; }

qsub_script="qsub_script.sh"
[ -s $qsub_script ] || svn export http://biodiversity/svn/source/AssemblyPipeline/qsub_script.sh
assembly_script="assembly.sh"

assembly_script_args=
config_file=
func=
exp_cov_file=
holdid=
while getopts "a:c:e:f:h:" opt; do
    case "${opt}" in
        a)
            # pass any number of args to the assembly.sh script, enclosed in quotes.
            assembly_script_args=${OPTARG}
            ;;
        c)
            config_file=${OPTARG}
            ;;
        e)
            exp_cov_file=${OPTARG}
            ;;
        f)
            func=${OPTARG}
            ;;
        h)
            holdid=${OPTARG}
            ;;
    esac
done

# Pass a function name and a set of job ids for the function to hold on.
submit_job() 
{
    fn=$1
    min_count=0
    if [[ $fn = "run_velvetg_raw" || $fn = "run_velvetg_trim" ]]; then
        fn="$1 $2"
        min_count=1
    fi
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

# Submit all velvetg jobs in range.
submit_velvetg()
{
    fn=$1
    if [ -z $exp_cov_file ]; then
        echo "Error: Must pass a file listing kmers/exp_cov values when running velvetg. Ending."
        exit 1;
    fi
    hold_jid_str=
    if [ ! -z $holdid ]; then
        hold_jid_str=" -hold_jid_str $holdid "
    fi
    count=0
    hostnum=1
    subhost="biocomp-0-1"
    while read line; do
        kmer=`echo $line | awk '{print $1}'`
        #until [ ! -z `qhost | grep $subhost` ]; do
        #    count=$((count+1))
        #    hostnum=$((count%11+1))
        #    subhost="biocomp-0-$hostnum"
        #done
        echo qsub -N vg${kmer} -l h=$subhost $hold_jid_str $qsub_script "$assembly_script $assembly_script_args -c $config_file -f $fn -k $kmer"
        qsub -N vg${kmer} -l h=$subhost $hold_jid_str $qsub_script "$assembly_script $assembly_script_args -c $config_file -f $fn -k $kmer"
        count=$((count+1))
        hostnum=$((count%11+1)) # Add the one at the end so we can hit biocomp-0-11
        subhost="biocomp-0-$hostnum"        
    done < $exp_cov_file
}

# Submit all tasks via qsub.
submit_all() 
{
    raw_fastqc_jid=`submit_job setup_raw_fastqc`
    velvetkh_raw_jid=`submit_job velvetkh_raw $raw_fastqc_jid`
    velvetg_raw_jid=`submit_job run_velvetg_raw $velvetkh_raw_jid`
    
    trim_fastqc_jid=`submit_job trim_fastqc $raw_fastqc_jid`
    velvetkh_trim_jid=`submit_job velvetkh_trim $trim_fastqc_jid`
    velvetg_trim_jid=`submit_job run_velvetg_trim $velvetkh_trim_jid`
}

if [[ "$func" = "all" ]]; then
    submit_all
elif [[ "$func" = "run_velvetg_raw" || "$func" = "run_velvetg_trim" ]]; then
    submit_velvetg;
else
    submit_job $func $holdid
fi


