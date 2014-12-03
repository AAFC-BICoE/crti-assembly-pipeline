#reads_R1_in= # not to be set here - set these in user-defined config before sourcing this file.
#reads_R2_in=
reads_prefix="@M01696" # should be able to determine this based on head -1
est_genome_size=25000000
insert_length=301
readlen=${insert_length}

velvetk_cov=30
velvetk_path="velvetk.pl"
kmer_start=65
kmer_end=127
kmer_step=4
kmer_rad= # Only used if we're using the best kmer from velvetk
velveth_bin="velveth_127"
velvetg_bin="velvetg_127"

velvetg_params=
kmer_index=
velvetg_kmer=
func=
config_file=
while getopts "c:f:i:k:" opt; do
    case "${opt}" in
        c)
            config_file=${OPTARG}
            ;;
        f)
            func=${OPTARG}
            ;;
        i)
            kmer_index=${OPTARG}
            ;;
        k)
            velvetg_kmer=${OPTARG} # Not implemented
            ;;
    esac
done

velvet_dir="velvet"

[[ ! -z $config_file && -s $config_file ]] && source $config_file

reads_R1=`pwd`/`basename $reads_R1_in`
reads_R2=`pwd`/`basename $reads_R2_in`
velvetk_outfile="$velvet_dir/velvetk_best_kmer.txt"

get_numreads_fname()
{
    reads_fq=$1
    reads_base_fq=`basename $reads_fq`
    reads_base="${reads_base_fq%.*}"
    numreads_fname="${reads_base}_numreads.txt"
    echo $numreads_fname
}

numreads_R1_file=`get_numreads_fname $reads_R1`
numreads_R2_file=`get_numreads_fname $reads_R2`

build_kmer_range()
{
    end_kmer=$1
    step_size=$2
    start_kmer=$((best_kmer-10*step_size))
    echo "$start_kmer,$end_kmer,$step_size"
}

build_kmer_range()
{
    # 1. if start/end def. in config, use those
    # 2. if best kmer set in file, use that
    # 3. else: not set.
    
    if [ -s $velvetk_outfile ]; then
        kmer_best=`cat $velvetk_outfile`
        kmer_start=$((kmer_best-kmer_rad*kmer_step))
        kmer_end=$((kmer_best+kmer_rad*kmer_step))
        kmer_range="${kmer_start},${kmer_end},${kmer_step}"
    fi
}
    
kmer_range="${kmer_start},${kmer_end},${kmer_step}"
tmprange=`echo $kmer_range | tr , -` # replace commas with dashes in filename
exp_cov_file="$velvet_dir/exp_cov_kmer_${tmprange}.txt"

dir_setup()
{
    ln -s $reads_R1_in $reads_R1
    ln -s $reads_R2_in $reads_R2
    mkdir $velvet_dir
    ln -s $reads_R1_in $velvet_dir/
    ln -s $reads_R2_in $velvet_dir/
    svn export -q http://svn.biodiversity.agr.gc.ca/repo/source/AssemblyPipeline/ExpKmerCov.pl
}

read_counts()
{
    reads_fq=$1
    numreads_fname=$2
    [ ! -z $numreads_fname ] || numreads_fname="/dev/null"
    if [ -z $numreads_fname ]; then
        numreads_fname="/dev/null"
    fi
    count=""
    if [[ -e $numreads_fname && -s $numreads_fname ]]; then
        count=`cat $numreads_fname`
    fi
    
    if [[ -z $count || $count -eq 0 ]]; then
        reads_base_fq=`basename $reads_fq`
        ext="${reads_base_fq##*.}"
        if [ $ext = "gz" ]; then
            count=`gunzip -c $reads_fq | grep -c "^${reads_prefix}"`
        else
            count=`grep -c "^${reads_prefix}" $reads_fq`
        fi
        echo "$count" >${numreads_fname}
    fi
    echo $count
}

run_read_counts()
{
    >&2 echo "Read counts for file $reads_R1:"
    read_counts $reads_R1 $numreads_R1_file
    >&2 echo "Read counts for file $reads_R2:"
    read_counts $reads_R2 $numreads_R2_file
}


run_velvetk()
{
    R1=$reads_R1
    R2=$reads_R2
    best_kmer=`$velvetk_path --size $est_genome_size --cov $velvetk_cov --best $R1 $R2`
    >&2 echo $best_kmer
    echo "$best_kmer" >$velvetk_outfile
}

run_velveth()
{
    R1=$reads_R1
    R2=$reads_R2
    cd $velvet_dir
    R1_base=`basename $R1`
    R1_ext="${R1_base##*.}"
    R2_base=`basename $R2`
    R2_ext="${R2_base##*.}"
    fq_type=" -fastq "
    if [[ $R1_ext = "gz" && $R2_ext = "gz" ]]; then
        fq_type=" -fastq.gz "
    fi
    vh_cmd="$velveth_bin velvet $kmer_range -shortPaired -separate $fq_type $R1 $R2"
    >&2 echo $vh_cmd
    eval $vh_cmd
}

run_exp_cov()
{
    R1_rlen=$readlen
    R2_rlen=$readlen
    R1_numreads=`read_counts $reads_R1 $numreads_R1_file`
    R2_numreads=`read_counts $reads_R2 $numreads_R2_file`
    outfile=$exp_cov_file
    ./ExpKmerCov.pl --read_length $R1_rlen --num_reads $R1_numreads --read_length $R2_rlen --num_reads $R2_numreads --genome_size $est_genome_size --kmer_range $kmer_range >$outfile
}

run_velvetg()
{
    kmer=
    exp_cov=
    if [[ ! -z $kmer_index && -s $exp_cov_file ]]; then
        line=`awk -v idx=$kmer_index '{if(FNR==idx) { print $0; }}' $exp_cov_file`
        kmer=`echo $line | awk '{print $1}'`
        exp_cov=`echo $line | awk '{print $2}'`
    fi
    if [ -z $exp_cov ]; then
        exp_cov="auto"
    fi
    if [ ! -z $kmer ]; then
        cd $velvet_dir
        velvetg_cmd="$velvetg_bin velvet_$kmer -very_clean yes -amos_file no -cov_cutoff auto -exp_cov $exp_cov -unused_reads yes -scaffolding yes -ins_length $insert_length $velvetg_params"
        >&2 echo "$velvetg_cmd"
        eval $velvetg_cmd
        cd ..
    fi
}

# Dummy function for assembly_qs.sh so it can perform jobid holding
dummy_velvetg()
{
    x=
}

assembly_stats()
{
    velvet_dir=$1
    start_dir=`pwd`
    if [[ -s "$velvet_dir/contigs.fa" && -s "$velvet_dir/Log" ]]; then
        cd $velvet_dir
        mkdir -p assembly_stats
        cd assembly_stats
        if [ ! -d quast_results ]; then
            quast.py "../contigs.fa" 2>&1 >/dev/null
        fi
        assembly_stats.pl --contig_file "../contigs.fa" --velvet_log "../Log" --quast_report_tsv "quast_results/latest/report.tsv" --name "$velvet_dir"
        cd $start_dir
    fi 
}

combine_assembly_stats()
{
    base_dir=$velvet_dir
    stats_outfile=$velvet_dir/assembly_stats_ts_combined.tab
    dir_list=`find $base_dir -name "contigs.fa" -exec dirname {} \; | sort | uniq`
    #printf "" >$stats_outfile # clear/touch the stats file
    assembly_stats.pl -h >$stats_outfile # Add just the header to the stats file.
    for d in $dir_list; do
        assembly_stats $d
        asm_stats_file="$d/assembly_stats/assembly_stats_transpose.tab"
        if [ -s $asm_stats_file ]; then
            cat $asm_stats_file >>$stats_outfile
        fi
    done
}


# run dir_setup + read counts
run_setup()
{
    dir_setup
    run_read_counts
}

# view fastqc results
# view_fastqc_html

# Velvetk + velveth
velvetkh()
{
    run_velvetk
    run_exp_cov
    run_velveth
}

# Running velvetg on multiple kmers handled by assembly_qs.sh

# Finally, get stats for all output assemblies

if [ ! -z `echo $func | awk '{print $1;}'` ]; then
    eval "$func";
fi
