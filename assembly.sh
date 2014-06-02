#reads_R1_in= # not to be set here - set these in user-defined config before sourcing this file.
#reads_R2_in=
reads_prefix="@M01696"
est_genome_size=20000000
insert_length=301
raw_readlen=${insert_length}
trim_start=5
trim_stop=175
trim_readlen=$((trim_stop-trim_start+1))

velvetk_cov=30
velvetk_path="velvetk.pl"
trim_kmer_start=77
trim_kmer_end=127
trim_kmer_step=4
raw_kmer_start=105
raw_kmer_end=145
raw_kmer_step=4
trim_velveth_bin="velveth_127"
trim_velvetg_bin="velvetg_127"
raw_velveth_bin="velveth_145"
raw_velvetg_bin="velvetg_145"

kmer_index_raw=
kmer_index_trim=
velvetg_kmer=
func=
config_file=
while getopts "c:f:i:j:k:" opt; do
    case "${opt}" in
        c)
            config_file=${OPTARG}
            ;;
        f)
            func=${OPTARG}
            ;;
        i)
            kmer_index_raw=${OPTARG}
            ;;
        j)
            kmer_index_trim=${OPTARG}
            ;;
        k)
            velvetg_kmer=${OPTARG} # Not implemented
            ;;
    esac
done

[[ ! -z $config_file && -s $config_file ]] && source $config_file

reads_R1=`pwd`/`basename $reads_R1_in`
reads_R2=`pwd`/`basename $reads_R2_in`
trim_range="${trim_start}-${trim_stop}"
raw_velvet_dir="velvet_raw"
trim_velvet_dir="velvet_trim_${trim_range}"
velvetk_raw_outfile="$raw_velvet_dir/velvetk_raw_best_kmer.txt"
velvetk_trim_outfile="$trim_velvet_dir/velvetk_trim_${trim_range}_best_kmer.txt"

get_trim_reads_fname()
{
    reads_fq=$1
    reads_base_fq=`basename $reads_fq`
    reads_base="${reads_base_fq%%.*}"
    trim_reads="`pwd`/${reads_base}_trim_${trim_range}.fq.gz" 
    echo $trim_reads
} 

reads_R1_trim=`get_trim_reads_fname $reads_R1`
reads_R2_trim=`get_trim_reads_fname $reads_R2`

get_numreads_fname()
{
    reads_fq=$1
    reads_base_fq=`basename $reads_fq`
    reads_base="${reads_base_fq%%.*}"
    numreads_fname="${reads_base}_numreads.txt"
    echo $numreads_fname
}

numreads_R1_file=`get_numreads_fname $reads_R1`
numreads_R2_file=`get_numreads_fname $reads_R2`
numreads_R1_trim_file=`get_numreads_fname $reads_R1_trim`
numreads_R2_trim_file=`get_numreads_fname $reads_R2_trim`

build_kmer_range()
{
    end_kmer=$1
    step_size=$2
    start_kmer=$((best_kmer-10*step_size))
    echo "$start_kmer,$end_kmer,$step_size"
}

build_raw_kmer_range()
{
    # 1. if start/end def. in config, use those
    # 2. if best kmer set in file, use that
    # 3. else: not set.
    
    if [ -s $velvetk_raw_outfile ]; then
        best_kmer=`cat $velvetk_raw_outfile`
        raw_kmer_start
        raw_kmer_range=`build_kmer_range $best_kmer`
    fi
}
    
#raw_kmer_start=
#raw_kmer_end=
#raw_kmer_step=
#raw_kmer_range=`build_raw_kmer_range`

raw_kmer_range="${raw_kmer_start},${raw_kmer_end},${raw_kmer_step}"
trim_kmer_range="${trim_kmer_start},${trim_kmer_end},${trim_kmer_step}"
tmprange=`echo $raw_kmer_range | tr , -` # avoid commas in filename
exp_cov_raw_file="$raw_velvet_dir/exp_cov_raw_kmer_${tmprange}.txt"
tmprange=`echo $trim_kmer_range | tr , -` 
exp_cov_trim_file="$trim_velvet_dir/exp_cov_trim_${trim_range}_kmer_${tmprange}.txt"

dir_setup()
{
    ln -s $reads_R1_in $reads_R1
    ln -s $reads_R2_in $reads_R2
    mkdir $raw_velvet_dir
    ln -s $reads_R1_in $raw_velvet_dir/
    ln -s $reads_R2_in $raw_velvet_dir/
    mkdir $trim_velvet_dir
    svn export http://biodiversity/svn/source/AssemblyPipeline/ExpKmerCov.pl
}

run_fastqc()
{
    reads=$1
    /opt/bio/FastQC/fastqc $reads
}

get_fastqc_dir()
{
    reads=$1
    reads_base=`basename $reads`
    fastqc_dir1="${reads_base%.*}_fastqc" # just remove the .gz
    fastqc_dir2="${reads_base%%.*}_fastqc" # remove the .fq as well as .gz
    if [ -e $fastqc_dir1 ]; then
        echo $fastqc_dir1;
    elif [ -e $fastqc_dir2 ]; then
        echo $fastqc_dir2;
    else
        echo ""
    fi
}

run_fastqc_raw()
{
    run_fastqc $reads_R1
    run_fastqc $reads_R2
}

view_fastqc_html_raw()
{
    fastqc_dir_R1=`get_fastqc_dir $reads_R1`
    fastqc_dir_R2=`get_fastqc_dir $reads_R2`
    report_html_R1="${fastqc_dir_R1}/fastqc_report.html"
    report_txt_R1="${fastqc_dir_R1}/fastqc_data.txt"
    #firefox $report_html_R1
    report_html_R2="${fastqc_dir_R2}/fastqc_report.html"
    report_txt_R2="${fastqc_dir_R2}/fastqc_data.txt"
    #firefox $report_html_R2
    firefox $report_html_R1 $report_html_R2
}

run_fastx_trim()
{
    reads_in_fq=$1
    reads_out_fqgz=$2
    reads_base_fq=`basename ${reads_in_fq}`
    ext="${reads_base_fq##*.}"
    #trim_reads=`get_trim_reads_fname $reads_fq`
    echo "Starting fastx_trim for $reads_in_fq"
    if [ $ext = "gz" ]; then
        gunzip -c $reads_in_fq | /opt/bio/fastx/bin/fastx_trimmer -f $trim_start -l $trim_stop -z -o $reads_out_fqgz -Q33
    else
        /opt/bio/fastx/bin/fastx_trimmer -f $trim_start -l $trim_stop -z -i $reads_in_fq -o $reads_out_fqgz -Q33
    fi
    echo "Done fastx_trim for $reads_fq"
}

run_fastx_trim_raw()
{
    run_fastx_trim $reads_R1 $reads_R1_trim
    run_fastx_trim $reads_R2 $reads_R2_trim
}

run_fastqc_trim()
{
    #reads_R1_trim=`get_trim_reads_fname $reads_R1`
    #reads_R2_trim=`get_trim_reads_fname $reads_R2`
    # above vars now set using function at start of script.
    # this allows for subsequent runs to not have to manually set these vars here.
    run_fastqc $reads_R1_trim
    run_fastqc $reads_R2_trim
}

view_fastqc_html_trim()
{
    fastqc_dir_R1=`get_fastqc_dir $reads_R1_trim`
    fastqc_dir_R2=`get_fastqc_dir $reads_R2_trim`
    report_html_R1=$fastqc_dir_R1/fastqc_report.html
    report_txt_R1=$fastqc_dir_R1/fastqc_data.txt
    #firefox $report_html_R1
    report_html_R2=$fastqc_dir_R2/fastqc_report.html
    report_txt_R2=$fastqc_dir_R2/fastqc_data.txt
    #firefox $report_html_R2
    firefox $report_html_R1 $report_html_R2
}

# Quickly pull up the first few reads to see what the read id prefix is for counting.
read_head()
{
    reads_fq=$1
    reads_base_fq=`basename $reads_fq`
    ext="${reads_base_fq##*.}"    
    if [ $ext = "gz" ]; then
        gunzip -c $reads_fq | head -10
    else
        head -10 $reads_fq
    fi
}

read_head_raw()
{
    echo ""
    echo "Read head for raw file: head -10 $reads_R1"
    read_head $reads_R1
    echo ""
    echo "Read head for raw file: head -10 $reads_R1"
    read_head $reads_R2
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
    else
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

read_counts_raw()
{
    echo "Read counts for file $reads_R1:"
    read_counts $reads_R1 $numreads_R1_file
    echo "Read counts for file $reads_R2:"
    read_counts $reads_R2 $numreads_R2_file
}

read_counts_trim()
{
    echo "Read counts for file $reads_R1_trim:"
    read_counts $reads_R1_trim $numreads_R1_trim_file
    echo "Read counts for file $reads_R2_trim:"
    read_counts $reads_R2_trim $numreads_R2_trim_file
}

# currently read lengths are static based on trim size 
# and therefore a simple subtraction can be done in 
# initial variable settings, so this func. is not implemented yet.
#get_read_lengths()
#{
#    reads_fq=$1
#}

run_velvetk()
{
    R1=$1
    R2=$2
    best=`$velvetk_path --size $est_genome_size --cov $velvetk_cov --best $R1 $R2`
    echo $best
}

run_velvetk_raw()
{
    best_kmer=`run_velvetk $reads_R1 $reads_R2`
    echo "$best_kmer"
    echo "$best_kmer" >$velvetk_raw_outfile
}

run_velvetk_trim()
{
    best_kmer=`run_velvetk $reads_R1_trim $reads_R2_trim`
    echo "$best_kmer"
    echo "$best_kmer" >$velvetk_trim_outfile
}

run_velveth()
{
    velveth_bin=$1
    kmer_range=$2
    R1=$3
    R2=$4
    R1_base=`basename $R1`
    R1_ext="${R1_base##*.}"
    R2_base=`basename $R2`
    R2_ext="${R2_base##*.}"
    fq_type=" -fastq "
    if [[ $R1_ext = "gz" && $R2_ext = "gz" ]]; then
        fq_type=" -fastq.gz "
    fi
    vh_cmd="$velveth_bin velvet $kmer_range -shortPaired -separate $fq_type $R1 $R2"
    echo $vh_cmd
    eval $vh_cmd
}

run_velveth_raw()
{
    cd $raw_velvet_dir
    run_velveth $raw_velveth_bin $raw_kmer_range $reads_R1 $reads_R2
    cd ..
}

run_velveth_trim()
{
    cd $trim_velvet_dir
    run_velveth $trim_velveth_bin $trim_kmer_range $reads_R1_trim $reads_R2_trim
    cd ..
}

run_exp_cov_raw()
{
    R1_rlen=$raw_readlen
    R2_rlen=$raw_readlen
    R1_numreads=`read_counts $reads_R1 $numreads_R1_file`
    R2_numreads=`read_counts $reads_R2 $numreads_R2_file`
    kmer_range=$raw_kmer_range
    outfile=$exp_cov_raw_file
    ./ExpKmerCov.pl --read_length $R1_rlen --num_reads $R1_numreads --read_length $R2_rlen --num_reads $R2_numreads --genome_size $est_genome_size --kmer_range $kmer_range >$outfile
}

run_exp_cov_trim()
{
    R1_rlen=$trim_readlen
    R2_rlen=$trim_readlen
    R1_numreads=`read_counts $reads_R1_trim $numreads_R1_trim_file`
    R2_numreads=`read_counts $reads_R2_trim $numreads_R2_trim_file`
    kmer_range=$trim_kmer_range
    outfile=$exp_cov_trim_file
    ./ExpKmerCov.pl --read_length $R1_rlen --num_reads $R1_numreads --read_length $R2_rlen --num_reads $R2_numreads --genome_size $est_genome_size --kmer_range $kmer_range >$outfile
}

run_velvetg()
{
    velvetg_bin=$1
    kmer=$2
    exp_cov=$3
    velvetg_cmd="$velvetg_bin velvet_$kmer -amos_file yes -cov_cutoff auto -exp_cov $exp_cov -unused_reads yes -scaffolding yes -ins_length $insert_length"
    echo "$velvetg_cmd"
    eval $velvetg_cmd
}

run_velvetg_raw()
{
    kmer=
    exp_cov=
    if [[ ! -z $kmer_index_raw && -s $exp_cov_raw_file ]]; then
        #exp_cov=`awk -v kmer=$kmer '{if ($1==kmer) { print $2; }}' $exp_cov_raw_file`
        line=`awk -v idx=$kmer_index_raw '{if(FNR==idx) { print $0; }}' $exp_cov_raw_file`
        kmer=`echo $line | awk '{print $1}'`
        exp_cov=`echo $line | awk '{print $2}'`
    fi
    if [ -z $exp_cov ]; then
        exp_cov="auto"
    fi
    if [ ! -z $kmer ]; then
        cd $raw_velvet_dir
        run_velvetg $raw_velvetg_bin $kmer $exp_cov
        cd ..
    fi
}

run_velvetg_trim()
{
    kmer=
    exp_cov=
    if [[ -z $exp_cov && -s $exp_cov_trim_file ]]; then
        #exp_cov=`awk -v kmer=$kmer '{if ($1==kmer) { print $2; }}' $exp_cov_trim_file`
        line=`awk -v idx=$kmer_index_trim '{if(FNR==idx) { print $0; }}' $exp_cov_trim_file`
        kmer=`echo $line | awk '{print $1}'`
        exp_cov=`echo $line | awk '{print $2}'`        
    fi
    if [ -z $exp_cov ]; then
        exp_cov="auto"
    fi
    if [ ! -z $kmer ]; then
        cd $trim_velvet_dir
        run_velvetg $trim_velvetg_bin $kmer $exp_cov
        cd ..
    fi
}

assembly_stats()
{
    kmer=$1
    vdir="velvet_$kmer"
    log="$vdir/Log"
    ctg="$vdir/contigs.fa"
    printf "\n$kmer:\n"
    if [[ -e $vdir && -e $log && -e $ctg ]]; then
        echo "$log:"
        tail -2 $log
        echo "quast:"
        quast.py $ctg | grep "N50"
        echo "contigs:"
        grep -c "^>" $ctg
    else
        echo "Some assembly outputs are missing for kmer $kmer"
    fi
}

# Some combined functions below to perform tasks in larger chunks with simpler qsub ordering

# run dir_setup + raw fastqc
setup_raw_fastqc()
{
    dir_setup
    read_counts_raw
    run_fastqc_raw
}

# view raw fastqc results
# view_fastqc_html_raw

# Run trimming and trim fastqc
trim_fastqc()
{
    run_fastx_trim_raw
    read_counts_trim
    run_fastqc_trim
}

# view_fastqc_html_trim

# Velvetk + velveth
velvetkh_raw()
{
    run_velvetk_raw
    run_exp_cov_raw
    run_velveth_raw
    
}

velvetkh_trim()
{
    run_velvetk_trim
    run_exp_cov_trim
    run_velveth_trim
}

# Running velvetg on multiple kmers handled by assembly_qs.sh

# Finally, get stats for all output assemblies

[ ! -z `echo $func | awk '{print $1;}'` ] && eval $func
    
    
    
    
    
    
    