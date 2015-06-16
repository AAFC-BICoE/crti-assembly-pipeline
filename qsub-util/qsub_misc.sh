SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

# bug 4094
run_spades_PE() {
    R1=$1
    R2=$2
    prefix=spades_out_PE
    submit_host=biocomp-0-9
    [ ! -z $3 ] && prefix=$3
    get_qsub_script
    # SPAdes uses 16 threads and max 250G by default.
    cmd="qsub -N SPAdes_PE -pe smp 16 -l h=$submit_host -l mem_free=250G qsub_script.sh \"spades.py --careful --pe1-1 $R1 --pe1-2 $R2 -o $prefix\""
    echo $cmd
    eval $cmd
}

#run_spades_PE WG2-Ov-LEV6574-BC6_S1_L001_R1_001.fastq.gz WG2-Ov-LEV6574-BC6_S1_L001_R2_001.fastq.gz Se_LEV6574_spades

# Run with mate-pair data
# Note that insert sizes are not specified.
run_spades_MP() {
    PE_R1=$1
    PE_R2=$2
    MP3_R1=$3
    MP3_R2=$4
    MP8_R1=$5
    MP8_R2=$6
    prefix=spades_out_P38
    submit_host=biocomp-0-10
    [ ! -z $7 ] && prefix=$7
    get_qsub_script
    spades_cmd="spades.py --careful --pe1-1 $PE_R1 --pe1-2 $PE_R2 --mp1-1 $MP3_R1 --mp1-2 $MP3_R2 --mp2-1 $MP8_R1 --mp2-2 $MP8_R2 -o $prefix"
    # SPAdes uses 16 threads and max 250G by default
    qsub_cmd="qsub -N SPAdes_MP -pe smp 16 -l h=$submit_host -l mem_free=250G qsub_script.sh \"$spades_cmd\""
    echo $qsub_cmd
    eval $qsub_cmd
}     

# Run dipSPAdes
run_dipspades_PE() {
    R1=$1
    R2=$2
    prefix=dipspades_out_PE
    submit_host=biocomp-0-9
    [ ! -z $3 ] && prefix=$3
    get_qsub_script
    # SPAdes uses 16 threads and max 250G by default.
    cmd="qsub -N SPAdes_PE -pe smp 16 qsub_script.sh \"dipspades.py --pe1-1 $R1 --pe1-2 $R2 -o $prefix\""
    echo $cmd
    eval $cmd
}

#run_dipspades_PE WG2-Ov-LEV6574-BC6_S1_L001_R1_001.fastq.gz WG2-Ov-LEV6574-BC6_S1_L001_R2_001.fastq.gz Se_LEV6574_spades

# Run with mate-pair data
# Note that insert sizes are not specified.
run_dipspades_MP() {
    PE_R1=$1
    PE_R2=$2
    MP3_R1=$3
    MP3_R2=$4
    MP8_R1=$5
    MP8_R2=$6
    prefix=dipspades_out_P38
    submit_host=biocomp-0-10
    [ ! -z $7 ] && prefix=$7
    get_qsub_script
    dipspades_cmd="dipspades.py --pe1-1 $PE_R1 --pe1-2 $PE_R2 --mp1-1 $MP3_R1 --mp1-2 $MP3_R2 --mp2-1 $MP8_R1 --mp2-2 $MP8_R2 -o $prefix"
    # SPAdes uses 16 threads and max 250G by default
    qsub_cmd="qsub -N SPAdes_MP -pe smp 16 qsub_script.sh \"$dipspades_cmd\""
    echo $qsub_cmd
    eval $qsub_cmd
}     
