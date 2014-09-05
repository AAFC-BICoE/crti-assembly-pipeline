# Miscellaneous functions, useful in a genome assembly context.
# Bash and R

get_qsub_script() {
    [ -e qsub_script.sh ] || svn export http://biodiversity/svn/source/AssemblyPipeline/qsub_script.sh
    sed -i 's/#$ -pe .*//' qsub_script.sh
}

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

# 
run_bowtie() {
    reads_R1=$1
    reads_R2=$2
    genome=$3
    samfile=$4
    bamfile=$samfile.bam
    nprocs=6
    get_qsub_script
    bowtie2_build_cmd="bowtie2-build $genome $genome"
    qsub_bowtie2_build_cmd="qsub -N bowtie_build -pe orte 1 qsub_script.sh \"${bowtie2_build_cmd}\""
    echo $qsub_bowtie2_build_cmd
    qsub_bowtie2_build_out=`eval ${qsub_bowtie2_build_cmd}`
    echo $qsub_bowtie2_build_out
    qsub_bowtie2_build_jobid=`echo $qsub_bowtie2_build_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    bowtie2_cmd="bowtie2 -x $genome -q -1 $reads_R1 -2 $reads_R2 -S $samfile --threads $nprocs"
    qsub_bowtie2_cmd="qsub -N bowtie2 -pe smp $nprocs -hold_jid ${qsub_bowtie2_build_jobid} qsub_script.sh \"${bowtie2_cmd}\""
    echo $qsub_bowtie2_cmd
    qsub_bowtie2_out=`eval $qsub_bowtie2_cmd`
    qsub_bowtie2_jobid=`echo $qsub_bowtie2_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    sam2bam_cmd="samtools view -S $samfile -b $bamfile"
    qsub_sam2bam_cmd="qsub -N sam2bam -pe orte 1 -hold_jid ${qsub_bowtie2_jobid} qsub_script.sh \"${sam2bam_cmd}\""
    echo $qsub_sam2bam_cmd
    qsub_sam2bam_out=`eval $qsub_sam2bam_cmd`
    echo $qsub_sam2bam_out
}

# sort bam
# index bam

# Run quake


# R function to get n50 score
# block-comment-out
: << 'END'
#!/usr/bin/env Rscript
calc_n50 <- function(counts, genome_size = sum(counts)) {
    h = ceiling (genome_size/2)
    k = 1
    ksum = 0
    while (ksum < h) {
     ksum = ksum + counts[k]
     k = k + 1
    }
    print (paste("N50:", counts[k-1], "at contig number", k-1))
}

# Read in velvet stats.txt
stats <- read.table ("stats.txt", sep="\t", header=TRUE)
stats[,2] <- stats[,2]+89

y <- rev (sort (stats[,2]))
calc_n50(y)
END
