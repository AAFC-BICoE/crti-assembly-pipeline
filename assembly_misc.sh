# Miscellaneous functions, useful in a genome assembly context.
# Bash and R

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
    qsub_cmd="qsub -N $jobname -pe smp $nprocs -hold_jid $qsub_holdid qsub_script.sh \"$cmd\""
    >&2 echo ${qsub_cmd}
    qsub_out=`eval ${qsub_cmd}`
    >&2 echo ${qsub_out}
    qsub_jobid=`echo $qsub_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_jobid}
} 

# Given arbitrary number of input jobids, submit a single dummy job that holds on them all.
qsub_dummy_hold() {
    hold_str=""
    for jobid in "$@"; do
        hold_str="${hold_str} -hold_jid $jobid"
    done
    qsub_cmd="qsub -N qsub_dummy -pe smp 1 $hold_str qsub_script.sh \"ls\""
    >&2 echo ${qsub_cmd}
    qsub_out=`eval ${qsub_cmd}`
    >&2 echo ${qsub_out}
    qsub_jobid=`echo $qsub_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_jobid}     
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

run_bowtie2_build() {
    genome=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    bowtie2_build_cmd="bowtie2-build $genome $genome"
    qsub_bowtie2_build_cmd="qsub -N bowtie2_build -pe orte 1 -hold_jid ${qsub_holdid} qsub_script.sh \"${bowtie2_build_cmd}\""
    >&2 echo $qsub_bowtie2_build_cmd
    qsub_bowtie2_build_out=`eval ${qsub_bowtie2_build_cmd}`
    >&2 echo $qsub_bowtie2_build_out
    qsub_bowtie2_build_jobid=`echo $qsub_bowtie2_build_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`    
    echo ${qsub_bowtie2_build_jobid}
}

# 
run_bowtie2() {
    reads_R1=$1
    reads_R2=$2
    genome=$3
    samfile=$4
    qsub_holdid=1
    [ ! -z $5 ] && qsub_holdid=$5
    nprocs=6
    get_qsub_script
    bowtie2_cmd="bowtie2 -x $genome -q -1 $reads_R1 -2 $reads_R2 -S $samfile --threads $nprocs"
    qsub_bowtie2_cmd="qsub -N bowtie2 -pe smp $nprocs -hold_jid ${qsub_holdid} qsub_script.sh \"${bowtie2_cmd}\""
    >&2 echo ${qsub_bowtie2_cmd}
    qsub_bowtie2_out=`eval ${qsub_bowtie2_cmd}`
    >&2 echo ${qsub_bowtie2_out}
    qsub_bowtie2_jobid=`echo $qsub_bowtie2_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_bowtie2_jobid}
}

run_sam2bam() {
    samfile=$1
    bamfile=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    sam2bam_cmd="samtools view -S $samfile -b  -o $bamfile"
    qsub_sam2bam_cmd="qsub -N sam2bam -pe orte 1 -hold_jid ${qsub_holdid} qsub_script.sh \"${sam2bam_cmd}\""
    >&2 echo ${qsub_sam2bam_cmd}
    qsub_sam2bam_out=`eval ${qsub_sam2bam_cmd}`
    >&2 echo ${qsub_sam2bam_out}
    qsub_sam2bam_jobid=`echo ${qsub_sam2bam_out} | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_sam2bam_jobid}
}

sort_bam() {
    bamfile_in=$1
    bamfile_out=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    sort_bam_cmd="samtools sort $bamfile_in $bamfile_out"
    qsub_sort_bam_cmd="qsub -N sort_bam -pe orte 1 -hold_jid ${qsub_holdid} qsub_script.sh \"${sort_bam_cmd}\""
    >&2 echo ${qsub_sort_bam_cmd}
    qsub_sort_bam_out=`eval ${qsub_sort_bam_cmd}`
    >&2 echo ${qsub_sort_bam_out}
    qsub_sort_bam_jobid=`echo $qsub_sort_bam_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_sort_bam_jobid}
}

index_bam() {
    bamfile=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    index_bam_cmd="samtools index $bamfile"
    qsub_index_bam_cmd="qsub -N index_bam -pe orte 1 -hold_jid ${qsub_holdid} qsub_script.sh \"${index_bam_cmd}\""
    >&2 echo ${qsub_index_bam_cmd}
    qsub_index_bam_out=`eval ${qsub_index_bam_cmd}`
    >&2 echo ${qsub_index_bam_out}
    qsub_index_bam_jobid=`echo $qsub_index_bam_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_index_bam_jobid}
}

insert_histogram() {
    sorted_bam=$1
    prefix=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    insert_hist_cmd="java -jar /opt/bio/picard-tools/CollectInsertSizeMetrics.jar I=${sorted_bam} O=${prefix}.insertmetrics HISTOGRAM_FILE=${prefix}.insert.pdf"
    qsub_insert_hist_cmd="qsub -N insert_hist -pe orte 1 -hold_jid $qsub_holdid qsub_script.sh \"${insert_hist_cmd}\""
    >&2 echo ${qsub_insert_hist_cmd}
    qsub_insert_hist_out=`eval ${qsub_insert_hist_cmd}`
    >&2 echo ${qsub_insert_hist_out}
    qsub_insert_hist_jobid=`echo ${qsub_insert_hist_out} | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_insert_hist_jobid}
}

# do it all
run_bowtie2_all() {
    reads_R1=$1
    reads_R2=$2
    genome=$3
    prefix=$4
    samfile=$prefix.sam
    bamfile=$prefix.bam
    bamfile_sort=${prefix}_sort.bam
    insert_hist_prefix=${prefix}_sort
    
    bowtie2_build_jid=`run_bowtie2_build $genome`
    bowtie2_jid=`run_bowtie2 $reads_R1 $reads_R2 $genome $samfile ${bowtie2_build_jid}`
    sam2bam_jid=`run_sam2bam $samfile $bamfile ${bowtie2_jid}`
    sort_bam_jid=`sort_bam $bamfile ${bamfile_sort} ${sam2bam_jid}`
    index_bam_jid=`index_bam ${bamfile_sort} ${sort_bam_jid}`
    insert_hist_jid=`insert_histogram ${bamfile_sort} ${insert_hist_prefix} ${index_bam_jid}`
}

run_tophat() {
    RNA_R1=$1
    RNA_R2=$2
    bowtie2_genome_index=$3
    prefix=$4
    qsub_holdid=1
    [ ! -z $5 ] && qsub_holdid=$5
    get_qsub_script
    tophat_cmd="/opt/bio/tophat/bin/tophat -p 12 -o $prefix --mate-inner-dist 100  ${bowtie2_genome_index} $RNA_R1 $RNA_R2"
    qsub_tophat_cmd="qsub -N tophat -pe smp 12 -hold_jid ${qsub_holdid} qsub_script.sh \"${tophat_cmd}\""
    >&2 echo ${qsub_tophat_cmd}
    qsub_tophat_out=`eval ${qsub_tophat_cmd}`
    >&2 echo ${qsub_tophat_out}
    qsub_tophat_jobid=`echo ${qsub_tophat_out} | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_tophat_jobid}
}

run_cufflinks() {
    accepted_hits_bam=$1
    prefix=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    get_qsub_script
    cufflinks_cmd="cufflinks -p 12 -o $prefix ${accepted_hits_bam}"
    qsub_cufflinks_cmd="qsub -N cufflinks -pe smp 12 -hold_jid ${qsub_holdid} qsub_script.sh \"${cufflinks_cmd}\""
    >&2 echo ${qsub_cufflinks_cmd}
    qsub_cufflinks_out=`eval ${qsub_cufflinks_cmd}`
    >&2 echo ${qsub_cufflinks_out}
    qsub_cufflinks_jobid=`echo $qsub_cufflinks_out} | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_cufflinks_jobid}
}

run_bowtie_tophat_cufflinks()
{
    RNA_R1=$1
    RNA_R2=$2
    genome=$3
    prefix=$4
    get_qsub_script
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    bowtie2_build_jid=`run_bowtie2_build $genome`
    tophat_jid=`run_tophat ${RNA_R1} ${RNA_R2} $genome ${prefix}_tophat ${bowtie2_build_jid}`
    cufflinks_jid=`run_cufflinks ${prefix}_tophat/accepted_hits.bam ${prefix}_cufflinks ${tophat_jid}`
}

# Additional bowtie-alignment funcs added in context of bug 4260
run_bowtie2_lib() {
    reads_R1=$1
    reads_R2=$2
    lib_type=$3
    genome=$4
    samfile=$5
    qsub_holdid=1
    [ ! -z $6 ] && qsub_holdid=$6
    nprocs=6
    get_qsub_script
    lib_args=""
    if [ $lib_type = "PE" ]; then
        lib_args=" --fr --minins 200 --maxins 400 "
    elif [ $lib_type = "MSPE" ]; then
        lib_args=" --fr --minins 200 --maxins 400 "
    elif [ $lib_type = "MP3" ]; then
        lib_args="  --rf --minins 2000 --maxins 4000 "
    elif [ $lib_type = "MP8" ]; then
        lib_args="  --rf --minins 6000 --maxins 10000 "
    fi
    bowtie2_cmd="bowtie2 -x $genome -q -1 $reads_R1 -2 $reads_R2 ${lib_args} -S $samfile --threads $nprocs"
    qsub_bowtie2_cmd="qsub -N bowtie2_${lib_type} -pe smp $nprocs -hold_jid ${qsub_holdid} qsub_script.sh \"${bowtie2_cmd}\""
    >&2 echo ${qsub_bowtie2_cmd}
    qsub_bowtie2_out=`eval ${qsub_bowtie2_cmd}`
    >&2 echo ${qsub_bowtie2_out}
    qsub_bowtie2_jobid=`echo $qsub_bowtie2_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_bowtie2_jobid}
}

run_bowtie2_all_lib() {
    reads_R1=$1
    reads_R2=$2
    lib_type=$3 # options are PE, MSPE, MP3, MP8 (MSPE=MiSeq PE, though params are same as PE).
    genome=$4
    prefix=$5
    samfile=$prefix.sam
    bamfile=$prefix.bam
    bamfile_sort=${prefix}_sort
    insert_hist_prefix=${prefix}_sort.bam
    #bowtie2_build_jid=`run_bowtie2_build $genome`
    bowtie2_jid=`run_bowtie2_lib $reads_R1 $reads_R2 $lib_type $genome $samfile ${bowtie2_build_jid}`
    sam2bam_jid=`run_sam2bam $samfile $bamfile ${bowtie2_jid}`
    sort_bam_jid=`sort_bam $bamfile ${bamfile_sort} ${sam2bam_jid}`
    index_bam_jid=`index_bam ${bamfile_sort} ${sort_bam_jid}`
    insert_hist_jid=`insert_histogram ${bamfile_sort} ${insert_hist_prefix} ${index_bam_jid}`
}

run_interpro() {
    maker_proteins_fasta=$1 # usually maker2/<genome_name>.all.maker.proteins.fasta
    out_prefix=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    get_qsub_script
    IPSDIR=interproscan-42.0
    mkdir $IPSDIR
    interpro_cmd="/isilon/biodiversity/pipelines/interproscan-5/interproscan.sh -b ${out_prefix} -f TSV,XML,GFF3,HTML -goterms -iprlookup -pa -i ${maker_proteins_fasta} --seqtype p"
    qsub_interpro_cmd="qsub -N interpro -pe orte 1 -hold_jid $qsub_holdid qsub_script.sh \"${interpro_cmd}\""
    >&2 echo $qsub_interpro_cmd
    qsub_interpro_out=`eval ${qsub_interpro_cmd}`
    >&2 echo $qsub_interpro_out
    qsub_interpro_jobid=`echo $qsub_interpro_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_interpro_jobid}
}

# Needed to get bowtie aligned reads .bam file into format usable by GATK toolkit.
add_readgroups() {
    bamfile=$1
    qsub_holdid=$1
    [ ! -z $2 ] && qsub_holdid=$2
    get_qsub_script
    bamfile_short="${bamfile%.*}"
    outfile="${bamfile_short}.addRG.bam"
    addrg_cmd="java -jar /opt/bio/picard-tools/AddOrReplaceReadGroups.jar I=$bamfile O=$outfile LB=LB PL=illumina PU=PU SM=SM"
    qsub_addrg_cmd="qsub -N addrg -pe smp 1 -hold_jid ${qsub_holdid} qsub_script.sh \"${addrg_cmd}\""
    >&2 echo ${qsub_addrg_cmd}
    qsub_addrg_out=`eval ${qsub_addrg_cmd}`
    >&2 echo ${qsub_addrg_out}
    qsub_addrg_jobid=`echo ${qsub_addrg_out} | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_addrg_jobid}
}

# index bam between func above and func below.

run_depthcov() {
    bamfile=$1
    ref_genome=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    depthcov_cmd="java -jar /opt/bio/GenomeAnalysisTK/GenomeAnalysisTK.jar -T DepthOfCoverage -I $bamfile -R ${ref_genome}"
    jobname="depthcov"
    jobid=`run_qsub 1 $qsub_holdid "$depthcov_cmd" $jobname`
    echo $jobid
}

# Note - need to be a bit careful about paired-end reads
# See 'Trimming paired-end reads' section of docs: https://cutadapt.readthedocs.org/en/latest/guide.html#basic-usage
# Forward and reverse adapter strings should be different.
run_cutadapt() {
    reads_file_in=$1
    reads_file_out=$2
    adapt_str=$3
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    cutadapt_cmd="/home/AAFC-AAC/cullisj/software/pyenv/versions/2.7.6/bin/cutadapt -a ${adapt_str} -o ${reads_file_out} ${reads_file_in}"
    jobname="cutadapt"
    jobid=`run_qsub 1 $qsub_holdid "$cutadapt_cmd" $jobname`
    echo $jobid
}

run_fastqc()
{
    reads=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    fastqc_cmd="/opt/bio/FastQC/fastqc $reads"
    jobname="fastqc"
    jobid=`run_qsub 1 $qsub_holdid "$fastqc_cmd" $jobname`
    echo $jobid
}

# Note: reads in must be unzipped reads
run_fastx_trimmer()
{
    reads_in=$1
    reads_out=$2
    first_base=$3
    last_base=$4
    qsub_holdid=1
    [ ! -z $5 ] && qsub_holdid=$5
    fastx_trimmer_cmd="/opt/bio/fastx/bin/fastx_trimmer -Q 33 -i $reads_in -o $reads_out -f $first_base -l $last_base"
    jobname="fastx_trimmer"
    jobid=`run_qsub 1 $qsub_holdid "$fastx_trimmer_cmd" $jobname`
    echo $jobid
}

# Note - reads must be gunzipped before running potrim.
run_potrim() {
    reads_r1=$1
    reads_r2=$2
    outname=$3
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    potrim_cmd="perl /opt/bio/popoolation/basic-pipeline/trim-fastq.pl --fastq-type sanger --input1 $reads_r1 --input2 $reads_r2 --quality-threshold 20 --min-length 50 --output $outname.trim"
    jobname="potrim"
    jobid=`run_qsub 1 $qsub_holdid "$potrim_cmd" $jobname`
    echo $jobid
}


run_gunzip() {
    gzipped_in=$1
    gunzipped_out=$2
    qsub_holdid=1
    [ ! -z $3 ] && qsub_holdid=$3
    gunzip_cmd="/home/AAFC-AAC/cullisj/scripts/gunzip_keep.sh $gzipped_in $gunzipped_out"
    jobname="gunzip"
    jobid=`run_qsub 1 $qsub_holdid "$gunzip_cmd" $jobname`
    echo $jobid
}

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
    [ $zipbin = "bzip2" ] && zext="bz2"
    tarzip_cmd="tar -cf $dir.tar.$zipext --use-compress-prog=$zipbin $dir/ $delstr"
    jobname="tar$zipbin"
    jobid=`run_qsub 1 $qsub_holdid "$tarzip_cmd" $jobname`
    echo $jobid
}
    

# Convert cuff GTF to GFF + release renaming steps.
    
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
