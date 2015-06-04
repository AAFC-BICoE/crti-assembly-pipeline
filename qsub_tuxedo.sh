SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

run_bowtie2_build() {
    genome=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    get_qsub_script
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

# From bug 4255
run_mpileup() {
    ref_genome=$1
    bamfile=$2
    outfile=$3
    qsub_holdid=1
    [ ! -z $4 ] && qsub_holdid=$4
    
    mpileup_cmd="samtools mpileup -Q0 -f ${ref_genome} $bamfile >$outfile"
    qsub_mpileup_cmd="qsub -N mpileup -pe orte 1 -hold_jid ${qsub_holdid} qsub_script.sh \"{mpileup_cmd}\""
    >&2 echo ${qsub_mpileup_cmd}
    qsub_mpileup_out=`eval ${qsub_mpileup_cmd}`
    >&2 echo ${qsub_mpileup_out}
    qsub_mpileup_jobid=`echo $qsub_index_bam_out | perl -ne 'if (/Your job ([0-9]+)/) { print $1 }'`
    echo ${qsub_mpileup_jobid}
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

