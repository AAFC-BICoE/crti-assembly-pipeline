SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

run_quast() {
    assembly_fasta=$1
    qsub_holdid=1
    [ ! -z $2 ] && qsub_holdid=$2
    quast_cmd="quast.py $assembly_fasta"
    jobname="quast.py"
    jobid=`run_qsub 1 $qsub_holdid "$quast_cmd" $jobname`
    echo $jobid
}

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

