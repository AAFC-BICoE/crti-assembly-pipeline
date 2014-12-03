SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
source $SCRIPTDIR/qsub_utils.sh

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

