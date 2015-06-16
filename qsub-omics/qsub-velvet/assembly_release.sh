#!/bin/bash

# Do some steps to create the release.

# For now, it's pretty basic: just add the release name to the contigs file.

release_name=
contigs=
contigs_out=
func=
min_contig_len=
median_contig_len=
contig_stats="contigStats.txt"
quast_report="quast_results/latest/report.tsv"
genome_assembly_yaml="genome_assembly.yml"
while getopts "f:c:n:o:" opt; do
    case "${opt}" in
        c)
            contigs=${OPTARG}
            ;;
        o)
            contigs_out=${OPTARG}
            ;;       
        n)
            release_name=${OPTARG}
            ;;
        f)
            func=${OPTARG}
            ;;
    esac
done

usage() {
    echo "Usage: $0 -f <function> -c <contigs file> ..."
    exit 1;
}

#[[ ! -z $func && ! -z $contigs ]] || usage

rename_contigs() {
    rename_contigs.pl -i $contigs -n $release_name -o $contigs_out
}

run_quast() {
    [ -e ./qsub_script.sh ] || svn export http://biodiversity/svn/source/AssemblyPipeline/qsub_script.sh
    qsub -N quast qsub_script.sh "quast.py $contigs"
}

get_contig_stats() {
    [ -e ./contigStats.pl ] || svn export http://biodiversity/svn/source/misc_scripts/contigStats.pl
    ./contigStats.pl -c $contigs >$contig_stats
}

get_assembly_stats() {
    ./velvetReleaseStats.pl -q $quast_report -c $contig_stats -y $genome_assembly_yaml
}

$func;