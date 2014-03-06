#!/bin/bash
qsub_path=/opt/gridengine/bin/lx26-amd64/qsub

# g is your contigs file
# p is number of processors per submitted job
# s is number of files to generate.
num_procs=8
num_splits=100
while getopts g:c:n: flag; do
    case $flag in
      g)
        genome_file=$OPTARG;
        ;;
      c)
        num_procs=$OPTARG;
        ;;
      n)
        num_splits=$OPTARG;
        ;;
      ?)
        exit;
        ;;
    esac
done
blast_working_dir=`echo $genome_file | sed 's/\..*/_work/'`
mkdir -p $blast_working_dir

# Setup and run the blast identification, given an input contigs file.
[ -e qsubParseBlast.sh ] || svn export http://biodiversity/svn/source/AssemblyPipeline/Decontaminate/qsubParseBlast.sh
[ -e parallel-blast.pl ] || svn export http://biodiversity/svn/source/BlastParallel/parallel-blast.pl
[ -e parseBlastOutput.pl] || svn export http://biodiversity/svn/source/misc_scripts/parseBlastOutput.pl

chmod +x parallel-blast.pl
chmod +x parseBlastOutput.pl

# Run parallel-blast.pl using the contigs file as input. 
merge_jobid=`./parallel-blast.pl -i $genome_file -d $blast_working_dir -p \
/opt/bio/ncbi-blast+/bin/blastn -db \
/isilon/biodiversity/reference/ncbi/blastdb/reference/nt/nt -o \
$genome_file.blast_out -c $num_procs -n $num_splits | grep 'Merge_results' \
| perl -ne 'if (/Your\sjob\s(\d+)/) { print $1; }'`

parse_jobid=`$qsub_path -N Parse_blast -hold_jid $merge_jobid ./qsubParseBlast.sh $genome_file`


