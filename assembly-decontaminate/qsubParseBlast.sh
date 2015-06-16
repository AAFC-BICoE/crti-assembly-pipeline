#!/bin/bash

#$ -S /bin/bash
#$ -V
#$ -M $EMAIL
#$ -pe orte 1
#$ -cwd

genome_file=$1
blast_working_dir=`echo $genome_file | sed 's/\..*/_work/'`

blast_out_file=$genome_file.blast_out.txt
parsed_blast_file=$blast_out_file.parsed
sorted_blast_file=$parsed_blast_file.sorted

cat $blast_working_dir/tmp/*-result.txt > $blast_out_file
cat $genome_file.blast_out.txt | ./parseBlastOutput.pl -s -f -p > $parsed_blast_file
cat $parsed_blast_file | cut -f 4,14 | sort -k 2 -t$'\t' -nr > $sorted_blast_file
