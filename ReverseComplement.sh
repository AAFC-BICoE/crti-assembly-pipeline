#!/bin/bash
datafile=$1
outfile=$2
out_nogz=`echo $outfile | sed 's/\.gz//'`
/bin/gunzip -c $datafile | /opt/bio/fastx/bin/fastx_reverse_complement -Q 33 -o $out_nogz
/bin/gzip $out_nogz



