#!/bin/bash

rawfile=$1
Qval=$2
fval=$3
outfile=$4

out_nogz=`echo $outfile | sed 's/\.gz//'`
/bin/gunzip -c $rawfile | /opt/bio/fastx/bin/fastx_trimmer -Q $Qval -f $fval -o $out_nogz
/bin/gzip $out_nogz



