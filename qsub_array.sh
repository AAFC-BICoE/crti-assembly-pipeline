#!/bin/bash

#$ -S /bin/bash
#$ -V
# #$ -M $EMAIL
#$ -pe orte 1
#$ -cwd
export PATH=/usr/java/latest/bin/:$PATH
CMD=`/bin/sed -n ${SGE_TASK_ID}p $1`

echo Running on host: `hostname`.
echo In directory: `pwd`
echo Starting on: `date`

echo "$CMD"
$CMD

/bin/echo Finished on: `date`
