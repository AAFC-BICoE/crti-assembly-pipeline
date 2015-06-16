#!/bin/bash

#$ -S /bin/bash
#$ -V
#$ -M $EMAIL
#$ -cwd
CMD=$1

/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

/bin/echo "Running command: ${CMD}"
$CMD

/bin/echo Finished on: `date`
