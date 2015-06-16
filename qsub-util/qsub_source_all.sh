#!/bin/bash

SCRIPT=`readlink -f $BASH_SOURCE`
SCRIPTDIR=`dirname $SCRIPT`
SCRIPTLIST="qsub_utils.sh
qsub_assembly_qc.sh
qsub_misc.sh
qsub_read_qc.sh
qsub_tuxedo.sh
qsub_velvet.sh
qsub_zip_utils.sh"

for script in $SCRIPTLIST; do
	source $SCRIPTDIR/$script
done

