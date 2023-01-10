#!/bin/bash -l
# -l option is interactive login shell
#$ -pe omp 4
# Load the newest version of matlab on SCC
module load matlab
# Additional qsub options here . . .
matlab -nodisplay -r "runBatchfNIRSJob($NSLOTS); exit"
# matlab -nodisplay -r "runBatchfNIRSJob($1); exit"
