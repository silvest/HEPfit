#!/bin/sh -f
#BSUB -a openmp
#BSUB -q plong
#BSUB -J name
#BSUB -o output%J.txt
#BSUB -e error%J.txt
#          Specify the number of nodes requested and the
#          number of processors per node.
########BSUB -n 6
##########################################
# COMMAND to EXECUTE:
/teo/soft/silvest/bin/mpiexec -n 10 /teo/soft/silvest/HEPfit/dist/bin/analysis StandardModel.conf MonteCarlo.conf
