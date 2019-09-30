#!/bin/bash

# Load the necessary modules
module load intel/2019
module load mpi/openmpi-x86_64-intel

# $1 and $2 pick up the two arguments specified in submit_condor.sh
mpiexec -n 16 /afs/desy.de/user/a/apaul/Github/HiggsEW/WORK/WORK_DIR/analysis $1 $2
# mpiexec -n 16 ./analysis $1 $2

# Get the Observables directory back
cp -r Observables /afs/desy.de/user/a/apaul/Github/HiggsEW/WORK/WORK_DIR/.
tar zcvf Observables.tar.gz Observables
