#!/bin/bash
#SBATCH --nodes=3                       # number of nodes requested
#SBATCH --ntasks=240                    # number of tasks (default: 1 per core)
#SBATCH --ntasks-per-node=80            # number of tasks per node (default: whole node)
#SBATCH --partition=all                 # partition to run in
#SBATCH --job-name=NAME                 # job name
#SBATCH --output=NAME-%j-%N.out         # output file name
#SBATCH --error=NAME-%j-%N.err          # error file name
#SBATCH --time=48:00:00                 # runtime requested
#SBATCH --mail-user=ayan.paul@desy.de   # notification email
#SBATCH --mail-type=END,FAIL            # notification type
export LD_PRELOAD=""

# Define BAT
export BATINSTALLDIR=/home/apaul/NetBeansProjects/HEPfit-B1/BAT_parallel
export PATH="${BATINSTALLDIR}/bin:$PATH"
export LD_LIBRARY_PATH="${BATINSTALLDIR}/lib:$LD_LIBRARY_PATH"
export CPATH="${BATINSTALLDIR}/include:$CPATH"
# End Define BAT


# load modules:
module load mpi/openmpi-x86_64

# run the application:
mpiexec ./analysis config/model.conf config/MonteCarlo.conf
