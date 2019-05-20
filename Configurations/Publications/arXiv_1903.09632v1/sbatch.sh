#!/bin/bash
#SBATCH --nodes=6                       # number of nodes requested
#SBATCH --ntasks=240                    # number of tasks (default: 1)
#SBATCH --ntasks-per-node=40            # number of tasks per node (default: whole node)
#SBATCH --partition=maxwell             # partition to run in (all or maxwell)
#SBATCH --job-name=NAME                 # job name
#SBATCH --output=NAME-%N-%j.out         # output file name
#SBATCH --error=NAME-%N-%j.err          # error file name
#SBATCH --time=96:00:00                 # runtime requested
#SBATCH --mail-user=ayan.paul@desy.de   # notification email
#SBATCH --mail-type=END,FAIL            # notification type
#SBATCH --export=ALL

# load modules:
module load mpi/openmpi-x86_64

# run the application:
mpiexec ./analysis config/StandardModel.conf config/MonteCarlo.conf
