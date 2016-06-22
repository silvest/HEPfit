#!/bin/bash

#$ -M apaul2@nd.edu	 # Email address for job notification
#$ -m ae		 # Send mail when job begins, ends and aborts
#$ -pe mpi-* 48 	 # Specify parallel environment and legal core size
#$ -q long	         # Specify queue
#$ -N FOLDER_TAG             # Specify job name

#module load ompi/1.6.5-gcc # Required modules
#module load gsl          # Required modules
#module load boost         # Required modules

mpiexec -n $NSLOTS ./analysis conf_files/StandardModel.conf conf_files/MonteCarlo.conf --job_tag _FOLDER_TAG # Application to execute
