#!/bin/sh -f 
#PBS -q sl7_mpi  
#          Specify the queue  
#PBS -N NAME  
#          Specify the number of nodes requested and the  
#          number of processors per node.  
#          PBS -l nodes=2:ppn=1,ncpus=12  
#PBS -l nodes=1:ppn=8,ncpus=8  
#          The directive below directs that the standard output and  
#          error streams are to be merged, intermixed, as standard  
#          output.  
#PBS -j oe  
#          Specify the walltime  
#PBS -l walltime=2400:00:00  
#          Specify the maximum amount of physical memory required per process.  
#          kb for kilobytes, mb for megabytes, gb for gigabytes.  
#          Take some care in setting this value.  Setting it too large  
#          can result in your job waiting in the queue for sufficient  
#          resources to become available.  
#PBS -l pmem=512mb  
  
########################################## 
  
#export MV2_ENABLE_AFFINITY=0  
NCPU=`wc -l < $PBS_NODEFILE`  
NNODES=`uniq $PBS_NODEFILE | wc -l`  
WORKDIR=$PBS_O_WORKDIR  
EXEFILE=/storage/local/home/theorm3/deblasmateo/storage/Software/HEPfit/Hf_Snowmass_2022/HEPfit/Analysis/dist/Debug_MPI/GNU-Linux/analysis  
RUNFILES="model.conf ../../MonteCarlo-HIGH.conf"  
  
echo ------------------------------------------------------  
echo ' This job is allocated on '${NCPU}' cpu(s)'  
echo 'Job is running on node(s): '  
cat $PBS_NODEFILE  
echo ------------------------------------------------------  
echo PBS: qsub is running on $PBS_O_HOST  
echo PBS: executing queue is $PBS_QUEUE  
echo PBS: working directory is $PBS_O_WORKDIR  
echo PBS: execution mode is $PBS_ENVIRONMENT  
echo PBS: job identifier is $PBS_JOBID  
echo PBS: job name is $PBS_JOBNAME  
echo PBS: node file is $PBS_NODEFILE  
echo PBS: number of nodes is $NNODES  
echo PBS: current home directory is $PBS_O_HOME  
echo PBS: PATH = $PBS_O_PATH  
echo ------------------------------------------------------  
  
  
SERVER=$PBS_O_HOST  
MACHINES=${WORKDIR}/NODEFILE  
LAUNCH="mpiexec -machinefile $PBS_NODEFILE -n $NCPU"  
  
echo server is $SERVER  
echo ------------------------------------------------------  
echo 'Job is running on node(s): '  
cat $PBS_NODEFILE  
echo ------------------------------------------------------  
echo ' '  
echo ' '  
  
 cd ${WORKDIR};  
  
# COMMAND to EXECUTE:  
 ${LAUNCH} ${EXEFILE} ${RUNFILES}  
  
exit
