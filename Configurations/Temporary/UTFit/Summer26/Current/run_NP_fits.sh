#!/bin/sh
for dir0 in NPWC
do
    if [ -d ${dir0} ]
    then  
	cd ${dir0}
	for dir in `ls -d */`
	do
	    if [ ${dir} != "input/" ] && [ ${dir} != "Pred_Obs/" ] 
	    then  
		cd ${dir}  
#		if [ ! -f MonteCarlo_plots.pdf ]  
#		then
		    #-------------------------------------------------------------  
		    echo "submitting job in $PWD"  
		    qsub -vARGS="NPWC.conf ../MonteCarlo.conf" -N ${dir0}_${dir} ../submit_job.sh 
#		    bsub < submit_job.sh  
#		fi  
		cd ..  
	    fi  
	done
	cd ..  
    fi  
done
