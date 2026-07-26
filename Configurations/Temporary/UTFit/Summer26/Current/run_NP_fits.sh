#!/bin/sh
for dir0 in $1
do
    if [ -d ${dir0} ]
    then  
		cd ${dir0}
		if [[ "$dir0" == "NPDF2"* ]]; then
	    	echo "submitting job in $PWD"  
	    	qsub -vARGS="NPDF2.conf MonteCarlo.conf" -N ${dir0} submit_job.sh 
		else
			for dir in `ls -d */`
			do
	    		if [ ${dir} != "input/" ] && [ ${dir} != "Pred_Obs/" ] 
	    		then  
					cd ${dir}  
#		if [ ! -f MonteCarlo_plots.pdf ]  
#		then
		    #-------------------------------------------------------------  
		    		echo "submitting job in $PWD"  
		    		qsub -vARGS="NPWC.conf MonteCarlo.conf" -N ${dir0}_${dir} submit_job.sh 
#		    bsub < submit_job.sh  
#		fi  
					cd ..  
	    		fi  
			done
		fi
		cd ..  
    fi  
done
