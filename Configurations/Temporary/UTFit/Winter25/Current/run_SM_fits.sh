#!/bin/sh  
  
#    for dir0 in `ls -d Fits_*`  
#    do  
#  
#	cd $dir0  
# SM fits  
#for dir0 in Standard_Model Standard_Model_DM
for dir0 in Standard_Model_DM
do
    if [ -d ${dir0} ]  
    then  
	cd ${dir0}
	
	# Full fits, tree
	for dir in `ls -d */`
	do
	    if [ ${dir} != "input/" ] && [ ${dir} != "Pred_Obs/" ] 
	    then  
		cd ${dir}  
		if [ ! -f MonteCarlo_plots.pdf ]  
		then
		    #-------------------------------------------------------------  
		    echo "submitting job in $PWD"  
		    #qsub -w `pwd` -vARGS="model.conf ../MonteCarlo-HIGH.conf" -N SM_FullFit /storage/local/home/theorm3/silvestrini/storage/HEPfit/Analysis_EW/submit_analysis_v1.sh  
		    rm -rf analysis core.*  
		    #ln -s ~/storage/HEPfit/Analysis/dist/Debug_MPI/GNU-Linux/analysis .  
		    bsub < submit_job.sh  
		else  
		    rm -f *.o* analysis core.*  
		    #-------------------------------------------------------------  
		fi  
		cd ..  
	    fi  
	done
	
	# All the rest: Indirect determinations, parametric uncertainties, predictions, MW_vs_X plots  
	
	#    for dir1 in Ind_determ Param_Unc Pred_Obs MW_mt_fits MW_sin2Eff_fits  
	for dir1 in Pred_Obs
	do  
	    
	    if [ -d ${dir1} ]  
            then  
		cd ${dir1}  
		for dir2 in `ls -d */`  
		do  
		    cd ${dir2}  
		    if [ ! -f MonteCarlo_plots.pdf ]  
		    then  
			#-------------------------------------------------------------  
			echo "submitting job in $PWD"  
			#qsub -w `pwd` -vARGS="model.conf ../../MonteCarlo-HIGH.conf" -N ${dir0}_${dir1}_${dir2} /storage/local/home/theorm3/silvestrini/storage/HEPfit/Analysis_EW/submit_analysis_v1.sh  
			rm -rf analysis core.*  
			#ln -s ~/storage/HEPfit/Analysis/dist/Debug_MPI/GNU-Linux/analysis .  
			bsub < submit_job.sh  
                    else  
			rm -f *.o* analysis core.*  
			#-------------------------------------------------------------  
		    fi  
		    cd ..  
		done  
		cd ..  
            fi  
	    
	done  
	cd ..  
    fi  
done
  
#cd ..  
  
#done  
#----------
