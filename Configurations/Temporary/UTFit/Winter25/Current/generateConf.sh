#!/bin/bash
for Model in Standard_Model Standard_Model_DM
do
    if [ ! -d "${Model}" ] ; then
	mkdir ${Model}
    fi
    if [[ $Model == "Standard_Model_DM" ]] ; then
	cp StandardModel.conf Flavour.conf MonteCarlo.conf submit_job.sh UTfit_inputs.conf ${Model}/
	cp UTfit_DM.conf ${Model}/UTfit.conf
    else
	cp StandardModel.conf Flavour.conf UTfit.conf UTfit_inputs.conf MonteCarlo.conf submit_job.sh ${Model}/
    fi
    if [ ! -L "${Model}/input" ] ; then
	ln -s `pwd`/input ${Model}/input
    fi
    cd ${Model}
    if [[ $Model == "Standard_Model" ]] || [[  $Model == "Standard_Model_DM" ]]  ; then
	sed -i 's/ModelParameter\s\+'BBs'\([2-5]\s\+[0-9\.]\+\s\+\)[0-9\.]\+/ModelParameter  BBs\1 0./' UTfit.conf
	sed -i 's/ModelParameter\s\+'BBd'\([2-5]\s\+[0-9\.]\+\s\+\)[0-9\.]\+/ModelParameter  BBd\1 0./' UTfit.conf
	sed -i 's/ModelParameter\s\+'BK'\([2-5]\s\+[0-9\.]\+\s\+\)[0-9\.]\+/ModelParameter  BK\1 0./' UTfit.conf
    fi
    
    # Full Fit
    for Fit in Full_Fit Full_Fit_epe_th
    do
	if [ ! -d "${Fit}" ] ; then
	    mkdir ${Fit}
	fi
	cp StandardModel.conf Flavour.conf UTfit.conf MonteCarlo.conf submit_job.sh ${Fit}/
	if [ ! -L "${Fit}/input" ] ; then
	    ln -s `pwd`/input ${Fit}/input
	fi
	if [[ $Fit == "Full_Fit_epe_th" ]] ; then
	    sed -i 's/Observable\s\+EpsilonP_O_Epsilon \(.*\)MCMC\s\+\(weight\|file\)/Observable EpsilonP_O_Epsilon \1noMCMC noweight/' ${Fit}/UTfit.conf
	    sed -i 's/Observable\s\+EpsilonP_O_Epsilon_TH \(.*\)noMCMC\s\+\(noweight\|file\)/Observable EpsilonP_O_Epsilon_TH \1MCMC weight/' ${Fit}/UTfit.conf
	fi
	sed -i 's/-J name/-J '${Model}_${Fit}'/' ${Fit}/submit_job.sh
    done
    
    # Predictions for observables
    allObs="Dmd Dms EpsilonK alpha SJPsiK C2beta btaunu BRbar_Bsmumu EpsilonP_O_Epsilon beta"
    for Obs in $allObs
    do
	if [ ! -d "Pred_Obs/${Obs}" ] ; then
	    if [ ! -d "Pred_Obs" ] ; then
		mkdir Pred_Obs
	    fi
	    mkdir Pred_Obs/${Obs}
	fi
	cp StandardModel.conf Flavour.conf UTfit.conf MonteCarlo.conf submit_job.sh Pred_Obs/${Obs}/
	if [ ! -L "Pred_Obs/${Obs}/input" ] ; then
	    ln -s `pwd`/input Pred_Obs/${Obs}/input
	fi
	if [[ $Obs == "beta" ]] ; then
	    sed -i 's/Observable\s\+'SJPsiK'\(.*\)MCMC\s\+\(weight\|file\)/Observable 'SKPsiK'\1noMCMC noweight/' Pred_Obs/beta/UTfit.conf
	    sed -i 's/Observable\s\+'C2beta'\(.*\)MCMC\s\+\(weight\|file\)/Observable 'C2beta'\1noMCMC noweight/' Pred_Obs/beta/UTfit.conf
	else
	    sed -i 's/Observable\s\+'${Obs}'\(.*\)MCMC\s\+\(weight\|file\)/Observable '${Obs}'\1noMCMC noweight/' Pred_Obs/${Obs}/UTfit.conf
	fi
	sed -i 's/-J name/-J '${Model}_no${Obs}'/' Pred_Obs/${Obs}/submit_job.sh
    done
    # Predictions for CKM elements used as inputs
    declare -A CKMInputVals=( ["V_ud"]=0.97375 ["V_cb"]=0.0425 ["V_ub"]=0.0037 ["gamma"]=1.13)
    declare -A CKMInputErrs=( ["V_ud"]=0.004 ["V_cb"]=0.0025 ["V_ub"]=0.0005 ["gamma"]=.2)
    for Obs in "${!CKMInputErrs[@]}"
    do
	if [ ! -d "Pred_Obs/${Obs}" ] ; then
	    if [ ! -d "Pred_Obs" ] ; then
		mkdir Pred_Obs
	    fi
	    mkdir Pred_Obs/${Obs}
	fi
	cp StandardModel.conf Flavour.conf UTfit.conf MonteCarlo.conf submit_job.sh Pred_Obs/${Obs}/
	if [ ! -L "Pred_Obs/${Obs}/input" ] ; then
	    ln -s `pwd`/input Pred_Obs/${Obs}/input
	fi
	if [[ ${Obs} == "V_ub" ]] || [[ ${Obs} == "V_cb" ]] ; then
	    sed -i 's/CorrelatedGaussianParameters VubVcb 2/#/' Pred_Obs/${Obs}/UTfit.conf
	    sed -i 's/1\. 0\.11/#/' Pred_Obs/${Obs}/UTfit.conf
	    sed -i 's/0\.11 1\./#/' Pred_Obs/${Obs}/UTfit.conf
	fi
	sed -i 's/ModelParameter\s\+'${Obs}'\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter '${Obs}$'\t'${CKMInputVals[$Obs]}$'\t''0. '$'\t'${CKMInputErrs[$Obs]}'/' Pred_Obs/${Obs}/UTfit.conf
	sed -i 's/-J name/-J '${Model}_no${Obs}'/' Pred_Obs/${Obs}/submit_job.sh
    done
    # Predictions for lattice and other parameters
    allFits="BBs1BBsoBBdOnly FBsoFBdBBsoBBdOnly BKOnly"
    declare -A LatInpVals=( ["FBs"]=0.225 ["FBsoFBd"]=1.22 ["BBsoBBd"]=1.06 ["BBs1"]=0.84 ["BK1"]=0.62 )
    declare -A LatInpErrs=( ["FBs"]=0.025 ["FBsoFBd"]=0.15 ["BBsoBBd"]=0.15 ["BBs1"]=0.1 ["BK1"]=0.15 )
    declare -A BBs1BErrs=( ["FBs"]=0.035 ["FBsoFBd"]=0.15 ["BK1"]=0.25 )
    declare -A FBsBErrs=( ["FBs"]=0.05 ["BBs1"]=0.3 ["BK1"]=0.25)
    declare -A BKOnlyVals=( ["FBs"]=0.225 ["FBsoFBd"]=1.2 ["BBsoBBd"]=1.3 ["BBs1"]=0.9 ["BK1"]=0.62 )
    declare -A BKOnlyErrs=( ["FBs"]=0.1 ["FBsoFBd"]=0.5 ["BBsoBBd"]=1.0 ["BBs1"]=0.3 )
    for lfit in ${allFits}
    do
	if [ ! -d "Pred_Obs/${lfit}" ] ; then
            if [ ! -d "Pred_Obs" ] ; then
		mkdir Pred_Obs
            fi
            mkdir Pred_Obs/${lfit}
	fi
	cp StandardModel.conf Flavour.conf UTfit.conf MonteCarlo.conf submit_job.sh Pred_Obs/${lfit}/
	if [ ! -L "Pred_Obs/${lfit}/input" ] ; then
	    ln -s `pwd`/input Pred_Obs/${lfit}/input
	fi
    done
    for Obs in "${!LatInpErrs[@]}"
    do
	if [ ! -d "Pred_Obs/${Obs}" ] ; then
	    if [ ! -d "Pred_Obs" ] ; then
		mkdir Pred_Obs
	    fi
	    mkdir Pred_Obs/${Obs}
	fi
	cp StandardModel.conf Flavour.conf UTfit.conf MonteCarlo.conf submit_job.sh Pred_Obs/${Obs}/
	if [ ! -L "Pred_Obs/${Obs}/input" ] ; then
	    ln -s `pwd`/input Pred_Obs/${Obs}/input
	fi
	sed -i 's/ModelParameter\s\+'${Obs}'\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter '${Obs}$'\t'${LatInpVals[$Obs]}$'\t''0. '$'\t'${LatInpErrs[$Obs]}'/' Pred_Obs/${Obs}/UTfit.conf
	sed -i 's/-J name/-J '${Model}_no${Obs}'/' Pred_Obs/${Obs}/submit_job.sh
    done
    for Obs in  "${!BBs1BErrs[@]}"
    do
	sed -i 's/ModelParameter\s\+'${Obs}'\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter '${Obs}$'\t'' \1 '$'\t''0. '$'\t'${BBs1BErrs[$Obs]}'/' Pred_Obs/BBs1BBsoBBdOnly/UTfit.conf
    done
    sed -i 's/-J name/-J '${Model}_BBs1BBsoBBdOnly'/' Pred_Obs/BBs1BBsoBBdOnly/submit_job.sh
    for Obs in  "${!FBsBErrs[@]}"
    do
	sed -i 's/ModelParameter\s\+'${Obs}'\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter '${Obs}$'\t'' \1 '$'\t''0. '$'\t'${FBsBErrs[$Obs]}'/' Pred_Obs/FBsoFBdBBsoBBdOnly/UTfit.conf
    done
    sed -i 's/-J name/-J '${Model}_FBsoFBdBBsoBBdOnly'/' Pred_Obs/FBsoFBdBBsoBBdOnly/submit_job.sh
    for Obs in  "${!BKOnlyErrs[@]}"
    do
	sed -i 's/ModelParameter\s\+'${Obs}'\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter '${Obs}$'\t'${BKOnlyVals[$Obs]}$'\t''0. '$'\t'${BKOnlyErrs[$Obs]}'/' Pred_Obs/BKOnly/UTfit.conf
    done
    sed -i 's/-J name/-J '${Model}_BKOnly'/' Pred_Obs/BKOnly/submit_job.sh
#    for Obs in  "${!noLatErrs[@]}"
#    do
#	sed -i 's/ModelParameter\s\+'${Obs}'\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter '${Obs}$'\t'' \1 '$'\t''0. '$'\t'${noLatErrs[$Obs]}'/' Pred_Obs/noLattice/UTfit.conf
#    done
#    sed -i 's/-J name/-J '${Model}_noLattice'/' Pred_Obs/noLattice/submit_job.sh    
    # tree level and other fits
    allFits="tree UTfit-like angles noangles"
    for Fit in $allFits
    do
	if [ ! -d "${Fit}" ] ; then
	    mkdir ${Fit}
	fi
	cp StandardModel.conf Flavour.conf UTfit.conf MonteCarlo.conf submit_job.sh ${Fit}/
	if [ ! -L "${Fit}/input" ] ; then
	    ln -s `pwd`/input ${Fit}/input
	fi
	if [[ $Fit == "tree" ]] ; then
	    sed -i 's/Observable\s\+\(.*\)MCMC\s\+\(weight\|file\)/Observable \1noMCMC noweight/' ${Fit}/UTfit.conf
	elif  [[ $Fit == "UTfit-like" ]] ; then
	    allobs="BRbar_Bsmumu EpsilonP_O_Epsilon"
	    for Obs in ${allobs}
	    do
		sed -i 's/Observable\s\+'${Obs}' \(.*\)MCMC\s\+\(weight\|file\)/Observable '${Obs}' \1noMCMC noweight/' ${Fit}/UTfit.conf
	    done
	elif  [[ $Fit == "angles" ]] ; then
	    cp UTfit_inputs.conf ${Fit}/UTfit.conf
	    allobs="alpha_pipi alpha_rhopi alpha_rhorho gamma SJPsiK C2beta Phis_JPsiPhi"
	    for Obs in ${allobs}
	    do
		sed -i 's/Observable\s\+'${Obs}'\(.*\)noMCMC\s\+noweight/Observable '${Obs}'\1MCMC weight/' ${Fit}/UTfit.conf
		sed -i 's/Observable\s\+'${Obs}'\(.*\)noMCMC\s\+file/Observable '${Obs}'\1MCMC file/' ${Fit}/UTfit.conf
	    done
	elif  [[ $Fit == "noangles" ]] ; then
	    allobs="alpha_pipi alpha_rhopi alpha_rhorho SJPsiK C2beta Phis_JPsiPhi"
	    for Obs in ${allobs}
	    do
		sed -i 's/Observable\s\+'${Obs}'\(.*\)MCMC\s\+weight/Observable '${Obs}'\1noMCMC noweight/' ${Fit}/UTfit.conf
		sed -i 's/Observable\s\+'${Obs}'\(.*\)MCMC\s\+file/Observable '${Obs}'\1noMCMC noweight/' ${Fit}/UTfit.conf
	    done
	    sed -i 's/ModelParameter\s\+gamma\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)\s\+\([0-9\.\-]\+\)/ModelParameter gamma 1.570796 0. 1.570796/' ${Fit}/UTfit.conf  
	fi
	sed -i 's/-J name/-J '${Fit}'/' ${Fit}/submit_job.sh
    done
    # auxiliary fits for 2D UT plot
    allFits="Dmd Dms EpsilonK alpha SJPsiK gamma"
    for Fit in $allFits
    do
	if [ ! -d "Aux_${Fit}" ] ; then
	    mkdir Aux_${Fit}
	fi
	cp StandardModel.conf Flavour.conf UTfit_inputs.conf MonteCarlo.conf submit_job.sh Aux_${Fit}/
	if [ ! -L "Aux_${Fit}/input" ] ; then
	    ln -s `pwd`/input Aux_${Fit}/input
	fi
	if [[ $Fit == "Dmd" ]] ; then
	    sed -i 's/Observable\s\+'${Fit}'\(.*\)noMCMC\s\+\(noweight\|file\)/Observable '${Fit}'\1MCMC weight/' Aux_${Fit}/UTfit_inputs.conf
	elif [[ $Fit == "Dms" ]] ; then
	    sed -i 's/Observable\s\+'${Fit}'\(.*\)noMCMC\s\+\(noweight\|file\)/Observable '${Fit}'\1MCMC weight/' Aux_${Fit}/UTfit_inputs.conf
	elif [[ $Fit == "EpsilonK" ]] ; then
	    sed -i 's/Observable\s\+'${Fit}'\(.*\)noMCMC\s\+\(noweight\|file\)/Observable '${Fit}'\1MCMC weight/' Aux_${Fit}/UTfit_inputs.conf
	elif  [[ $Fit == "alpha" ]] ; then
	    allobs="alpha_pipi alpha_rhopi alpha_rhorho"
	    for Obs_i in ${allobs}
	    do
		sed -i 's/Observable\s\+'${Obs_i}'\(.*\)noMCMC\s\+\(weight\|file\)/Observable '${Obs_i}'\1MCMC file/' Aux_${Fit}/UTfit_inputs.conf
	    done
	elif [[ $Fit == "gamma" ]] ; then
	    sed -i 's/Observable\s\+'${Fit}'\(.*\)noMCMC\s\+\(noweight\|file\)/Observable '${Fit}'\1MCMC weight/' Aux_${Fit}/UTfit_inputs.conf
	fi
	sed -i 's/-J name/-J 'Aux_${Fit}'/' Aux_${Fit}/submit_job.sh
	cp Aux_${Fit}/UTfit_inputs.conf Aux_${Fit}/UTfit.conf
    done
    cd ..
done
