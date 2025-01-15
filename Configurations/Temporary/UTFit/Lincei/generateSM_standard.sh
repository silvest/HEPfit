#!/bin/bash

# IMPORTANT: For the SM fits it is important to remove from the observables any of the "delta" observables, which are
# only defined in NP models

# Script needed to create all the fits needed to extract all the relevant
# info of the SM fit

# Runs from "SM_fits" folder containing:
# model.conf, ObservablesEW.conf, FlavourFixed.conf
# model.conf must call ObservablesEW.conf, FlavourFixed.conf at the same level

# WARNING: If correlation matrices change the corresponding lines in Section 4
# for the correlated Gaussian observables must be updated accordingly

# Run from Current_Future: The fits Folder will be created in Current_Future/Fits_Current with
# the name "SM_fits"

#------------------------------------------------------------------------------

# Fits to generate

# SM full fit
GenFullfit='Yes'

# Indirect determinations of SM inputs
GenIndDet='Yes'

# Parametric uncertainties
GenParUnc='Yes'

# Predictions for observables
GenPredObs='Yes'

# Special plots: MW, mt, sinEff2
GenMWmtsinEff='Yes'

# Lepton Flavor universal: 'Yes' use LFU for predictions. It also selects the right set of observables for the fit.
LFUobs='Yes'

#------------------------------------------------------------------------------

# EW observable files

ObsEWFU=ObservablesEW_Current_LFU.conf
ObsEWnoFU=ObservablesEW_Current_noLFU.conf

#------------------------------------------------------------------------------

# Inputs:

# Experimental and SM best fit values (mode) for the SM input parameters

alSexp=0.1177
alSSM=0.117606

dal5hexp=0.02766
dal5hSM=0.0275361

mtexp=171.79
mtSM=172.361

mHexp=125.21
mHSM=125.2055

MZexp=91.1875
MZSM=91.191099

# Experimental errors for the SM input parameters

alSexpErr=0.0010

dal5hexpErr=0.00010

mtexpErr=0.38

mHexpErr=0.12

MZexpErr=0.0021

#------------------------------------------------------------------------------

# Options

#     Include SM theoretical uncertainties
#     No -> No theory uncertainties
#     Current -> Use current theory uncertainties
#     Future -> Use future projections on theory uncertainties
ThU='Current'

# Choose if the intrinsic theory uncertainties must be used with 'Gaussian' of 'Flat' priors
IntUPrior='Gaussian'

# M_W [GeV]: Initialize with no th. unc.
MWthU=0.

# sin^2{\theta_{Eff}^\ell}: Initialize with no th. unc.
s2lEffthU=0.

# sin^2{\theta_{Eff}^q}: Initialize with no th. unc.
s2qEffthU=0.

# sin^2{\theta_{Eff}^b}: Initialize with no th. unc.
s2bEffthU=0.

# \Gamma_Z [GeV]: Initialize with no th. unc.
GammaZthU=0.

# \sigma_H [nb]: Initialize with no th. unc.
SigmaHthU=0.

# R_l^0: Initialize with no th. unc.
RlthU=0.

# R_c^0: Initialize with no th. unc.
RcthU=0.

# R_b^0: Initialize with no th. unc.
RbthU=0.

if [ $ThU == 'Current' ] || [ $ThU == 'current' ]
then

# M_W [GeV]: Current: 0.004      Future: 0.001
      MWthU=0.004

# sin^2{\theta_{Eff}^\ell}: Current: 0.000047      Future: 0.000015
      s2lEffthU=0.000047

# sin^2{\theta_{Eff}^q}: Current: 0.00005      Future: 0.000015?
      s2qEffthU=0.00005

# sin^2{\theta_{Eff}^b}: Current: 0.00005      Future: 0.000015?
      s2bEffthU=0.00005

# \Gamma_Z [GeV]: Current: 0.0004      Future: 0.00015
      GammaZthU=0.0004

# \sigma_H [nb]: Current: 0.006      Future: 0.000?
      SigmaHthU=0.006

# R_l^0: Current: 0.006      Future: 0.0015
      RlthU=0.006

# R_c^0: Current: 0.00005      Future: 0.0000?
      RcthU=0.00005

# R_b^0: Current: 0.00010      Future: 0.00005
      RbthU=0.0001

elif [ $ThU == 'Future' ] || [ $ThU == 'future' ]
then

# M_W [GeV]: Current: 0.004      Future: 0.001
      MWthU=0.001

# sin^2{\theta_{Eff}^\ell}: Current: 0.000047      Future: 0.000015
      s2lEffthU=0.000015

# sin^2{\theta_{Eff}^q}: Current: 0.00005      Future: 0.000015?
      s2qEffthU=0.000015

# sin^2{\theta_{Eff}^b}: Current: 0.00005      Future: 0.000015?
      s2bEffthU=0.000015

# \Gamma_Z [GeV]: Current: 0.0004      Future: 0.00015
      GammaZthU=0.00015

# \sigma_H [nb]: Current: 0.006      Future: 0.000?
      SigmaHthU=0.000

# R_l^0: Current: 0.006      Future: 0.0015
      RlthU=0.0015

# R_c^0: Current: 0.00005      Future: 0.0000?
      RcthU=0.0000

# R_b^0: Current: 0.00010      Future: 0.00005
      RbthU=0.00005

# Run also the script that modifies the top quark mass at FCCee-tt and ILC: NOT ANYMORE
#      sh fixmt_therr_FCC_ILC.sh

fi

#------------------------------------------------------------------------------

# Path of the HEPfit installation in the cluster (if going ro run there)
# This will go into a sed command so add the esc characters if using especial symbols

HEPfitPATH='\/storage\/local\/home\/theorm3\/deblasmateo\/storage\/Software\/HEPfit\/Hf_Current\/HEPfit\/Analysis\/dist\/Debug_MPI\/GNU-Linux\/analysis'

# Other Cluster settings for short and long fits

# Queue
QUEUESH='mpi_ib_short'
QUEUELG='mpi_ib'

# Wall time
WALLSH='450:00:00'
WALLLG='5000:00:00'

# Memory (512 is enough)
MEMSH=512
MEMLG=512

# Monte Carlo file to use (Adjust properties in Section 2)
MCCONFSH=MonteCarlo.conf
MCCONFLG=MonteCarlo-HIGH.conf

#------------------------------------------------------------------------------

# MonteCarlo file

# Contents of the MonteCarlo file (DO NOT TOUCH THIS. ADJUST ITERATIONS BELOW!!!!!!!!)
MCfile="## Number of chains \n\
NChains                    12 \n\
## Max iterations in prerun \n\
PrerunMaxIter              1000000 \n\
## Analysis iterations \n\
Iterations                 2000000 \n\
## Write Markov Chain \n\
WriteChain                 false \n\
## Use a particular seed \n\
#Seed                      0 \n\
## Find mode with Minuit \n\
FindModeWithMinuit         true \n\
## Calculate the evidence (total normalization) \n\
CalculateNormalization     false \n\
## Print all marginalized plots \n\
PrintAllMarginalized       true \n\
## Print correlation matrix \n\
PrintCorrelationMatrix     true \n\
## Print knowledge update plots \n\
PrintKnowledgeUpdatePlots  false \n\
## Print parameter plots \n\
PrintParameterPlot         true \n\
## Print triangle parameter plots \n\
PrintTrianglePlot          false \n\
## Set proposal to multivariate \n\
MultivariateProposal       true \n\
## Adjust minimum efficiency of chains (between 0. and 1.; def=0.15) \n\
MinimumEfficiency	   0.15 \n\
## Adjust maximum efficiency of chains (between 0. and 1.; def=0.5) \n\
MaximumEfficiency	   0.5 \n\
## Adjust R value for convergence of chains (> 1.0); def=1.1 \n\
RValueForConvergence    1.005 \n\
## Set the number of iterations between updates \n\
NIterationsUpdateMax       1000 \n\
## Smooth 1D histogram (true, false, 0-5) \n\
Histogram1DSmooth          false \n\
## Type of 2D Histogram 1001 -> Pixelated, 101 -> Filled, 1 -> Contour. \n\
Histogram2DType            			101 \n\
## Toggle to print legends for the histograms true -> no legend printed, false -> legend printed \n\
NoHistogramLegend          			true \n\
## Toggle to print HEPfit logo for the histograms true -> logo printed, false -> no logo printed \n\
PrintLogo                  			false \n\
## Adjust the opacity of 2D histograms (between 0. and 1.) \n\
Histogram2DAlpha            		1. \n\
## Use the average of each parameter as initial point \n\
MCMCInitialPosition    Center \n\
## Number of significant digits in Statistics file \n\
SignificantDigits 5 \n\
## The buffer size for histograms (0: default or positive integer) \n\
HistogramBufferSize           1000000 \n\
## Number of bins used for 1D and 2D observables \n\
NBinsHistogram1D	    400 \n\
NBinsHistogram2D	    400"

#------------------------------------------------------------------------------

# Contents of the Cluster submission file (DO NOT TOUCH THIS. ADJUST PARAMETERS BELOW DEPENDING ON FIT)
# Add esc characters for the special symbols
CSfile="#!/bin/sh -f \n\
#PBS -q mpi_ib  \n\
#          Specify the queue  \n\
#PBS -N NAME  \n\
#          Specify the number of nodes requested and the  \n\
#          number of processors per node.  \n\
#          PBS -l nodes=2:ppn=1,ncpus=12  \n\
#PBS -l nodes=1:ppn=12,ncpus=12  \n\
#          The directive below directs that the standard output and  \n\
#          error streams are to be merged, intermixed, as standard  \n\
#          output.  \n\
#PBS -j oe  \n\
#          Specify the walltime  \n\
#PBS -l walltime=2400:00:00  \n\
#          Specify the maximum amount of physical memory required per process.  \n\
#          kb for kilobytes, mb for megabytes, gb for gigabytes.  \n\
#          Take some care in setting this value.  Setting it too large  \n\
#          can result in your job waiting in the queue for sufficient  \n\
#          resources to become available.  \n\
#PBS -l pmem=512mb  \n\
  \n\
########################################## \n\
  \n\
#export MV2_ENABLE_AFFINITY=0  \n\
NCPU=\`wc -l < \$PBS_NODEFILE\`  \n\
NNODES=\`uniq \$PBS_NODEFILE | wc -l\`  \n\
WORKDIR=\$PBS_O_WORKDIR  \n\
EXEFILE=/storage/local/home/theorm3/deblasmateo/storage/Software/HEPfit/Master/HEPfit/Analysis/dist/Debug_MPI/GNU-Linux/analysis  \n\
RUNFILES=\"model.conf ../../MonteCarlo-HIGH.conf\"  \n\
  \n\
echo ------------------------------------------------------  \n\
echo ' This job is allocated on '\${NCPU}' cpu(s)'  \n\
echo 'Job is running on node(s): '  \n\
cat \$PBS_NODEFILE  \n\
echo ------------------------------------------------------  \n\
echo PBS: qsub is running on \$PBS_O_HOST  \n\
echo PBS: executing queue is \$PBS_QUEUE  \n\
echo PBS: working directory is \$PBS_O_WORKDIR  \n\
echo PBS: execution mode is \$PBS_ENVIRONMENT  \n\
echo PBS: job identifier is \$PBS_JOBID  \n\
echo PBS: job name is \$PBS_JOBNAME  \n\
echo PBS: node file is \$PBS_NODEFILE  \n\
echo PBS: number of nodes is \$NNODES  \n\
echo PBS: current home directory is \$PBS_O_HOME  \n\
echo PBS: PATH = \$PBS_O_PATH  \n\
echo ------------------------------------------------------  \n\
  \n\
  \n\
SERVER=\$PBS_O_HOST  \n\
MACHINES=\${WORKDIR}/NODEFILE  \n\
LAUNCH=\"mpiexec -machinefile \$PBS_NODEFILE -n \$NCPU\"  \n\
  \n\
echo server is \$SERVER  \n\
echo ------------------------------------------------------  \n\
echo 'Job is running on node(s): '  \n\
cat \$PBS_NODEFILE  \n\
echo ------------------------------------------------------  \n\
echo ' '  \n\
echo ' '  \n\
  \n\
 cd \${WORKDIR};  \n\
  \n\
# COMMAND to EXECUTE:  \n\
 \${LAUNCH} \${EXEFILE} \${RUNFILES}  \n\
  \n\
exit"


# Contents of the script to run/submit all fits to the cluster
RFfile="#!/bin/sh  \n\
  \n\
    for dir0 in \`ls -d Fits_*\`  \n\
    do  \n\
  \n\
	cd \$dir0  \n\
# SM fits  \n\
	if [ -d SM_fits ]  \n\
	then  \n\
	    cd SM_fits  \n\
  \n\
# Full fit  \n\
	    if [ -d Full_Fit ]  \n\
          then  \n\
          cd Full_Fit  \n\
            if [ ! -f MonteCarlo_plots.pdf ]  \n\
            then  \n\
#-------------------------------------------------------------  \n\
#echo \"submitting job in \$PWD\"  \n\
#qsub -w \`pwd\` -vARGS=\"model.conf ../MonteCarlo-HIGH.conf\" -N SM_FullFit /storage/local/home/theorm3/silvestrini/storage/HEPfit/Analysis_EW/submit_analysis_v1.sh  \n\
                    rm -rf analysis core.*  \n\
                    #ln -s ~/storage/HEPfit/Analysis/dist/Debug_MPI/GNU-Linux/analysis .  \n\
                    sed -i -e \"s#NAME#SM_FullFit#\" -e \"4s#/##g\" submit_analysis_IB.sh  \n\
                    sed -i -e \"s#../../MonteCarlo-HIGH.conf#../MonteCarlo-HIGH.conf#\" submit_analysis_IB.sh  \n\
                    qsub submit_analysis_IB.sh  \n\
                else  \n\
                    rm -f *.o* analysis core.*  \n\
#-------------------------------------------------------------  \n\
            fi  \n\
          cd ..  \n\
          fi  \n\
  \n\
# All the rest: Indirect determinations, parametric uncertainties, predictions, MW_vs_X plots  \n\
  \n\
          for dir1 in Ind_determ Param_Unc Pred_Obs MW_mt_fits MW_sin2Eff_fits  \n\
          do  \n\
  \n\
	    if [ -d \${dir1} ]  \n\
          then  \n\
          cd \${dir1}  \n\
	    for dir2 in \`ls -d */\`  \n\
	    do  \n\
		cd \${dir2}  \n\
            if [ ! -f MonteCarlo_plots.pdf ]  \n\
            then  \n\
#-------------------------------------------------------------  \n\
#echo \"submitting job in \$PWD\"  \n\
#qsub -w \`pwd\` -vARGS=\"model.conf ../../MonteCarlo-HIGH.conf\" -N \${dir0}_\${dir1}_\${dir2} /storage/local/home/theorm3/silvestrini/storage/HEPfit/Analysis_EW/submit_analysis_v1.sh  \n\
                    rm -rf analysis core.*  \n\
                    #ln -s ~/storage/HEPfit/Analysis/dist/Debug_MPI/GNU-Linux/analysis .  \n\
                    sed -i -e \"s#NAME#\${dir0}_\${dir1}_\${dir2}#\" -e \"4s#/##g\" submit_analysis_IB.sh  \n\
                    qsub submit_analysis_IB.sh  \n\
                else  \n\
                    rm -f *.o* analysis core.*  \n\
#-------------------------------------------------------------  \n\
            fi  \n\
		cd ..  \n\
	    done  \n\
	    cd ..  \n\
          fi  \n\
  \n\
          done  \n\
          cd ..  \n\
  \n\
	fi  \n\
  \n\
	cd ..  \n\
  \n\
    done  \n\
#----------"

#------------------------------------------------------------------------------

# Before star creating the config files for the fits place the script to submit them to the cluster here (root folder)
# Remove any previous version. (The script is universal anyway so it should be compatible with whatever done in any other previous version)

rm run_SM_fits.sh

(echo -e "$RFfile")>run_SM_fits.sh

#------------------------------------------------------------------------------

# The script

# 0. Initial settings

# 0.1 Create folder and put there all the necessary files

for dir0 in Fits_Current
do

# Get the type of fit to perform from the value in dir
# Remove the "Fits_" part and the "/" at the end
FitType=${dir0/#"Fits_"}
FitType=$(echo $FitType | cut -d'/' -f-1)
echo ' Generating fits for '$FitType
echo ' '

cd ${dir0}

mkdir SM_fits

cp  model.conf SM_fits/.
cp  ../FlavourFixed.conf SM_fits/.
#     (In some cases there are several ObservablesEW*.conf files)
cp  ObservablesEW* SM_fits/.


# 0.2 Go inside the SM_fits folder and adjust the global files
cd SM_fits

# Select the file with or without LFU assumptions
if [ $LFUobs != 'Yes' ] && [ $LFUobs != 'yes' ]
then
      sed -i "" 's/IncludeFile ObservablesEW_Current_LFU.conf/#IncludeFile ObservablesEW_Current_LFU.conf/g' ObservablesEW.conf
      sed -i "" 's/#IncludeFile ObservablesEW_Current_noLFU.conf/IncludeFile ObservablesEW_Current_noLFU.conf/g' ObservablesEW.conf
      
      sed -i "" 's/IncludeFile ObservablesEW_Current_SM_LFU.conf/#IncludeFile ObservablesEW_Current_SM_LFU.conf/g' ObservablesEW.conf
      sed -i "" 's/#IncludeFile ObservablesEW_Current_SM_noLFU.conf/IncludeFile ObservablesEW_Current_SM_noLFU.conf/g' ObservablesEW.conf
else
      sed -i "" 's/#IncludeFile ObservablesEW_Current_LFU.conf/IncludeFile ObservablesEW_Current_LFU.conf/g' ObservablesEW.conf
      sed -i "" 's/IncludeFile ObservablesEW_Current_noLFU.conf/#IncludeFile ObservablesEW_Current_noLFU.conf/g' ObservablesEW.conf
      
      sed -i "" 's/#IncludeFile ObservablesEW_Current_SM_LFU.conf/IncludeFile ObservablesEW_Current_SM_LFU.conf/g' ObservablesEW.conf
      sed -i "" 's/IncludeFile ObservablesEW_Current_SM_noLFU.conf/#IncludeFile ObservablesEW_Current_SM_noLFU.conf/g' ObservablesEW.conf
fi

# IMPORTANT: For the SM fits it is important to remove from the EW observables any of the "delta" or "AuxObsNP" observables, which are
# only defined in NP models
for obsEW in `ls ObservablesEW*`
do

# Comment out delta and AuxObsNP obs.
      sed -i "" 's/\(Observable.*deltag.*\)/##\1 /g' $obsEW
      sed -i "" 's/\(Observable.*deltaV.*\)/##\1 /g' $obsEW
      sed -i "" 's/\(Observable.*AuxObsNP.*\)/##\1 /g' $obsEW

# Comment out correlated gaussian obs for D0 and the correlations
      sed -i "" 's/CorrelatedGaussianObservables ZqqD0 4/##CorrelatedGaussianObservables ZqqD0 4/g' $obsEW
      sed -i "" 's/1.000  0.470 -0.201 -0.217/##1.000  0.470 -0.201 -0.217/g' $obsEW
      sed -i "" 's/0.470  1.000  0.606 -0.925/##0.470  1.000  0.606 -0.925/g' $obsEW
      sed -i "" 's/-0.201  0.606  1.000 -0.813/##-0.201  0.606  1.000 -0.813/g' $obsEW
      sed -i "" 's/-0.217 -0.925 -0.813  1.000/##-0.217 -0.925 -0.813  1.000/g' $obsEW

done

# Create the global Monte Carlo config file

#     Put the MonteCarlo file here (use the " " to preserve the blanks)
(echo -e "$MCfile")>MonteCarlo-HIGH.conf
sed -i "" 's/PrerunMaxIter              1000000/PrerunMaxIter              10000000/' MonteCarlo-HIGH.conf
sed -i "" 's/Iterations                 2000000/Iterations                 5000000/' MonteCarlo-HIGH.conf

#     Put the basic cluster submission file here and adjust configuration
(echo -e "$CSfile")>submit_analysis_IB.sh
sed -i "" 's/EXEFILE=\/storage\/local\/home\/theorm3\/deblasmateo\/storage\/Software\/HEPfit\/Master\/HEPfit\/Analysis\/dist\/Debug_MPI\/GNU-Linux\/analysis/EXEFILE='$HEPfitPATH'/' submit_analysis_IB.sh

# Adjust the model file

# Add the model name at the beginning (keep the line structure for sed to work properly in non-GNU sed)
sed -i "" '1 i\
StandardModel
' model.conf

# Apply settings about theory errors to the global model file

# Add the SM theoretical uncertainties specified at the beginning of this file
if [ $IntUPrior == 'Gaussian' ]
then
#     Gaussian priors for intrinsic uncertainties
            sed -i "" 's/ModelParameter  delMw       0.          0.          0./ModelParameter  delMw       0.          '$MWthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delSin2th_l 0.          0.          0./ModelParameter  delSin2th_l 0.          '$s2lEffthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delSin2th_q 0.          0.          0./ModelParameter  delSin2th_q 0.          '$s2qEffthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delSin2th_b 0.          0.          0./ModelParameter  delSin2th_b 0.          '$s2bEffthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delGammaZ   0.          0.          0./ModelParameter  delGammaZ   0.          '$GammaZthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delsigma0H  0.          0.          0./ModelParameter  delsigma0H  0.          '$SigmaHthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delR0l      0.          0.          0./ModelParameter  delR0l      0.          '$RlthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delR0c      0.          0.          0./ModelParameter  delR0c      0.          '$RcthU'          0./' model.conf
            sed -i "" 's/ModelParameter  delR0b      0.          0.          0./ModelParameter  delR0b      0.          '$RbthU'          0./' model.conf

elif [ $IntUPrior == 'Flat' ]
then
#     Flat priors for intrinsic uncertainties
            sed -i "" 's/ModelParameter  delMw       0.          0.          0./ModelParameter  delMw       0.          0.          '$MWthU'/' model.conf
            sed -i "" 's/ModelParameter  delSin2th_l 0.          0.          0./ModelParameter  delSin2th_l 0.          0.          '$s2lEffthU'/' model.conf
            sed -i "" 's/ModelParameter  delSin2th_q 0.          0.          0./ModelParameter  delSin2th_q 0.          0.          '$s2qEffthU'/' model.conf
            sed -i "" 's/ModelParameter  delSin2th_b 0.          0.          0./ModelParameter  delSin2th_b 0.          0.          '$s2bEffthU'/' model.conf
            sed -i "" 's/ModelParameter  delGammaZ   0.          0.          0./ModelParameter  delGammaZ   0.          0.          '$GammaZthU'/' model.conf
            sed -i "" 's/ModelParameter  delsigma0H  0.          0.          0./ModelParameter  delsigma0H  0.          0.          '$SigmaHthU'/' model.conf
            sed -i "" 's/ModelParameter  delR0l      0.          0.          0./ModelParameter  delR0l      0.          0.          '$RlthU'/' model.conf
            sed -i "" 's/ModelParameter  delR0c      0.          0.          0./ModelParameter  delR0c      0.          0.          '$RcthU'/' model.conf
            sed -i "" 's/ModelParameter  delR0b      0.          0.          0./ModelParameter  delR0b      0.          0.          '$RbthU'/' model.conf
fi

# Remove unnecessary lines
sed -i "" '/#MZExp ModelParameter  Mz          91.1875     0.0021    0./d' model.conf

# Add the observables for the MW-mt, MW-sin2Eff plots at the end of the model.conf file
(echo "#")>>model.conf
(echo "# Extra 2D Observables")>>model.conf
(echo "Observable2D  MW_vs_mt Mw Mw 1. -1. noMCMC noweight mtop mtop 1. -1.")>>model.conf
(echo "Observable2D  MW_vs_sin2Eff Mw Mw 1. -1. noMCMC noweight sin2thetaEff sin2thetaEff 1. -1.")>>model.conf

(echo "#")>>model.conf
(echo "# Extra Correlated Observables")>>model.conf

(echo "#")>>model.conf
(echo "CorrelatedObservables MWmtcorr 2")>>model.conf
(echo "Observable  Mw_corr1 Mw M_{W}  0.  0.  noMCMC noweight")>>model.conf
(echo "Observable  mtop_corr1 mtop m_{t}  0.  0.  noMCMC noweight")>>model.conf

(echo "#")>>model.conf
(echo "CorrelatedObservables MWsin2Effcorr 2")>>model.conf
(echo "Observable  Mw_corr2 Mw M_{W}  0.  0.  noMCMC noweight")>>model.conf
(echo "Observable  sin2thetaEff_corr2 sin2thetaEff sin^{2}#theta_{eff}  0.  0.  noMCMC noweight")>>model.conf

#------------------------------------------------------------------------------

# 1. Global fit

if [ $GenFullfit == 'Yes' ] || [ $GenFullfit == 'yes' ]
then
      mkdir Full_Fit

# Copy the model and observables files
      cp  model.conf Full_Fit/.
      cp  FlavourFixed.conf Full_Fit/.
#     (In some cases there are several ObservablesEW*.conf files)
      cp  ObservablesEW* Full_Fit/.

#     Copy the cluster submission script and make adjustments if needed
      cp  submit_analysis_IB.sh Full_Fit/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' Full_Fit/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' Full_Fit/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' Full_Fit/submit_analysis_IB.sh
      sed -i "" 's/model.conf ..\/..\/MonteCarlo-HIGH.conf/model.conf ..\/'$MCCONFLG'/' Full_Fit/submit_analysis_IB.sh
fi

#------------------------------------------------------------------------------

# 2. Indirect determinations of SM input parameters

if [ $GenIndDet == 'Yes' ] || [ $GenIndDet == 'yes' ]
then

mkdir Ind_determ

cd Ind_determ

mkdir alphaS

mkdir Dalpha5h

mkdir mtop

mkdir mHiggs

mkdir MZ

for dir in `ls -d */`
do

# Copy the model and observables files
     cp  ../model.conf ${dir}/.
     cp  ../FlavourFixed.conf ${dir}/.
#     (In some cases there are several ObservablesEW*.conf files)
     cp  ../ObservablesEW* ${dir}/.

#     Copy the cluster submission script and make adjustments if needed
      cp  ../submit_analysis_IB.sh ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/MonteCarlo-HIGH.conf/'$MCCONFLG'/' ${dir}/submit_analysis_IB.sh

done

# Adjust the model and config files for each case:
# For each fit:   model.conf: Assign large flat prior to the corresponding SM model parameter
#                 ObservablesEW.conf: For MZ, switch the observable corresponding to MZ to noMCMC noweight

# alpha_S (Don't go beyond +-0.02 or otherwise it runs into numerical problems with the running of alpha_S)
sed -i "" 's/ModelParameter  AlsMz .*/ModelParameter  AlsMz       '$alSexp'      0.          0.02/' alphaS/model.conf

# Dalpha_5H
sed -i "" 's/ModelParameter  dAle5Mz .*/ModelParameter  dAle5Mz     '$dal5hexp'     0.          0.005/' Dalpha5h/model.conf

# mt
sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtexp'      0.          20./' mtop/model.conf

# mH
sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         505.      0.          494./' mHiggs/model.conf

# MZ
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZexp'     0.          0.2/' MZ/model.conf
#sed -i "" 's/Observable  Mz           Mz           M_{Z} 1. -1. MCMC weight/Observable  Mz           Mz           M_{Z} 1. -1. noMCMC noweight/g' MZ/$ObsEWFU
#sed -i "" 's/Observable  Mz             Mz             M_{Z}            1. -1. MCMC weight/Observable  Mz             Mz             M_{Z}            1. -1. noMCMC noweight/g' MZ/$ObsEWnoFU

#sed -i "" 's/Observable  Mz_C           Mz           M_{Z} 1. -1. MCMC weight/Observable  Mz_C           Mz           M_{Z} 1. -1. noMCMC noweight/g' MZ/$ObsEWFU
#sed -i "" 's/Observable  Mz_C             Mz             M_{Z}            1. -1. MCMC weight/Observable  Mz_C             Mz             M_{Z}            1. -1. noMCMC noweight/g' MZ/$ObsEWnoFU

for obsEW in `ls MZ/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*Mz.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# Go back to main folder (Current_Future/Fits_Current/SM_fits)
cd ..

fi

#------------------------------------------------------------------------------

# 3. Parametric uncertainties

if [ $GenParUnc == 'Yes' ] || [ $GenParUnc == 'yes' ]
then

mkdir Param_Unc

cd Param_Unc

mkdir Total

mkdir alphaS

mkdir Dalpha5h

mkdir mtop

mkdir mHiggs

mkdir MZ

# Copy the model file
cp  ../model.conf .

# and make sure to remove SM theory uncertainties (if these were included in the calculation)
# so only the effects of the parametric errors are shown
sed -i "" 's/ModelParameter  delMw .*/ModelParameter  delMw       0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delSin2th_l .*/ModelParameter  delSin2th_l 0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delSin2th_q .*/ModelParameter  delSin2th_q 0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delSin2th_b .*/ModelParameter  delSin2th_b 0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delGammaZ .*/ModelParameter  delGammaZ   0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delsigma0H .*/ModelParameter  delsigma0H   0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delR0l .*/ModelParameter  delR0l   0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delR0c .*/ModelParameter  delR0c   0.          0.          0./' model.conf
sed -i "" 's/ModelParameter  delR0b .*/ModelParameter  delR0b      0.          0.          0./' model.conf

for dir in `ls -d */`
do

# Copy the model and observables files
# (Use the modified model.con without theory uncertainties in /Param_Unc
     cp  model.conf ${dir}/.
     cp  ../FlavourFixed.conf ${dir}/.
#     (In some cases there are several ObservablesEW*.conf files)
     cp  ../ObservablesEW* ${dir}/.

#     Copy the cluster submission script and make adjustments if needed
      cp  ../submit_analysis_IB.sh ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/MonteCarlo-HIGH.conf/'$MCCONFLG'/' ${dir}/submit_analysis_IB.sh

done

# Remove the modified model file without SM theory unc (not needed anymore)
rm model.conf

# Adjust the model and config files for each case:
# For Total fit:  model.conf: Use the experimental error for MZ (Gaussian)
#                 ObservablesEW.conf: Switch all observables to noMCMC noweight

# For each fit:   model.conf: Fix all the SM parameter to the best fit value (mode) of global fit,
#                             except for the relevant parameter which is given exactly the exp. value and error (Gaussian)
#                 ObservablesEW.conf: Switch all observables to noMCMC noweight

# Total par. unc
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZexp'     '$MZexpErr'          0./' Total/model.conf
#sed -i "" 's/MCMC weight/noMCMC noweight/g' Total/$ObsEWFU
#sed -i "" 's/MCMC weight/noMCMC noweight/g' Total/$ObsEWnoFU
for obsEW in `ls Total/ObservablesEW*`
do
      sed -i "" 's/MCMC weight/noMCMC noweight/g' $obsEW
done

# alpha_S par. unc
sed -i "" 's/ModelParameter  AlsMz .*/ModelParameter  AlsMz       '$alSexp'      '$alSexpErr'          0./' alphaS/model.conf
sed -i "" 's/ModelParameter  dAle5Mz .*/ModelParameter  dAle5Mz     '$dal5hSM'     0.          0./' alphaS/model.conf
sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtSM'      0.          0./' alphaS/model.conf
sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         '$mHSM'      0.          0./' alphaS/model.conf
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZSM'     0.          0./' alphaS/model.conf
#sed -i "" 's/MCMC weight/noMCMC noweight/g' alphaS/$ObsEWFU
#sed -i "" 's/MCMC weight/noMCMC noweight/g' alphaS/$ObsEWnoFU
for obsEW in `ls alphaS/ObservablesEW*`
do
      sed -i "" 's/MCMC weight/noMCMC noweight/g' $obsEW
done

# DAlpha_5H par. unc
sed -i "" 's/ModelParameter  AlsMz .*/ModelParameter  AlsMz       '$alSSM'      0.          0./' Dalpha5h/model.conf
sed -i "" 's/ModelParameter  dAle5Mz .*/ModelParameter  dAle5Mz     '$dal5hexp'     '$dal5hexpErr'          0./' Dalpha5h/model.conf
sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtSM'      0.          0./' Dalpha5h/model.conf
sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         '$mHSM'      0.          0./' Dalpha5h/model.conf
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZSM'     0.          0./' Dalpha5h/model.conf
#sed -i "" 's/MCMC weight/noMCMC noweight/g' Dalpha5h/$ObsEWFU
#sed -i "" 's/MCMC weight/noMCMC noweight/g' Dalpha5h/$ObsEWnoFU
for obsEW in `ls Dalpha5h/ObservablesEW*`
do
      sed -i "" 's/MCMC weight/noMCMC noweight/g' $obsEW
done

# mtop par. unc
sed -i "" 's/ModelParameter  AlsMz .*/ModelParameter  AlsMz       '$alSSM'      0.          0./' mtop/model.conf
sed -i "" 's/ModelParameter  dAle5Mz .*/ModelParameter  dAle5Mz     '$dal5hSM'     0.          0./' mtop/model.conf
sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtexp'      '$mtexpErr'          0./' mtop/model.conf
sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         '$mHSM'      0.          0./' mtop/model.conf
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZSM'     0.          0./' mtop/model.conf
#sed -i "" 's/MCMC weight/noMCMC noweight/g' mtop/$ObsEWFU
#sed -i "" 's/MCMC weight/noMCMC noweight/g' mtop/$ObsEWnoFU
for obsEW in `ls mtop/ObservablesEW*`
do
      sed -i "" 's/MCMC weight/noMCMC noweight/g' $obsEW
done

# mH par. unc
sed -i "" 's/ModelParameter  AlsMz .*/ModelParameter  AlsMz       '$alSSM'      0.          0./' mHiggs/model.conf
sed -i "" 's/ModelParameter  dAle5Mz .*/ModelParameter  dAle5Mz     '$dal5hSM'     0.          0./' mHiggs/model.conf
sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtSM'      0.          0./' mHiggs/model.conf
sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         '$mHexp'      '$mHexpErr'          0./' mHiggs/model.conf
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZSM'     0.          0./' mHiggs/model.conf
#sed -i "" 's/MCMC weight/noMCMC noweight/g' mHiggs/$ObsEWFU
#sed -i "" 's/MCMC weight/noMCMC noweight/g' mHiggs/$ObsEWnoFU
for obsEW in `ls mHiggs/ObservablesEW*`
do
      sed -i "" 's/MCMC weight/noMCMC noweight/g' $obsEW
done

# MZ par. unc
sed -i "" 's/ModelParameter  AlsMz .*/ModelParameter  AlsMz       '$alSSM'      0.          0./' MZ/model.conf
sed -i "" 's/ModelParameter  dAle5Mz .*/ModelParameter  dAle5Mz     '$dal5hSM'     0.          0./' MZ/model.conf
sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtSM'      0.          0./' MZ/model.conf
sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         '$mHSM'      0.          0./' MZ/model.conf
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZexp'     '$MZexpErr'          0./' MZ/model.conf
#sed -i "" 's/MCMC weight/noMCMC noweight/g' MZ/$ObsEWFU
#sed -i "" 's/MCMC weight/noMCMC noweight/g' MZ/$ObsEWnoFU
for obsEW in `ls MZ/ObservablesEW*`
do
      sed -i "" 's/MCMC weight/noMCMC noweight/g' $obsEW
done

# Go back to main folder (Current_Future/Fits_Current/SM_fits)
cd ..

fi

#------------------------------------------------------------------------------

# 4. Predictions Obs

if [ $GenPredObs == 'Yes' ] || [ $GenPredObs == 'yes' ]
then

mkdir Pred_Obs

cd Pred_Obs

if [ $LFUobs == 'Yes' ] || [ $LFUobs == 'yes' ]
then

# LFU case
# --------

mkdir Mw

mkdir GammaW

mkdir PtauPol

mkdir sin2EfflLEPHC

mkdir sin2EfflHC

mkdir sin2EfflLEP

mkdir As

mkdir Ruc

mkdir RWc

mkdir Zpole1

mkdir Zpole2

# mkdir BRWhad  # Removed from standard fits. Not independent from BRWlept

mkdir BRWlept

for dir in `ls -d */`
do

# Copy the model and observables files
     cp  ../model.conf ${dir}/.
     cp  ../FlavourFixed.conf ${dir}/.
#     (In some cases there are several ObservablesEW*.conf files)
     cp  ../ObservablesEW* ${dir}/.

#     Copy the cluster submission script and make adjustments if needed
      cp  ../submit_analysis_IB.sh ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/MonteCarlo-HIGH.conf/'$MCCONFLG'/' ${dir}/submit_analysis_IB.sh

done

# Adjust the model and config files for each case:

# For each fit:   model.conf: Only for Zpole1, change the error of MZ from flat to the experimental one (Gaussian)
#                 ObservablesEW.conf: Switch the corresponding observable (set of correllated observables) to noMCMC noweight

# NOTE that several independent measurements of the same observable can be enteres.
# Match the theory name of the Obs and replace globally, i.e.
# Match Observable[Whatever]NameObs[Whatever]MCMC weight and replace it by Observable[Whatever]NameObs[Whatever]noMCMC noweight
# Use parenthesis \(  \) to keep the part of the matched expr that won't be replaced (its on the \1

# No MW
#sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' Mw/$ObsEWFU
for obsEW in `ls Mw/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No GammaW
#sed -i "" 's/\(Observable.*GammaW.*\) MCMC weight/\1 noMCMC noweight/g' GammaW/$ObsEWFU
for obsEW in `ls GammaW/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*GammaW.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No PtauPol
#sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' PtauPol/$ObsEWFU
for obsEW in `ls PtauPol/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No sinethetaEffl (All)
#sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' sin2EfflLEPHC/$ObsEWFU
for obsEW in `ls sin2EfflLEPHC/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# The following 2 are special because I want to remove only some measurements for the same observable so I need to match the
# name of the exp meas too

# No sinethetaEffl (Hadron coll)
#sed -i "" 's/\(Observable  sW2eff_.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' sin2EfflHC/$ObsEWFU
for obsEW in `ls sin2EfflHC/ObservablesEW*`
do
      sed -i "" 's/\(Observable  sW2eff_.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No sinethetaEffl (LEP)
#sed -i "" 's/\(Observable  sin2thetaEff sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' sin2EfflLEP/$ObsEWFU
for obsEW in `ls sin2EfflLEP/ObservablesEW*`
do
      sed -i "" 's/\(Observable  sin2thetaEff sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No As
for obsEW in `ls As/ObservablesEW*`
do
      sed -i "" 's/\(Observable  As       Astrange.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Ruc
for obsEW in `ls Ruc/ObservablesEW*`
do
      sed -i "" 's/\(Observable  Ruc      Ruc.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No RWc
for obsEW in `ls RWc/ObservablesEW*`
do
      sed -i "" 's/\(Observable  RWc      RWc.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Zpole 1
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZexp'     '$MZexpErr'          0./' Zpole1/model.conf
#sed -i "" 's/CorrelatedGaussianObservables Zpole1 5/CorrelatedObservables Zpole1 5/g' Zpole1/$ObsEWFU
#sed -i "" 's/\(Observable.*Mz.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWFU
#sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWFU
#sed -i "" 's/\(Observable.*sigmaHadron.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWFU
#sed -i "" 's/\(Observable.*Rlepton.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWFU
#sed -i "" 's/\(Observable.*AFBlepton.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWFU
## WARNING: The next 5 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/1.000   -0.023  -0.045   0.033   0.055/# /' Zpole1/$ObsEWFU
#sed -i "" '/-0.023   1.000  -0.297   0.004   0.003/d' Zpole1/$ObsEWFU
#sed -i "" '/-0.045  -0.297   1.000   0.183   0.006/d' Zpole1/$ObsEWFU
#sed -i "" '/0.033    0.004   0.183   1.000  -0.056/d' Zpole1/$ObsEWFU
#sed -i "" '/0.055    0.003   0.006  -0.056   1.000/d' Zpole1/$ObsEWFU
for obsEW in `ls Zpole1/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables Zpole1 5/CorrelatedObservables Zpole1 5/g' $obsEW
      sed -i "" 's/\(Observable.*Mz.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*sigmaHadron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rlepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBlepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 5 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/1.0000   -0.0228  -0.0521   0.0332   0.0549/# /' $obsEW
      sed -i "" '/-0.0228   1.0000  -0.3248   0.0037   0.0033/d' $obsEW
      sed -i "" '/-0.0521  -0.3248   1.0000   0.1960   0.0069/d' $obsEW
      sed -i "" '/0.0332    0.0037   0.1960   1.0000  -0.0560/d' $obsEW
      sed -i "" '/0.0549    0.0033   0.0069  -0.0560   1.0000/d' $obsEW

#     Treat also observables with no LFU
      sed -i "" 's/\(Observable.*Relectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rmuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBmuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Zpole 2
#sed -i "" 's/CorrelatedGaussianObservables Zpole2 7/CorrelatedObservables Zpole2 7/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*Rbottom.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*Rcharm.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*AFBbottom.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*AFBcharm.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*Abottom.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
#sed -i "" 's/\(Observable.*Acharm.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWFU
## WARNING: The next 7 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/1.00   0.00   0.00   0.00   0.00   0.09   0.05/# /' Zpole2/$ObsEWFU
#sed -i "" '/0.00   1.00  -0.18  -0.10   0.07  -0.08   0.04/d' Zpole2/$ObsEWFU
#sed -i "" '/0.00  -0.18   1.00   0.04  -0.06   0.04  -0.06/d' Zpole2/$ObsEWFU
#sed -i "" '/0.00  -0.10   0.04   1.00   0.15   0.06   0.01/d' Zpole2/$ObsEWFU
#sed -i "" '/0.00   0.07  -0.06   0.15   1.00  -0.02   0.04/d' Zpole2/$ObsEWFU
#sed -i "" '/0.09  -0.08   0.04   0.06  -0.02   1.00   0.11/d' Zpole2/$ObsEWFU
#sed -i "" '/0.05   0.04  -0.06   0.01   0.04   0.11   1.00/d' Zpole2/$ObsEWFU
for obsEW in `ls Zpole2/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables Zpole2 7/CorrelatedObservables Zpole2 7/g' $obsEW
      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rbottom.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rcharm.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBbottom.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBcharm.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Abottom.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Acharm.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 7 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/1.00   0.00   0.00   0.00   0.00   0.09   0.05/# /' $obsEW
      sed -i "" '/0.00   1.00  -0.18  -0.10   0.07  -0.08   0.04/d' $obsEW
      sed -i "" '/0.00  -0.18   1.00   0.04  -0.06   0.04  -0.06/d' $obsEW
      sed -i "" '/0.00  -0.10   0.04   1.00   0.15   0.06   0.01/d' $obsEW
      sed -i "" '/0.00   0.07  -0.06   0.15   1.00  -0.02   0.04/d' $obsEW
      sed -i "" '/0.09  -0.08   0.04   0.06  -0.02   1.00   0.11/d' $obsEW
      sed -i "" '/0.05   0.04  -0.06   0.01   0.04   0.11   1.00/d' $obsEW

#     Treat also observables with no LFU
      sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No BRWhad  # Removed from standard fits. Not independent from BRWlept
#sed -i "" 's/\(Observable.*BrWhadrons.*\) MCMC weight/\1 noMCMC noweight/g' BRWhad/$ObsEWFU
#for obsEW in `ls BRWhad/ObservablesEW*`
#do
#      sed -i "" 's/\(Observable.*BrWhadrons.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#done

# No BRWlept
#sed -i "" 's/\(Observable.*BrWlepton.*\) MCMC weight/\1 noMCMC noweight/g' BRWlept/$ObsEWFU
for obsEW in `ls BRWlept/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*BrWlepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW

#     Treat also observables with no LFU
      sed -i "" 's/\(Observable.*BrWelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*BrWmuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*BrWtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

fi

#------------------------------------------------------------------------------

if [ $LFUobs != 'Yes' ] && [ $LFUobs != 'yes' ]
then

# Non-LFU case
# ------------

mkdir Mw

mkdir GammaW

mkdir sin2EfflLEPHC

mkdir sin2EfflHC

mkdir sin2EfflLEP

mkdir As

mkdir Ruc

mkdir RWc

mkdir Zpole1

mkdir Zpole2

mkdir Zpole3

mkdir Zpole4

# mkdir BRWhad  # Removed from standard fits. Not independent from BRWlept

mkdir BRWlept

# Ignore the LFU tests in the prediction section. The predictions for those are fixed in the SM to ~1 (up to lepton mass effects, which are not floated anyway)

for dir in `ls -d */`
do

# Copy the model and observables files
     cp  ../model.conf ${dir}/.
     cp  ../FlavourFixed.conf ${dir}/.
#     (In some cases there are several ObservablesEW*.conf files)
     cp  ../ObservablesEW* ${dir}/.

#     Copy the cluster submission script and make adjustments if needed
      cp  ../submit_analysis_IB.sh ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/MonteCarlo-HIGH.conf/'$MCCONFLG'/' ${dir}/submit_analysis_IB.sh

done

# Adjust the model and config files for each case:

# For each fit:   model.conf: Only for Zpole1, change the error of MZ from flat to the experimental one (Gaussian)
#                 ObservablesEW.conf: Switch the corresponding observable (set of correllated observables) to noMCMC noweight

# NOTE that several independent measurements of the same observable can be enteres.
# Match the theory name of the Obs and replace globally, i.e.
# Match Observable[Whatever]NameObs[Whatever]MCMC weight and replace it by Observable[Whatever]NameObs[Whatever]noMCMC noweight
# Use parenthesis \(  \) to keep the part of the matched expr that won't be replaced (its on the \1

# No MW
#sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' Mw/$ObsEWnoFU
for obsEW in `ls Mw/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No GammaW
#sed -i "" 's/\(Observable.*GammaW.*\) MCMC weight/\1 noMCMC noweight/g' GammaW/$ObsEWnoFU
for obsEW in `ls GammaW/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*GammaW.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No sinethetaEffl (All)
#sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' sin2EfflLEPHC/$ObsEWnoFU
for obsEW in `ls sin2EfflLEPHC/ObservablesEW*`
do
      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# The following 2 are special because I want to remove only some measurements for the same observable so I need to match the
# name of the exp meas too

# No sinethetaEffl (Hadron coll)
#sed -i "" 's/\(Observable  sW2eff_.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' sin2EfflHC/$ObsEWnoFU
for obsEW in `ls sin2EfflHC/ObservablesEW*`
do
      sed -i "" 's/\(Observable  sW2eff_.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No sinethetaEffl (LEP)
#sed -i "" 's/\(Observable  sin2thetaEff sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' sin2EfflLEP/$ObsEWnoFU
for obsEW in `ls sin2EfflLEP/ObservablesEW*`
do
      sed -i "" 's/\(Observable  sin2thetaEff sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No As
for obsEW in `ls As/ObservablesEW*`
do
      sed -i "" 's/\(Observable  As       Astrange.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Ruc
for obsEW in `ls Ruc/ObservablesEW*`
do
      sed -i "" 's/\(Observable  Ruc      Ruc.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No RWc
for obsEW in `ls RWc/ObservablesEW*`
do
      sed -i "" 's/\(Observable  RWc      RWc.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Zpole 1
sed -i "" 's/ModelParameter  Mz .*/ModelParameter  Mz          '$MZexp'     '$MZexpErr'          0./' Zpole1/model.conf
#sed -i "" 's/CorrelatedGaussianObservables Zpole1 9/CorrelatedObservables Zpole1 9/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Mz.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*sigmaHadron.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Relectron.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Rmuon.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Rtau.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AFBelectron.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AFBmuon.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AFBtau.*\) MCMC weight/\1 noMCMC noweight/g' Zpole1/$ObsEWnoFU
## WARNING: The next 9 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/ 1.000 -0.024 -0.044  0.078  0.000  0.002 -0.014  0.046  0.035/# /' Zpole1/$ObsEWnoFU
#sed -i "" '/-0.024  1.000 -0.297 -0.011  0.008  0.006  0.007  0.002  0.001/d' Zpole1/$ObsEWnoFU
#sed -i "" '/-0.044 -0.297  1.000  0.105  0.131  0.092  0.001  0.003  0.002/d' Zpole1/$ObsEWnoFU
#sed -i "" '/ 0.078 -0.011  0.105  1.000  0.069  0.046 -0.371  0.020  0.013/d' Zpole1/$ObsEWnoFU
#sed -i "" '/ 0.000  0.008  0.131  0.069  1.000  0.069  0.001  0.012 -0.003/d' Zpole1/$ObsEWnoFU
#sed -i "" '/ 0.002  0.006  0.092  0.046  0.069  1.000  0.003  0.001  0.009/d' Zpole1/$ObsEWnoFU
#sed -i "" '/-0.014  0.007  0.001 -0.371  0.001  0.003  1.000 -0.024 -0.020/d' Zpole1/$ObsEWnoFU
#sed -i "" '/ 0.046  0.002  0.003  0.020  0.012  0.001 -0.024  1.000  0.046/d' Zpole1/$ObsEWnoFU
#sed -i "" '/ 0.035  0.001  0.002  0.013 -0.003  0.009 -0.020  0.046  1.000/d' Zpole1/$ObsEWnoFU
for obsEW in `ls Zpole1/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables Zpole1 9/CorrelatedObservables Zpole1 9/g' $obsEW
      sed -i "" 's/\(Observable.*Mz.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*sigmaHadron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Relectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rmuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBmuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 9 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/ 1.0000 -0.0238 -0.0507  0.0783  0.0000  0.0019 -0.0139  0.0459  0.0346/# /' $obsEW
      sed -i "" '/-0.0238  1.0000 -0.3249 -0.0110  0.0079  0.0059  0.0071  0.0020  0.0013/d' $obsEW
      sed -i "" '/-0.0507 -0.3249  1.0000  0.1138  0.1391  0.0987  0.0015  0.0035  0.0018/d' $obsEW
      sed -i "" '/ 0.0783 -0.0110  0.1138  1.0000  0.0694  0.0464 -0.3704  0.0197  0.0132/d' $obsEW
      sed -i "" '/ 0.0000  0.0079  0.1391  0.0694  1.0000  0.0696  0.0013  0.0121 -0.0030/d' $obsEW
      sed -i "" '/ 0.0019  0.0059  0.0987  0.0464  0.0696  1.0000  0.0029  0.0012  0.0093/d' $obsEW
      sed -i "" '/-0.0139  0.0071  0.0015 -0.3704  0.0013  0.0029  1.0000 -0.0242 -0.0202/d' $obsEW
      sed -i "" '/ 0.0459  0.0020  0.0035  0.0197  0.0121  0.0012 -0.0242  1.0000  0.0464/d' $obsEW
      sed -i "" '/ 0.0346  0.0013  0.0018  0.0132 -0.0030  0.0093 -0.0202  0.0464  1.0000/d' $obsEW

#     Treat also observables with LFU
      sed -i "" 's/\(Observable.*Rlepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBlepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Zpole 2: Do a more precise match, to distinguish from the measurements in Zpole3
#sed -i "" 's/CorrelatedGaussianObservables Zpole2 3/CorrelatedObservables Zpole2 3/g' Zpole2/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Aelectron   Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Amuon       Amuon.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Atau        Atau.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWnoFU
# _C version
#sed -i "" 's/\(Observable.*Aelectron_C   Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Amuon_C       Amuon.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Atau_C        Atau.*\) MCMC weight/\1 noMCMC noweight/g' Zpole2/$ObsEWnoFU
# WARNING: The next 3 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/1.000  0.038  0.033/# /' Zpole2/$ObsEWnoFU
#sed -i "" '/0.038  1.000  0.007/d' Zpole2/$ObsEWnoFU
#sed -i "" '/0.033  0.007  1.000/d' Zpole2/$ObsEWnoFU
for obsEW in `ls Zpole2/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables Zpole2 3/CorrelatedObservables Zpole2 3/g' $obsEW
      sed -i "" 's/\(Observable.*Aelectron   Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Amuon       Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Atau        Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# _C version
      sed -i "" 's/\(Observable.*Aelectron_C   Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Amuon_C       Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Atau_C        Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 3 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/1.000  0.038  0.033/# /' $obsEW
      sed -i "" '/0.038  1.000  0.007/d' $obsEW
      sed -i "" '/0.033  0.007  1.000/d' $obsEW

#     Treat also observables with LFU
      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Zpole 3
#sed -i "" 's/CorrelatedGaussianObservables Zpole3 2/CorrelatedObservables Zpole3 2/g' Zpole3/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AelectronPtau.*\) MCMC weight/\1 noMCMC noweight/g' Zpole3/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AtauPtau.*\) MCMC weight/\1 noMCMC noweight/g' Zpole3/$ObsEWnoFU
# WARNING: The next 2 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/1.000  0.012/# /' Zpole3/$ObsEWnoFU
#sed -i "" '/0.012  1.000/d' Zpole3/$ObsEWnoFU
for obsEW in `ls Zpole3/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables Zpole3 2/CorrelatedObservables Zpole3 2/g' $obsEW
      sed -i "" 's/\(Observable.*AelectronPtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AtauPtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 2 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/1.000  0.012/# /' $obsEW
      sed -i "" '/0.012  1.000/d' $obsEW

#     Treat also observables with LFU
      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

# No Zpole 4
#sed -i "" 's/CorrelatedGaussianObservables Zpole4 6/CorrelatedObservables Zpole4 6/g' Zpole4/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Rbottom.*\) MCMC weight/\1 noMCMC noweight/g' Zpole4/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Rcharm.*\) MCMC weight/\1 noMCMC noweight/g' Zpole4/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AFBbottom.*\) MCMC weight/\1 noMCMC noweight/g' Zpole4/$ObsEWnoFU
#sed -i "" 's/\(Observable.*AFBcharm.*\) MCMC weight/\1 noMCMC noweight/g' Zpole4/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Abottom.*\) MCMC weight/\1 noMCMC noweight/g' Zpole4/$ObsEWnoFU
#sed -i "" 's/\(Observable.*Acharm.*\) MCMC weight/\1 noMCMC noweight/g' Zpole4/$ObsEWnoFU
# WARNING: The next 6 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/ 1.00  -0.18  -0.10   0.07  -0.08   0.04/# /' Zpole4/$ObsEWnoFU
#sed -i "" '/-0.18   1.00   0.04  -0.06   0.04  -0.06/d' Zpole4/$ObsEWnoFU
#sed -i "" '/-0.10   0.04   1.00   0.15   0.06   0.01/d' Zpole4/$ObsEWnoFU
#sed -i "" '/ 0.07  -0.06   0.15   1.00  -0.02   0.04/d' Zpole4/$ObsEWnoFU
#sed -i "" '/-0.08   0.04   0.06  -0.02   1.00   0.11/d' Zpole4/$ObsEWnoFU
#sed -i "" '/ 0.04  -0.06   0.01   0.04   0.11   1.00/d' Zpole4/$ObsEWnoFU
for obsEW in `ls Zpole4/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables Zpole4 6/CorrelatedObservables Zpole4 6/g' $obsEW
      sed -i "" 's/\(Observable.*Rbottom.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Rcharm.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBbottom.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*AFBcharm.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Abottom.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*Acharm.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 6 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/ 1.00  -0.18  -0.10   0.07  -0.08   0.04/# /' $obsEW
      sed -i "" '/-0.18   1.00   0.04  -0.06   0.04  -0.06/d' $obsEW
      sed -i "" '/-0.10   0.04   1.00   0.15   0.06   0.01/d' $obsEW
      sed -i "" '/ 0.07  -0.06   0.15   1.00  -0.02   0.04/d' $obsEW
      sed -i "" '/-0.08   0.04   0.06  -0.02   1.00   0.11/d' $obsEW
      sed -i "" '/ 0.04  -0.06   0.01   0.04   0.11   1.00/d' $obsEW
done

# No BRWhad  # Removed from standard fits. Not independent from BRWlept
#sed -i "" 's/\(Observable.*BrWhadrons.*\) MCMC weight/\1 noMCMC noweight/g' BRWhad/$ObsEWnoFU
#for obsEW in `ls BRWhad/ObservablesEW*`
#do
#      sed -i "" 's/\(Observable.*BrWhadrons.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#done

# No BRWlept
#sed -i "" 's/CorrelatedGaussianObservables BRWlept 3/CorrelatedObservables BRWlept 3/g' BRWlept/$ObsEWnoFU
#sed -i "" 's/\(Observable.*BrWelectron.*\) MCMC weight/\1 noMCMC noweight/g' BRWlept/$ObsEWnoFU
#sed -i "" 's/\(Observable.*BrWmuon.*\) MCMC weight/\1 noMCMC noweight/g' BRWlept/$ObsEWnoFU
#sed -i "" 's/\(Observable.*BrWtau.*\) MCMC weight/\1 noMCMC noweight/g' BRWlept/$ObsEWnoFU
# WARNING: The next 3 lines must match exactly the correlations in the observables file. Update these if the corr. change
#sed -i "" 's/ 1.000  0.136 -0.201/# /' BRWlept/$ObsEWnoFU
#sed -i "" '/ 0.136  1.000 -0.122/d' BRWlept/$ObsEWnoFU
#sed -i "" '/-0.201 -0.122  1.000/d' BRWlept/$ObsEWnoFU
for obsEW in `ls BRWlept/ObservablesEW*`
do
      sed -i "" 's/CorrelatedGaussianObservables BRWlept 3/CorrelatedObservables BRWlept 3/g' $obsEW
      sed -i "" 's/\(Observable.*BrWelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*BrWmuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      sed -i "" 's/\(Observable.*BrWtau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
# WARNING: The next 3 lines must match exactly the correlations in the observables file. Update these if the corr. change
      sed -i "" 's/ 1.000  0.136 -0.201/# /' $obsEW
      sed -i "" '/ 0.136  1.000 -0.122/d' $obsEW
      sed -i "" '/-0.201 -0.122  1.000/d' $obsEW

#     Treat also observables with LFU
      sed -i "" 's/\(Observable.*BrWlepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
done

fi


# Go back to main folder (Current_Future/Fits_Current/SM_fits)
cd ..

fi

#------------------------------------------------------------------------------

# 5. MW vs mt fits

if [ $GenMWmtsinEff == 'Yes' ] || [ $GenMWmtsinEff == 'yes' ]
then

mkdir MW_mt_fits

cd MW_mt_fits

# 5.1 Make directories for all the relevant cases.

mkdir SM_noMWmt
mkdir SM_noMWmtMH

# 5.2 General operations

for dir in `ls -d */`
do

#     Put the MonteCarlo file here (use the " " to preserve the blanks)
#      (echo -e "$MCfile")>MonteCarlo-HIGH.conf

#     Copy all the general config files
      cp  ../model.conf $dir/.
      cp  ../FlavourFixed.conf $dir/.
#     (In some cases there are several ObservablesEW*.conf files)
      cp  ../ObservablesEW* $dir/.

#     Copy the cluster submission script and make adjustments if needed
      cp  ../submit_analysis_IB.sh ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/MonteCarlo-HIGH.conf/'$MCCONFLG'/' ${dir}/submit_analysis_IB.sh

done

# 5.3 Remove the corresponding measurements from each fit
for dir in `ls -d *_noMWmt`
do

#     mt
      sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtexp'      0.          100. /' $dir/model.conf
#     MW
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
      for obsEW in `ls $dir/ObservablesEW*`
      do
            sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      done

done

for dir in `ls -d *_noMWmtMH`
do

#     mt
      sed -i "" 's/ModelParameter  mtop .*/ModelParameter  mtop        '$mtexp'      0.          100. /' $dir/model.conf
#     MH
      sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         505.      0.          494./' $dir/model.conf
#     MW
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
      for obsEW in `ls $dir/ObservablesEW*`
      do
            sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      done
done

# Go back to main folder (Current_Future/Fits_Current/SM_fits)
cd ..

#------------------------------------------------------------------------------

# 6. MW vs sin2Eff fits

mkdir MW_sin2Eff_fits

cd MW_sin2Eff_fits

# 6.1 Make directories for all the relevant cases.

mkdir SM_noMWs2Eff
mkdir SM_noMWs2EffMH
mkdir SM_noMWs2EffGammaZ
mkdir SM_noMWs2EffMHGammaZ

# 6.2 General operations

for dir in `ls -d */`
do

#     Put the MonteCarlo file here (use the " " to preserve the blanks)
#      (echo -e "$MCfile")>MonteCarlo-HIGH.conf

#     Copy all the general config files
      cp  ../model.conf $dir/.
      cp  ../FlavourFixed.conf $dir/.
#     (In some cases there are several ObservablesEW*.conf files)
      cp  ../ObservablesEW* $dir/.

#     Copy the cluster submission script and make adjustments if needed
      cp  ../submit_analysis_IB.sh ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -q mpi_ib/#PBS -q '$QUEUELG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l walltime=2400:00:00/#PBS -l walltime='$WALLLG'/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/#PBS -l pmem=512mb/#PBS -l pmem='$MEMLG'mb/' ${dir}/submit_analysis_IB.sh
      sed -i "" 's/MonteCarlo-HIGH.conf/'$MCCONFLG'/' ${dir}/submit_analysis_IB.sh

done

# 6.3 Remove the corresponding measurements from each fit
for dir in `ls -d *_noMWs2Eff`
do

#     MW
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU

#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
      for obsEW in `ls $dir/ObservablesEW*`
      do
#     MW
            sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
            sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW

            sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      done
done

for dir in `ls -d *_noMWs2EffMH`
do

#     MH
      sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         505.      0.          494./' $dir/model.conf

#     MW
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU

#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
      for obsEW in `ls $dir/ObservablesEW*`
      do
#     MW
            sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
            sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW

            sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      done

done

for dir in `ls -d *_noMWs2EffGammaZ`
do

#     MW
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU

#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#     GammaZ
#      sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU

      for obsEW in `ls $dir/ObservablesEW*`
      do
#     MW
            sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
            sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW

            sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#     GammaZ
            sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      done

done

for dir in `ls -d *_noMWs2EffMHGammaZ`
do

#     MH
      sed -i "" 's/ModelParameter  mHl .*/ModelParameter  mHl         505.      0.          494./' $dir/model.conf

#     MW
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU

#      sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#      sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
#     GammaZ
#      sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWFU
#      sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $dir/$ObsEWnoFU
      for obsEW in `ls $dir/ObservablesEW*`
      do
#     MW
            sed -i "" 's/\(Observable.*Mw.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#     LEP+SLC+Tevatron+LHC determinations of sin2Eff
            sed -i "" 's/\(Observable.*PtauPol.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*sin2thetaEff.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Alepton.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW

            sed -i "" 's/\(Observable.*Aelectron.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Amuon.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
            sed -i "" 's/\(Observable.*Atau.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
#     GammaZ
            sed -i "" 's/\(Observable.*GammaZ.*\) MCMC weight/\1 noMCMC noweight/g' $obsEW
      done

done

# Go back to main folder (Current_Future/Fits_Current/SM_fits)
cd ..

fi

#------------------------------------------------------------------------------

# Go back to Current_Future
cd ../..

done

#------------------------------------------------------------------------------
