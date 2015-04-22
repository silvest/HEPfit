#!/bin/sh

# Usage: 
#  > sh mod.sh ModelFactory.cpp EWObservables.h ThObsFactory.cpp StandardModel.cpp StandardModel.h

MFCPP=$1
EWOH=$2
TOFCPP=$3
SMCPP=$4
SMH=$5


#########  ModelFactory.cpp  ##########

ORG="\QDUMMY\E"
NEW="DUMMY"
perl -p -w -i.bak -e "s/${ORG}/${NEW}/g" $MFCPP

MODELLIST="NPSTUVWXY MFV GeneralSUSY pMSSM SUSY SUSYMassInsertion THDM"
for MODEL in $MODELLIST
do
    ORG="\Q#include <${MODEL}.h>\E"
    NEW="\/\/#include <${MODEL}.h>"
    perl -p -w -i -e "s/${ORG}/${NEW}/g" $MFCPP
    ORG="\QmodelFactory[\"${MODEL}\E"
    NEW="\/\/modelFactory[\"${MODEL}"
    perl -p -w -i -e "s/${ORG}/${NEW}/g" $MFCPP
done


#########  EWObservables.h  #########

ORG="\Q#include \"LEP2\E"
NEW="\/\/#include \"LEP2"
perl -p -w -i.bak -e "s/${ORG}/${NEW}/g" $EWOH


#########  ThObsFactory.cpp  #########

ORG="\Q#include <LeptonFlavourObservables.h>\E"
NEW="\/\/#include <LeptonFlavourObservables.h>"
perl -p -w -i.bak -e "s/${ORG}/${NEW}/g" $TOFCPP

ORG="\Q#include <SUSYObservables.h>\E"
NEW="\/\/#include <SUSYObservables.h>"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $TOFCPP

OBSLIST="M12D ArgD EpsiloP_o_Epsilon li_lj_gamma OutputSLHAfromFH MHl MHh MHa MHp Msu Msd Mch Mneu Mw_dRho sigmaqLEP2 sigmamuLEP2 sigmatauLEP2 AFBmuLEP2 AFBtauLEP2 AFBbottomLEP2 AFBcharmLEP2 RbottomLEP2 RcharmLEP2"
for OBS in $OBSLIST
do
    ORG="\QobsThFactory[\"${OBS}\E"
    NEW="\/\/obsThFactory[\"${OBS}"
    perl -p -w -i -e "s/${ORG}/${NEW}/g" $TOFCPP
done


#########  StandardModel.cpp  #########

ORG="\Q#include <LeptonFlavour.h>\E"
NEW="\/\/#include <LeptonFlavour.h>"
perl -p -w -i.bak -e "s/${ORG}/${NEW}/g" $SMCPP

ORG="\Q#include \"EWSMTwoFermionsLEP2.h\"\E"
NEW="\/\/#include \"EWSMTwoFermionsLEP2.h"\"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMCPP

ORG="\Qif (myLeptonFlavour\E"
NEW="\/\/if (myLeptonFlavour"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMCPP

ORG="\Qif (myTwoFermionsLEP2\E"
NEW="\/\/if (myTwoFermionsLEP2"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMCPP

ORG="\QmyLeptonFlavour\E"
NEW="\/\/myLeptonFlavour"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMCPP

ORG="\QmyTwoFermionsLEP2\E"
NEW="\/\/myTwoFermionsLEP2"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMCPP


#########  StandardModel.h  #########

ORG="\Qclass EWSMTwoFermionsLEP2\E"
NEW="\/\/class EWSMTwoFermionsLEP2"
perl -p -w -i.bak -e "s/${ORG}/${NEW}/g" $SMH

ORG="\Qclass LeptonFlavour\E"
NEW="\/\/class LeptonFlavour"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH

ORG="EWSMTwoFermionsLEP2\* getMyTwoFermionsLEP2\(\) const\n"
NEW="\/\/EWSMTwoFermionsLEP2\* getMyTwoFermionsLEP2\(\) const\n    \/\/"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH

ORG="    return myTwoFermionsLEP2;\n"
NEW="\/\/    return myTwoFermionsLEP2"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH

ORG="LeptonFlavour\* getMyLeptonFlavour\(\) const\n"
NEW="\/\/LeptonFlavour\* getMyLeptonFlavour\(\) const\n    \/\/"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH

ORG="    return myLeptonFlavour;\n"
NEW="\/\/    return myLeptonFlavour"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH

ORG="\QEWSMTwoFermionsLEP2*\E"
NEW="\/\/EWSMTwoFermionsLEP2*"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH

ORG="\QLeptonFlavour*\E"
NEW="\/\/LeptonFlavour*"
perl -p -w -i -e "s/${ORG}/${NEW}/g" $SMH
