/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MPllObservables.h"
#include "Flavour.h"

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MPll::BR_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    if (SM.ModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to P l+ l-: BR_MPll not implemented in: " + SM.ModelName() + " model, returning Standard Model value.\n" << std::endl;
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
}

double BR_MPll::computeBR_MPll(double qmin, double qmax, StandardModel::lepton lep) 
{
    double q_min = qmin;
    double q_max = qmax;
    StandardModel::lepton lep_i = lep;
    
    return (3.*SM.getMyFlavour()->getMPll(meson, pseudoscalar, lep_i)->integrateSigma(0,q_min,q_max) - SM.getMyFlavour()->getMPll(meson, pseudoscalar, lep_i)->integrateSigma(2,q_min,q_max))/(4. * SM.getMyFlavour()->getMPll(meson, pseudoscalar, lep_i)->width);
}

double BR_MPll::computeThValue()
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeBR_MPll(q_min, q_max, lep) / ( q_max - q_min );
}

R_MPll::R_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2) 
: BR_MPll(SM_i, meson_i, pseudoscalar_i, lep_1) 
{  
    if (SM.ModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to P l+ l-: R_MPll not implemented in: " + SM.ModelName() + " model, returning Standard Model value.\n" << std::endl;
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
}

double R_MPll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeBR_MPll(q_min, q_max, lep1)/computeBR_MPll(q_min, q_max, lep2);
}

ACP_MPll::ACP_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i) 
: BR_MPll(SM_i, meson_i, pseudoscalar_i, lep_i) 
{
    if (SM.ModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to P l+ l-: ACP not implemented in: " + SM.ModelName() + " model, returning Standard Model value.\n" << std::endl;
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
}

double ACP_MPll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (3.*SM.getMyFlavour()->getMPll(meson, pseudoscalar, lep)->integrateDelta(0, q_min, q_max) - SM.getMyFlavour()->getMPll(meson, pseudoscalar, lep)->integrateDelta(2, q_min, q_max))/(4.*computeBR_MPll(q_min, q_max, lep)* SM.getMyFlavour()->getMPll(meson, pseudoscalar, lep)->width);
}

