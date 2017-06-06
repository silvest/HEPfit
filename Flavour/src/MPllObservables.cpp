/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MPllObservables.h"
#include "StandardModel.h"

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MPll::BR_MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    if (SM.getModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to P l+ l-: BR_MPll not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep).initializeMPllParameters());
}

double BR_MPll::computeBR_MPll(double qmin, double qmax, QCD::lepton lep) 
{
    double q_min = qmin;
    double q_max = qmax;
    QCD::lepton lep_i = lep;
    
    return (3.*SM.getFlavour().getMPll(meson, pseudoscalar, lep_i).integrateSigma(0,q_min,q_max) - SM.getFlavour().getMPll(meson, pseudoscalar, lep_i).integrateSigma(2,q_min,q_max))/(4. * SM.getFlavour().getMPll(meson, pseudoscalar, lep_i).getwidth());
}

double BR_MPll::computeThValue()
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeBR_MPll(q_min, q_max, lep) / ( q_max - q_min );
}

R_MPll::R_MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2) 
: BR_MPll(SM_i, meson_i, pseudoscalar_i, lep_1) 
{  
    if (SM.getModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to P l+ l-: R_MPll not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep1).initializeMPllParameters());
}

double R_MPll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeBR_MPll(q_min, q_max, lep1)/computeBR_MPll(q_min, q_max, lep2);
}

ACP_MPll::ACP_MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) 
: BR_MPll(SM_i, meson_i, pseudoscalar_i, lep_i) 
{
    if (SM.getModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to P l+ l-: ACP not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep).initializeMPllParameters());
}

double ACP_MPll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (3.*SM.getFlavour().getMPll(meson, pseudoscalar, lep).integrateDelta(0, q_min, q_max) - SM.getFlavour().getMPll(meson, pseudoscalar, lep).integrateDelta(2, q_min, q_max))/(4.*computeBR_MPll(q_min, q_max, lep)* SM.getFlavour().getMPll(meson, pseudoscalar, lep).getwidth());
}

