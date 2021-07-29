/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MPllObservables.h"
#include "MPll.h"
#include "StandardModel.h"

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MPll::BR_MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
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
    
    return computeBR_MPll(q_min, q_max, lep);
}

dBR_MPll::dBR_MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) 
: BR_MPll(SM_i, meson_i, pseudoscalar_i, lep_i) 
{  
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep).initializeMPllParameters());
}

double dBR_MPll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeBR_MPll(q_min, q_max, lep) / ( q_max - q_min );
}

R_MPll::R_MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2) 
: BR_MPll(SM_i, meson_i, pseudoscalar_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep1).initializeMPllParameters());
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep2).initializeMPllParameters());
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

DC9_hlambda::DC9_hlambda(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPll(meson, pseudoscalar, lep).initializeMPllParameters());
}

double DC9_hlambda::computeThValue() 
{   
    double q2 = getBinMin();
    double sixteenM_PI2 = 16.*M_PI*M_PI;
    gslpp::complex hlambda = SM.getFlavour().getMPll(meson, pseudoscalar, lep).h_lambda(q2);
    double MM = SM.getMesons(meson).getMass();
    gslpp::complex VL = SM.getFlavour().getMPll(meson, pseudoscalar, lep).V_L(q2);
    gslpp::complex result = - sixteenM_PI2 * hlambda * MM * MM / q2 / VL;
    return result.abs();
}