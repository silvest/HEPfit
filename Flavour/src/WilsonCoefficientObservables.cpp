/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "WilsonCoefficientObservables.h"
#include "StandardModel.h"

/*******************************************************************************
 * Wilson Coefficients                                                         *
 * ****************************************************************************/


WC_C7g::WC_C7g(const StandardModel& SM_i, unsigned int part_i) 
: ThObservable(SM_i), mySM(SM_i) 
{  
    part = part_i;
    mu = mySM.getMub();
}

double WC_C7g::computeThValue()
{
    mu = getBinMin();
    allcoeff = mySM.getFlavour().ComputeCoeffBMll(mu, QCD::MU); /*QCD::MU does not matter here but necessary.*/
    
    if (part == 0) return ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6)).real();
    if (part == 1) return ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6)).imag();
    if (part == 2) return ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6)).abs();
    if (part == 3) return ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6)).arg();
    else throw std::runtime_error("WC_C7g: part specification can only be 0, 1 or 2");
}

WC_C9::WC_C9(const StandardModel& SM_i, unsigned int part_i, QCD::lepton lep_i) 
: ThObservable(SM_i), mySM(SM_i) 
{  
    part = part_i;
    lepton = lep_i;
    mu = mySM.getMub();
}

double WC_C9::computeThValue()
{
    mu = getBinMin();
    allcoeff = mySM.getFlavour().ComputeCoeffBMll(mu, lepton);
    
    if (part == 0) return ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8)).real();
    if (part == 1) return ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8)).imag();
    if (part == 2) return ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8)).abs();
    if (part == 3) return ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8)).arg();
    else throw std::runtime_error("WC_C7g: part specification can only be 0, 1 or 2");
}

WC_C10::WC_C10(const StandardModel& SM_i, unsigned int part_i, QCD::lepton lep_i) 
: ThObservable(SM_i), mySM(SM_i) 
{  
    part = part_i;
    lepton = lep_i;
    mu = mySM.getMub();
}

double WC_C10::computeThValue()
{
    mu = getBinMin();
    allcoeff = mySM.getFlavour().ComputeCoeffBMll(mu, lepton);
    
    if (part == 0) return ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9)).real();
    if (part == 1) return ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9)).imag();
    if (part == 2) return ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9)).abs();
    if (part == 3) return ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9)).arg();
    else throw std::runtime_error("WC_C7g: part specification can only be 0, 1 or 2");
}