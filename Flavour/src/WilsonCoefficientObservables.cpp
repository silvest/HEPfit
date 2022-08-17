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

WC_epspOeps::WC_epspOeps(const StandardModel& SM_i, unsigned int part_i, double mu_z) 
: ThObservable(SM_i), mySM(SM_i) 
{  
    part = part_i;
    mu = mu_z;
}

double WC_epspOeps::computeThValue()
{
    allcoeffv = mySM.getFlavour().ComputeCoeffDS1PPv(mu, NDR);
    allcoeffz = mySM.getFlavour().ComputeCoeffDS1PPz(mu, NDR);
    gslpp::vector<gslpp::complex> allcoeff_z = (*allcoeffz[LO]) + (*allcoeffz[LO_QED]) + (*allcoeffz[NLO]) + (*allcoeffz[NLO_QED11]);
    gslpp::vector<gslpp::complex> allcoeff_y = (*allcoeffv[LO]) + (*allcoeffv[LO_QED]) + (*allcoeffv[NLO]) + (*allcoeffv[NLO_QED11]);
    allcoeff_z.assign(0,allcoeff_y(0));
    allcoeff_z.assign(1,allcoeff_y(1));
    for(int i = 2; i<10; i++){
        allcoeff_y.assign(i,allcoeff_y(i)-allcoeff_z(i));
    }
    // Table 5 of 1909.05610
    if (part == 0) return allcoeff_z(0).real(); // z1
    if (part == 1) return allcoeff_z(1).real(); // z2
    if (part == 2) return allcoeff_y(2).real(); // y3
    if (part == 3) return allcoeff_y(3).real(); // y4
    if (part == 4) return allcoeff_y(4).real(); // y5
    if (part == 5) return allcoeff_y(5).real(); // y6
    if (part == 6) return allcoeff_y(6).real(); // y7
    if (part == 7) return allcoeff_y(7).real(); // y8
    if (part == 8) return allcoeff_y(8).real(); // y9
    if (part == 9) return allcoeff_y(9).real(); // y10
    else throw std::runtime_error("WC_epspOeps: part specification can only go from 0 to 9");
}
