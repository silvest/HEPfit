/* 
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficient_DF2Matching.h"
#include "FlavourWilsonCoefficient_DF2.h" 

FlavourWilsonCoefficient_DF2Matching::FlavourWilsonCoefficient_DF2Matching(const FlavourWilsonCoefficient_DF2& FWC_i) :
    StandardModelMatching(FWC_i), FWC(FWC_i), 
    mcdbd2(5, NDR, NLO), mcdbs2(5, NDR, NLO), 
    mcdc2(5, NDR, NLO), mcds2(5, NDR, NLO) {}   
        
std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdbd2() {
    vmcdbd2.clear();
    vmcdbd2 = StandardModelMatching::CMdbd2();
    
    mcdbd2.setMu(FWC.getWCscale_bd());
    for (int k = 0; k < 5; k++) 
        mcdbd2.setCoeff(k, FWC.getC_bd()(k), LO);
   
    vmcdbd2.push_back(mcdbd2);
    
    return(vmcdbd2);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdbs2() {
    vmcdbs2.clear();
    vmcdbs2 = StandardModelMatching::CMdbs2();
    
    mcdbs2.setMu(FWC.getWCscale_bs());
    for (int k = 0; k < 5; k++) 
        mcdbs2.setCoeff(k, FWC.getC_bs()(k), LO);
    
    vmcdbs2.push_back(mcdbs2);
    
    return(vmcdbs2);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdd2() {
    vmcdc2.clear();
    vmcdc2 = StandardModelMatching::CMdd2();
    
    mcdc2.setMu(FWC.getWCscale_c());
    for (int k = 0; k < 5; k++) 
        mcdc2.setCoeff(k, FWC.getC_c()(k), LO);
    
    vmcdc2.push_back(mcdc2);
    
    return(vmcdc2);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdk2() {
    vmcds2.clear();
    vmcds2 = StandardModelMatching::CMdk2();
    
    mcds2.setMu(FWC.getWCscale_s());
    for (int k = 0; k < 5; k++) 
        mcds2.setCoeff(k, FWC.getC_s()(k), LO);
    
    vmcds2.push_back(mcds2);
    
    return(vmcds2);
}
