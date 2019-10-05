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
    
    mcdbd2.setMu(FWC.GetWCscale_bd());
    mcdbd2.setCoeff(0, FWC.GetC1_bd(), LO);
    mcdbd2.setCoeff(1, FWC.GetC2_bd(), LO);
    mcdbd2.setCoeff(2, FWC.GetC3_bd(), LO);
    mcdbd2.setCoeff(3, FWC.GetC4_bd(), LO);
    mcdbd2.setCoeff(4, FWC.GetC5_bd(), LO);
    
    vmcdbd2.push_back(mcdbd2);
    
    return(vmcdbd2);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdbs2() {
    vmcdbs2.clear();
    vmcdbs2 = StandardModelMatching::CMdbs2();
    
    mcdbs2.setMu(FWC.GetWCscale_bs());
    mcdbs2.setCoeff(0, FWC.GetC1_bs(), LO);
    mcdbs2.setCoeff(1, FWC.GetC2_bs(), LO);
    mcdbs2.setCoeff(2, FWC.GetC3_bs(), LO);
    mcdbs2.setCoeff(3, FWC.GetC4_bs(), LO);
    mcdbs2.setCoeff(4, FWC.GetC5_bs(), LO);
    
    vmcdbs2.push_back(mcdbs2);
    
    return(vmcdbs2);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdd2() {
    vmcdc2.clear();
    vmcdc2 = StandardModelMatching::CMdd2();
    
    mcdc2.setMu(FWC.GetWCscale_c());
    mcdc2.setCoeff(0, FWC.GetC1_c(), LO);
    mcdc2.setCoeff(1, FWC.GetC2_c(), LO);
    mcdc2.setCoeff(2, FWC.GetC3_c(), LO);
    mcdc2.setCoeff(3, FWC.GetC4_c(), LO);
    mcdc2.setCoeff(4, FWC.GetC5_c(), LO);
    
    vmcdc2.push_back(mcdc2);
    
    return(vmcdc2);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficient_DF2Matching::CMdk2() {
    vmcds2.clear();
    vmcds2 = StandardModelMatching::CMdk2();
    
    mcds2.setMu(FWC.GetWCscale_s());
    mcds2.setCoeff(0, FWC.GetC1_s(), LO);
    mcds2.setCoeff(1, FWC.GetC2_s(), LO);
    mcds2.setCoeff(2, FWC.GetC3_s(), LO);
    mcds2.setCoeff(3, FWC.GetC4_s(), LO);
    mcds2.setCoeff(4, FWC.GetC5_s(), LO);
    
    vmcds2.push_back(mcds2);
    
    return(vmcds2);
}
