/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDL1.h"
#include "gslpp_complex.h"

using namespace gslpp;

HeffDL1::HeffDL1(const StandardModel & SM_i) :
        model(SM_i),
        coeffDL1(2, NDR , LO){
}

HeffDL1::~HeffDL1() {
}

vector<complex>** HeffDL1::ComputeCoeffDL1() {
    
    std::vector<WilsonCoefficient>& mcb = model.getMyMatching() -> CMDL1();
    
    orders ordDL1 = coeffDL1.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDL1; j++){
            coeffDL1.setCoeff(*coeffDL1.getCoeff(orders(j))
                                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }
     
    return coeffDL1.getCoeff(); 
} 




