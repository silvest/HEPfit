/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_lj_gamma.h"

li_lj_gamma::li_lj_gamma(const StandardModel& SM_i): ThObservable(SM_i), mySM(SM_i){
};

double li_lj_gamma::computeThValue(){
    
    gslpp::vector<complex> ** allcoeff = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma();
    
    return (1/10. * ((*(allcoeff[LO])) * (*(allcoeff[LO])).conjugate()).abs());
}