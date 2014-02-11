/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_lj_gamma.h"

li_lj_gamma::li_lj_gamma(LeptonFlavour& LeptonFlavour_i): ThObservable(LeptonFlavour_i), myLeptonFlavour(LeptonFlavour_i){
};

double li_lj_gamma::computeThValue(){
    
    gslpp::vector<complex> ** allcoeff = myLeptonFlavour.ComputeCoeffli_lj_gamma();
    
    return (1/10. * ((*(allcoeff[LO])) * (*(allcoeff[LO])).conjugate()).abs());
}