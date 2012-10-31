/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2AFBmu.h"


LEP2AFBmu::LEP2AFBmu(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myLEP2oblique(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2AFBmu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, bDP, bWEAK, bQED))
        SMresult_cache = SM.AFB_l_LEP2(StandardModel::MU, 
                                                 s, Mw, GammaZ, bDP, bWEAK, bQED);
    double AFB_mu = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_mu += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::MU, s);
    }
    
    return AFB_mu;
}
        
