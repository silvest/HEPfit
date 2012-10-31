/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2Rcharm.h"


LEP2Rcharm::LEP2Rcharm(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myLEP2oblique(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2Rcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, bDP, bWEAK, bQED)) {
        double sigma_c, sigma_had;
        sigma_c = SM.sigma_q_LEP2(StandardModel::CHARM, 
                                            s, Mw, GammaZ, bDP, bWEAK, bQED);
        sigma_had = SM.sigma_q_LEP2(StandardModel::UP, 
                                              s, Mw, GammaZ, bDP, bWEAK, bQED)
                  + SM.sigma_q_LEP2(StandardModel::DOWN, 
                                              s, Mw, GammaZ, bDP, bWEAK, bQED)
                  + sigma_c
                  + SM.sigma_q_LEP2(StandardModel::STRANGE, 
                                              s, Mw, GammaZ, bDP, bWEAK, bQED)
                  + SM.sigma_q_LEP2(StandardModel::BOTTOM,
                                              s, Mw, GammaZ, bDP, bWEAK, bQED);
        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return R_c;
}


