/* 
 * File:   LEP2AFBcharm.cpp
 * Author: mishima
 */

#include "LEP2AFBcharm.h"


LEP2AFBcharm::LEP2AFBcharm(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myLEP2oblique(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2AFBcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!myEW.getSM().getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, bDP, bWEAK, bQED))
        SMresult_cache = myEW.getSM().AFB_q_LEP2(StandardModel::CHARM, 
                                                 s, Mw, GammaZ, bDP, bWEAK, bQED);
    double AFB_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_c += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return AFB_c;
}
        


