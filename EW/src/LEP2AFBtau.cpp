/* 
 * File:   LEP2AFBtau.cpp
 * Author: mishima
 */

#include "LEP2AFBtau.h"


LEP2AFBtau::LEP2AFBtau(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myLEP2oblique(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2AFBtau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!myEW.getSM().getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, bDP, bWEAK, bQED))
        SMresult_cache = myEW.getSM().AFB_l_LEP2(StandardModel::TAU, 
                                                 s, Mw, GammaZ, bDP, bWEAK, bQED);
    double AFB_tau = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_tau += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::TAU, s);
    }
    
    return AFB_tau;
}
        