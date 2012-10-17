/* 
 * File:   LEP2AFBtau.cpp
 * Author: mishima
 */

#include "LEP2AFBtau.h"


double LEP2AFBtau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bRCs_cache, s, Mw, GammaZ, bRCs))
        SMresult_cache = SM.AFB_l_LEP2(StandardModel::TAU, s, Mw, GammaZ, bRCs);
    double AFB_tau = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_tau += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::TAU, s);
    }
    
    return AFB_tau;
}
        