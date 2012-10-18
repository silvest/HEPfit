/* 
 * File:   LEP2AFBmu.cpp
 * Author: mishima
 */

#include "LEP2AFBmu.h"


double LEP2AFBmu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ))
        SMresult_cache = myTwoFermions.AFB_l(StandardModel::MU, s, Mw, GammaZ, bRCs);
    double AFB_mu = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_mu += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::MU, s);
    }
    
    return AFB_mu;
}
        
