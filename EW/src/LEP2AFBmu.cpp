/* 
 * File:   LEP2AFBmu.cpp
 * Author: mishima
 */

#include "LEP2AFBmu.h"


double LEP2AFBmu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, Flags))
        SMresult_cache = SM.AFB_l_LEP2(StandardModel::MU, 
                                                 s, Mw, GammaZ, Flags);
    double AFB_mu = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_mu += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::MU, s);
    }
    
    return AFB_mu;
}
        
