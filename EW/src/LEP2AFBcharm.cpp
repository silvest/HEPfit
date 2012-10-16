/* 
 * File:   LEP2AFBcharm.cpp
 * Author: mishima
 */

#include "LEP2AFBcharm.h"


double LEP2AFBcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, Flags))
        SMresult_cache = SM.AFB_q_LEP2(StandardModel::CHARM, 
                                                 s, Mw, GammaZ, Flags);
    double AFB_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_c += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return AFB_c;
}
        


