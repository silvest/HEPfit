/* 
 * File:   LEP2AFBbottom.cpp
 * Author: mishima
 */

#include "LEP2AFBbottom.h"


double LEP2AFBbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bRCs_cache, s, Mw, GammaZ, bRCs))
        SMresult_cache = SM.AFB_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
    double AFB_b = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_b += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::BOTTOM, s);      
    }
    
    return AFB_b;
}
        


