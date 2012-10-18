/* 
 * File:   LEP2AFBbottom.cpp
 * Author: mishima
 */

#include "LEP2AFBbottom.h"


double LEP2AFBbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ))
        SMresult_cache = myTwoFermions.AFB_q(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
    double AFB_b = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_b += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::BOTTOM, s);      
    }
    
    return AFB_b;
}
        


