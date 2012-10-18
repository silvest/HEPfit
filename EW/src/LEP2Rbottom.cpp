/* 
 * File:   LEP2Rbottom.cpp
 * Author: mishima
 */

#include "LEP2Rbottom.h"


double LEP2Rbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        double sigma_b, sigma_had;
        sigma_b = myTwoFermions.sigma_q(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
        sigma_had = myTwoFermions.sigma_q(StandardModel::UP, s, Mw, GammaZ, bRCs)
                  + myTwoFermions.sigma_q(StandardModel::DOWN, s, Mw, GammaZ, bRCs)
                  + myTwoFermions.sigma_q(StandardModel::CHARM, s, Mw, GammaZ, bRCs)
                  + myTwoFermions.sigma_q(StandardModel::STRANGE, s, Mw, GammaZ, bRCs)
                  + sigma_b;
        SMresult_cache = sigma_b/sigma_had;
    }
    double R_b = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_b += myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s);
    }
    
    return R_b;
}
        




