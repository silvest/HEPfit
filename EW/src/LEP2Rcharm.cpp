/* 
 * File:   LEP2Rcharm.cpp
 * Author: mishima
 */

#include "LEP2Rcharm.h"


double LEP2Rcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        double sigma_c, sigma_had;
        sigma_c = myTwoFermions.sigma_q(StandardModel::CHARM, s, Mw, GammaZ, bRCs);
        sigma_had = myTwoFermions.sigma_q(StandardModel::UP, s, Mw, GammaZ, bRCs)
                  + myTwoFermions.sigma_q(StandardModel::DOWN, s, Mw, GammaZ, bRCs)
                  + sigma_c
                  + myTwoFermions.sigma_q(StandardModel::STRANGE, s, Mw, GammaZ, bRCs)
                  + myTwoFermions.sigma_q(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return R_c;
}


