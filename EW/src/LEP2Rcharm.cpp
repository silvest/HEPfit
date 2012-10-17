/* 
 * File:   LEP2Rcharm.cpp
 * Author: mishima
 */

#include "LEP2Rcharm.h"


double LEP2Rcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bRCs_cache, s, Mw, GammaZ, bRCs)) {
        double sigma_c, sigma_had;
        sigma_c = SM.sigma_q_LEP2(StandardModel::CHARM, s, Mw, GammaZ, bRCs);
        sigma_had = SM.sigma_q_LEP2(StandardModel::UP, s, Mw, GammaZ, bRCs)
                  + SM.sigma_q_LEP2(StandardModel::DOWN, s, Mw, GammaZ, bRCs)
                  + sigma_c
                  + SM.sigma_q_LEP2(StandardModel::STRANGE, s, Mw, GammaZ, bRCs)
                  + SM.sigma_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return R_c;
}


