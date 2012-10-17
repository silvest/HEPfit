/* 
 * File:   LEP2Rbottom.cpp
 * Author: mishima
 */

#include "LEP2Rbottom.h"


double LEP2Rbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bRCs_cache, s, Mw, GammaZ, bRCs)) {
        double sigma_b, sigma_had;
        sigma_b = SM.sigma_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
        sigma_had = SM.sigma_q_LEP2(StandardModel::UP, s, Mw, GammaZ, bRCs)
                  + SM.sigma_q_LEP2(StandardModel::DOWN, s, Mw, GammaZ, bRCs)
                  + SM.sigma_q_LEP2(StandardModel::CHARM, s, Mw, GammaZ, bRCs)
                  + SM.sigma_q_LEP2(StandardModel::STRANGE, s, Mw, GammaZ, bRCs)
                  + sigma_b;
        SMresult_cache = sigma_b/sigma_had;
    }
    double R_b = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_b += myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s);
    }
    
    return R_b;
}
        




