/* 
 * File:   LEP2Rcharm.cpp
 * Author: mishima
 */

#include "LEP2Rcharm.h"


double LEP2Rcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, Flags)) {
        double sigma_c, sigma_had;
        sigma_c = SM.sigma_q_LEP2(StandardModel::CHARM, 
                                            s, Mw, GammaZ, Flags);
        sigma_had = SM.sigma_q_LEP2(StandardModel::UP, 
                                              s, Mw, GammaZ, Flags)
                  + SM.sigma_q_LEP2(StandardModel::DOWN, 
                                              s, Mw, GammaZ, Flags)
                  + sigma_c
                  + SM.sigma_q_LEP2(StandardModel::STRANGE, 
                                              s, Mw, GammaZ, Flags)
                  + SM.sigma_q_LEP2(StandardModel::BOTTOM,
                                              s, Mw, GammaZ, Flags);
        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return R_c;
}


