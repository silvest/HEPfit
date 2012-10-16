/* 
 * File:   LEP2Rbottom.cpp
 * Author: mishima
 */

#include "LEP2Rbottom.h"


double LEP2Rbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, Flags)) {
        double sigma_b, sigma_had;
        sigma_b = SM.sigma_q_LEP2(StandardModel::BOTTOM, 
                                            s, Mw, GammaZ, Flags);
        sigma_had = SM.sigma_q_LEP2(StandardModel::UP, 
                                              s, Mw, GammaZ, Flags)
                  + SM.sigma_q_LEP2(StandardModel::DOWN, 
                                              s, Mw, GammaZ, Flags)
                  + SM.sigma_q_LEP2(StandardModel::CHARM, 
                                              s, Mw, GammaZ, Flags)
                  + SM.sigma_q_LEP2(StandardModel::STRANGE, 
                                              s, Mw, GammaZ, Flags)
                  + sigma_b;
        SMresult_cache = sigma_b/sigma_had;
    }
    double R_b = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_b += myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s);
    }
    
    return R_b;
}
        




