/* 
 * File:   LEP2sigmaTau.cpp
 * Author: mishima
 */

#include "LEP2sigmaTau.h"


double LEP2sigmaTau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, Flags))
        SMresult_cache = SM.sigma_l_LEP2(StandardModel::TAU, 
                                                   s, Mw, GammaZ, Flags);
    double sigma_tau = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        sigma_tau += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::TAU, s);
    }
    
    return ( sigma_tau*GeVminus2_to_nb*1000. );
}
        

