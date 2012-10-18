/* 
 * File:   LEP2sigmaTau.cpp
 * Author: mishima
 */

#include "LEP2sigmaTau.h"


double LEP2sigmaTau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ))
        SMresult_cache = myTwoFermions.sigma_l(StandardModel::TAU, s, Mw, GammaZ, bRCs);
    double sigma_tau = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        sigma_tau += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::TAU, s);
    }
    
    return ( sigma_tau*GeVminus2_to_nb*1000. );
}
        

