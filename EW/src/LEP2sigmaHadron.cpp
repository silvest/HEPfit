/* 
 * File:   LEP2sigmaHadron.cpp
 * Author: mishima
 */

#include "LEP2sigmaHadron.h"


double LEP2sigmaHadron::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, Flags))
        SMresult_cache = SM.sigma_q_LEP2(StandardModel::UP, 
                                                   s, Mw, GammaZ, Flags)
                       + SM.sigma_q_LEP2(StandardModel::DOWN, 
                                                   s, Mw, GammaZ, Flags)
                       + SM.sigma_q_LEP2(StandardModel::CHARM, 
                                                   s, Mw, GammaZ, Flags)
                       + SM.sigma_q_LEP2(StandardModel::STRANGE, 
                                                   s, Mw, GammaZ, Flags)
                       + SM.sigma_q_LEP2(StandardModel::BOTTOM, 
                                                   s, Mw, GammaZ, Flags);
    double sigmaH = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        sigmaH += myLEP2oblique.sigma_q_LEP2_NP(StandardModel::UP, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::DOWN, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::CHARM, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::STRANGE, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::BOTTOM, s);
    }
    
    return ( sigmaH*GeVminus2_to_nb*1000. );
}
        

