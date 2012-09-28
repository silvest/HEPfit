/* 
 * File:   LEP2sigmaHadron.cpp
 * Author: mishima
 */

#include "LEP2sigmaHadron.h"


LEP2sigmaHadron::LEP2sigmaHadron(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myLEP2oblique(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2sigmaHadron::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!myEW.getSM().getEWSM()->checkForLEP2(SMparams_cache, bool_cache,
                                              s, Mw, GammaZ, bDP, bWEAK, bQED))
        SMresult_cache = myEW.getSM().sigma_q_LEP2(StandardModel::UP, 
                                                   s, Mw, GammaZ, bDP, bWEAK, bQED)
                       + myEW.getSM().sigma_q_LEP2(StandardModel::DOWN, 
                                                   s, Mw, GammaZ, bDP, bWEAK, bQED)
                       + myEW.getSM().sigma_q_LEP2(StandardModel::CHARM, 
                                                   s, Mw, GammaZ, bDP, bWEAK, bQED)
                       + myEW.getSM().sigma_q_LEP2(StandardModel::STRANGE, 
                                                   s, Mw, GammaZ, bDP, bWEAK, bQED)
                       + myEW.getSM().sigma_q_LEP2(StandardModel::BOTTOM, 
                                                   s, Mw, GammaZ, bDP, bWEAK, bQED);
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
        

