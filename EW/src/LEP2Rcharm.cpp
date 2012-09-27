/* 
 * File:   LEP2Rcharm.cpp
 * Author: mishima
 */

#include "LEP2Rcharm.h"


LEP2Rcharm::LEP2Rcharm(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2Rcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double sigma_c, sigma_had;
    sigma_c = myEW.getSM().sigma_q_LEP2(StandardModel::CHARM, s, Mw, GammaZ, 
                                        bDP, bWEAK, bQED);
    sigma_had = myEW.getSM().sigma_q_LEP2(StandardModel::UP, s, Mw, GammaZ, 
                                          bDP, bWEAK, bQED)
                + myEW.getSM().sigma_q_LEP2(StandardModel::DOWN, s, Mw, GammaZ, 
                                            bDP, bWEAK, bQED)
                + myEW.getSM().sigma_q_LEP2(StandardModel::CHARM, s, Mw, GammaZ, 
                                            bDP, bWEAK, bQED)
                + myEW.getSM().sigma_q_LEP2(StandardModel::STRANGE, s, Mw, GammaZ, 
                                            bDP, bWEAK, bQED)
                + myEW.getSM().sigma_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ, 
                                            bDP, bWEAK, bQED);
    double R_c = sigma_c/sigma_had;
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        R_c += 0.0;       
    }
    
    return R_c;
}


