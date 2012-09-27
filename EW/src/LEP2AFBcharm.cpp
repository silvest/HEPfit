/* 
 * File:   LEP2AFBcharm.cpp
 * Author: mishima
 */

#include "LEP2AFBcharm.h"


LEP2AFBcharm::LEP2AFBcharm(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2AFBcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double AFB_c = myEW.getSM().AFB_q_LEP2(StandardModel::CHARM, s, Mw, GammaZ, 
                                           bDP, bWEAK, bQED);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        AFB_c += 0.0;       
    }
    
    return AFB_c;
}
        


