/* 
 * File:   LEP2AFBmu.cpp
 * Author: mishima
 */

#include "LEP2AFBmu.h"


LEP2AFBmu::LEP2AFBmu(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2AFBmu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double AFB_mu = myEW.getSM().AFB_l_LEP2(StandardModel::MU, s, Mw, GammaZ, 
                                            bDP, bWEAK, bQED);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        AFB_mu += 0.0;       
    }
    
    return AFB_mu;
}
        
