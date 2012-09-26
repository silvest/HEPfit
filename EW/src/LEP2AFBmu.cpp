/* 
 * File:   LEP2AFBmu.cpp
 * Author: mishima
 */

#include "LEP2AFBmu.h"


double LEP2AFBmu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double AFB_mu = myEW.getSM().AFB_l_LEP2(StandardModel::MU, s, Mw, GammaZ);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        AFB_mu += 0.0;       
    }
    
    return AFB_mu;
}
        
