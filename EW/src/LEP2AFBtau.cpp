/* 
 * File:   LEP2AFBtau.cpp
 * Author: mishima
 */

#include "LEP2AFBtau.h"


double LEP2AFBtau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double AFB_tau = myEW.getSM().AFB_l_LEP2(StandardModel::TAU, s, Mw, GammaZ);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        AFB_tau += 0.0;       
    }
    
    return AFB_tau;
}
        