/* 
 * File:   LEP2sigmaTau.cpp
 * Author: mishima
 */

#include "LEP2sigmaTau.h"


double LEP2sigmaTau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double sigma_tau = myEW.getSM().sigma_l_LEP2(StandardModel::TAU, s, Mw, GammaZ);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        sigma_tau += 0.0;       
    }
    
    return ( sigma_tau*GeVminus2_to_nb*1000. );
}
        

