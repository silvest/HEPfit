/* 
 * File:   LEP2sigmaMu.cpp
 * Author: mishima
 */

#include "LEP2sigmaMu.h"


double LEP2sigmaMu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double sigma_mu = myEW.getSM().sigma_l_LEP2(StandardModel::MU, s, Mw, GammaZ);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        sigma_mu += 0.0;       
    }
    
    return ( sigma_mu*GeVminus2_to_nb*1000. );
}
        

