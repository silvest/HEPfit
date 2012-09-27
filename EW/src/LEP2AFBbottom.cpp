/* 
 * File:   LEP2AFBbottom.cpp
 * Author: mishima
 */

#include "LEP2AFBbottom.h"


LEP2AFBbottom::LEP2AFBbottom(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2AFBbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double AFB_b = myEW.getSM().AFB_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ, 
                                           bDP, bWEAK, bQED);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        AFB_b += 0.0;       
    }
    
    return AFB_b;
}
        


