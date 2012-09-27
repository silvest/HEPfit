/* 
 * File:   LEP2Rbottom.cpp
 * Author: mishima
 */

#include "LEP2Rbottom.h"


LEP2Rbottom::LEP2Rbottom(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2Rbottom::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double sigma_b, sigma_had;
    sigma_b = myEW.getSM().sigma_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ, 
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
    double R_b = sigma_b/sigma_had;
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        R_b += 0.0;       
    }
    
    return R_b;
}
        




