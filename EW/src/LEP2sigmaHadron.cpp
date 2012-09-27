/* 
 * File:   LEP2sigmaHadron.cpp
 * Author: mishima
 */

#include "LEP2sigmaHadron.h"


LEP2sigmaHadron::LEP2sigmaHadron(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {
    bDP = true;
    bWEAK = true;
    bQED = true;
}


double LEP2sigmaHadron::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = myEW.getSM().Mw(); 
    double GammaZ = myEW.Gamma_Z();

    double sigmaH = myEW.getSM().sigma_q_LEP2(StandardModel::UP, s, Mw, GammaZ,
                                              bDP, bWEAK, bQED)
                    + myEW.getSM().sigma_q_LEP2(StandardModel::DOWN, s, Mw, GammaZ,
                                                bDP, bWEAK, bQED)
                    + myEW.getSM().sigma_q_LEP2(StandardModel::CHARM, s, Mw, GammaZ,
                                                bDP, bWEAK, bQED)
                    + myEW.getSM().sigma_q_LEP2(StandardModel::STRANGE, s, Mw, GammaZ,
                                                bDP, bWEAK, bQED)
                    + myEW.getSM().sigma_q_LEP2(StandardModel::BOTTOM, s, Mw, GammaZ,
                                                bDP, bWEAK, bQED);
    
    if ( myEW.checkModelForSTU() ) {
        // write codes!!
        sigmaH += 0.0;       
    }
    
    return ( sigmaH*GeVminus2_to_nb*1000. );
}
        

