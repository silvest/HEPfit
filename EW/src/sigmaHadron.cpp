/* 
 * File:   sigmaHadron.cpp
 * Author: mishima
 */

#include "sigmaHadron.h"


sigmaHadron::sigmaHadron(const EW& EW_i) : ThObservable(EW_i) {
    sigma_had = EW_i.sigma0_had();
}

double sigmaHadron::getThValue() {   
    return sigma_had;
}
        


