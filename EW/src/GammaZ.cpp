/* 
 * File:   GammaZ.cpp
 * Author: mishima
 */

#include "GammaZ.h"


GammaZ::GammaZ(const EW& EW_i) : ThObservable(EW_i) {
    Gamma_Z = EW_i.Gamma_Z();
}

double GammaZ::getThValue() {   
    return Gamma_Z;
}
        
