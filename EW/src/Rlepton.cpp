/* 
 * File:   Rlepton.cpp
 * Author: mishima
 */

#include "Rlepton.h"


Rlepton::Rlepton(const EW& EW_i) : ThObservable(EW_i) {
    R0_l = EW_i.Gamma_had()/EW_i.Gamma_l(SM.ELECTRON);
}

double Rlepton::getThValue() {   
    return R0_l;
}
        

