/* 
 * File:   AFBlepton.cpp
 * Author: mishima
 */

#include "AFBlepton.h"


AFBlepton::AFBlepton(const EW& EW_i) : ThObservable(EW_i) {
    AFB_l = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_l(SM.ELECTRON);
}

double AFBlepton::getThValue() {   
    return AFB_l;
}
        

