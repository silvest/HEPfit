/* 
 * File:   Alepton.cpp
 * Author: mishima
 */

#include "Alepton.h"


Alepton::Alepton(const EW& EW_i) : ThObservable(EW_i) {
    A_l = EW_i.A_l(SM.ELECTRON);
}

double Alepton::getThValue() {   
    return A_l;
}
        

