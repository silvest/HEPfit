/* 
 * File:   AFBcharm.cpp
 * Author: mishima
 */

#include "AFBcharm.h"


AFBcharm::AFBcharm(const EW& EW_i) : ThObservable(EW_i) {
    AFB_c = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_q(SM.CHARM);
}

double AFBcharm::getThValue() {   
    return AFB_c;
}
        
