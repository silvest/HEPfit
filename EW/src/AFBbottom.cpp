/* 
 * File:   AFBbottom.cpp
 * Author: mishima
 */

#include "AFBbottom.h"


AFBbottom::AFBbottom(const EW& EW_i) : ThObservable(EW_i) {
    AFB_b = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_q(SM.BOTTOM);
}

double AFBbottom::getThValue() {   
    return AFB_b;
}
        

