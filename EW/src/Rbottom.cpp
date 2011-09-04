/* 
 * File:   Rbottom.cpp
 * Author: mishima
 */

#include "Rbottom.h"


Rbottom::Rbottom(const EW& EW_i) : ThObservable(EW_i) {
    R0_b = EW_i.Gamma_q(SM.BOTTOM)/EW_i.Gamma_had();
}

double Rbottom::getThValue() {   
    return R0_b;
}
        

