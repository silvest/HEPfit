/* 
 * File:   Abottom.cpp
 * Author: mishima
 */

#include "Abottom.h"


Abottom::Abottom(const EW& EW_i) : ThObservable(EW_i) {
    A_b = EW_i.A_q(SM.BOTTOM);
}

double Abottom::getThValue() {   
    return A_b;
}
        

