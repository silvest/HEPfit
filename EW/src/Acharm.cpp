/* 
 * File:   Acharm.cpp
 * Author: mishima
 */

#include "Acharm.h"


Acharm::Acharm(const EW& EW_i) : ThObservable(EW_i) {
    A_c = EW_i.A_q(SM.CHARM);
}

double Acharm::getThValue() {   
    return A_c;
}
        


