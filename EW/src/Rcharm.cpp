/* 
 * File:   Rcharm.cpp
 * Author: mishima
 */

#include "Rcharm.h"


Rcharm::Rcharm(const EW& EW_i) : ThObservable(EW_i) {
    R0_c = EW_i.Gamma_q(SM.CHARM)/EW_i.Gamma_had();
}

double Rcharm::getThValue() {   
    return R0_c;
}
        

