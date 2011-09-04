/* 
 * File:   PtauPol.cpp
 * Author: mishima
 */

#include "PtauPol.h"


PtauPol::PtauPol(const EW& EW_i) : ThObservable(EW_i) {
    P_tau_pol = EW_i.A_l(SM.TAU);
}

double PtauPol::getThValue() {   
    return P_tau_pol;
}
        
