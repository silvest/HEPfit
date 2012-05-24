/* 
 * File:   PtauPol.cpp
 * Author: mishima
 */

#include "PtauPol.h"


PtauPol::PtauPol(const EW& EW_i) : ThObservable(EW_i) {
    P_tau_pol = EW_i.A_l(SM.TAU);

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
    
        P_tau_pol -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
                     *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double PtauPol::getThValue() {   
    return P_tau_pol;
}
        
