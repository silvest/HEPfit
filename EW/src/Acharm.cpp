/* 
 * File:   Acharm.cpp
 * Author: mishima
 */

#include "Acharm.h"


Acharm::Acharm(const EW& EW_i) : ThObservable(EW_i) {
    A_c = EW_i.A_q(SM.CHARM);

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        
        A_c -= 48.0*alpha*s2*(3.0-4.0*s2)/pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
               *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double Acharm::getThValue() {   
    return A_c;
}
        


