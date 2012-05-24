/* 
 * File:   Alepton.cpp
 * Author: mishima
 */

#include "Alepton.h"


Alepton::Alepton(const EW& EW_i) : ThObservable(EW_i) {
    A_l = EW_i.A_l(SM.ELECTRON);
        
    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        
        A_l -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
               *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double Alepton::getThValue() {   
    return A_l;
}
        

