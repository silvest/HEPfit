/* 
 * File:   Abottom.cpp
 * Author: mishima
 */

#include "Abottom.h"


Abottom::Abottom(const EW& EW_i) : ThObservable(EW_i) {
    A_b = EW_i.A_q(SM.BOTTOM);

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        
        A_b -= 12.0*alpha*s2*(3.0-2.0*s2)/pow(9.0-12.0*s2+8.0*s4, 2.0)/(c2-s2)
               *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double Abottom::getThValue() {   
    return A_b;
}
        

