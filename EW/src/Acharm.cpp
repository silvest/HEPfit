/* 
 * File:   Acharm.cpp
 * Author: mishima
 */

#include "Acharm.h"


double Acharm::getThValue() { 
    double A_c = myEW.A_q(SM.CHARM);

    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        double s4 = s2*s2;
        
        A_c -= 48.0*alpha*s2*(3.0-4.0*s2)/pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
               *( myEW.S() - 4.0*c2*s2*myEW.T() );
    }
    
    return A_c;
}
        


