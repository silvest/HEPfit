/* 
 * File:   Rbottom.cpp
 * Author: mishima
 */

#include "Rbottom.h"


Rbottom::Rbottom(const EW& EW_i) : ThObservable(EW_i) {
    R0_b = EW_i.Gamma_q(SM.BOTTOM)/EW_i.Gamma_had();
    
    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        
        R0_b += 6.0*alpha*(9.0-36.0*s2+16.0*s4)
                /pow(45.0-84.0*s2+88.0*s4, 2.0)/(c2-s2)
                 *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double Rbottom::getThValue() {   
    return R0_b;
}
        

