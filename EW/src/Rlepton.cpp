/* 
 * File:   Rlepton.cpp
 * Author: mishima
 */

#include "Rlepton.h"


Rlepton::Rlepton(const EW& EW_i) : ThObservable(EW_i) {
    R0_l = EW_i.Gamma_had()/EW_i.Gamma_l(SM.ELECTRON);
    
    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;

        R0_l += 8.0*alpha*(3.0-2.0*s2)*(1.0-5.0*s2)
                /3.0/pow(1.0-4.0*s2+8.0*s4, 2.0)/(c2-s2)
                *( EW_i.S() - 4.0*c2*s2*EW_i.T() );                
    }
}

double Rlepton::getThValue() {   
    return R0_l;
}
        

