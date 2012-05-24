/* 
 * File:   GammaZ.cpp
 * Author: mishima
 */

#include "GammaZ.h"


GammaZ::GammaZ(const EW& EW_i) : ThObservable(EW_i) {
    Gamma_Z = EW_i.Gamma_Z();

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        
        Gamma_Z += alpha*alpha*EW_i.getSM().getMz()/72.0/c2/s2/(c2-s2)
                   *( -10.0*(3.0-8.0*s2)*EW_i.S() 
                      + (63.0-126.0*s2-40.0*s4)*EW_i.T() );
    }
}

double GammaZ::getThValue() {   
    return Gamma_Z;
}
        
