/* 
 * File:   GammaW.cpp
 * Author: mishima
 */

#include "GammaW.h"


GammaW::GammaW(const EW& EW_i) : ThObservable(EW_i) {
    Gamma_W = EW_i.getSM().GammaW();
        
    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().alphaMz();
        double c = sqrt(EW_i.c2());
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
    
        Gamma_W -= 3.0*alpha*alpha*c*EW_i.getSM().getMz()/8.0/s2/(c2-s2)
                   *( EW_i.S() - 2.0*c2*EW_i.T() - (c2-s2)*EW_i.U()/2.0/s2 );
    }
}

double GammaW::getThValue() {   
    return Gamma_W;
}
        