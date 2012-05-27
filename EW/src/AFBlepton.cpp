/* 
 * File:   AFBlepton.cpp
 * Author: mishima
 */

#include "AFBlepton.h"


AFBlepton::AFBlepton(const EW& EW_i) : ThObservable(EW_i) {
    AFB_l = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_l(SM.ELECTRON);

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().alphaMz();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;

        AFB_l -= 6.0*alpha*s2*(1.0-4.0*s2)/pow(1.0-4.0*s2+8.0*s4, 3.0)
               *( EW_i.S() - 4.0*c2*s2*EW_i.T() );        
    }
}

double AFBlepton::getThValue() {   
    return AFB_l;
}
        

