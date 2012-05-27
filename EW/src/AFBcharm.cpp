/* 
 * File:   AFBcharm.cpp
 * Author: mishima
 */

#include "AFBcharm.h"


AFBcharm::AFBcharm(const EW& EW_i) : ThObservable(EW_i) {
    AFB_c = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_q(SM.CHARM);

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().alphaMz();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        double s6 = s4*s2;        
        double s8 = s6*s2;

        AFB_c -= 9.0*alpha*s2*(39.0-310.0*s2+992.0*s4-1600.0*s6+1024.0*s8)
                /pow(1.0-4.0*s2+8.0*s4, 2.0)
                /pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
                 *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double AFBcharm::getThValue() {   
    return AFB_c;
}
        
