/* 
 * File:   AFBbottom.cpp
 * Author: mishima
 */

#include "AFBbottom.h"


AFBbottom::AFBbottom(const EW& EW_i) : ThObservable(EW_i) {
    AFB_b = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_q(SM.BOTTOM);

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().alphaMz();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        double s6 = s4*s2;        
        double s8 = s6*s2;

        AFB_b -= 18.0*alpha*s2*(15.0-76.0*s2+152.0*s4-160.0*s6+64.0*s8)
                /pow(1.0-4.0*s2+8.0*s4, 2.0)
                /pow(9.0-12.0*s2+8.0*s4, 2.0)/(c2-s2)
                 *( EW_i.S() - 4.0*c2*s2*EW_i.T() ); 
    }
}

double AFBbottom::getThValue() {   
    return AFB_b;
}
        

