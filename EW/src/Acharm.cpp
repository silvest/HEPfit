/* 
 * File:   Acharm.cpp
 * Author: mishima
 */

#include "Acharm.h"


double Acharm::getThValue() { 
    double A_c;
    if (bCHMN)  
        A_c = myEW_CHMN.A_q(SM.CHARM);
    else {
        A_c = myEW.A_q(SM.CHARM);

        if ( myEW.checkModelForSTU() ) {
            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //A_c = 0.6680;

                double AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
                double delta_AFB_c = - 0.0147*myEW.S() + 0.0104*myEW.T();
                double delta_A_l = - 0.0284*myEW.S() + 0.0201*myEW.T();
                A_c *= 1.0 + delta_AFB_c/AFB_c - delta_A_l/myEW.A_l(SM.ELECTRON);
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                A_c -= 48.0*alpha*s2*(3.0-4.0*s2)/pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
                       *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
    }
    
    return A_c;
}
        


