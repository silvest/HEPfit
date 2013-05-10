/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Acharm.h"


double Acharm::getThValue() 
{ 
    double A_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        A_c = myEW.getMyEW_CHMN().A_q(SM.CHARM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        A_c = myEW.getMyEW_ABC().A_q(SM.CHARM,SM.epsilon1(),SM.epsilon3());
    else {
        A_c = myEW.A_q(SM.CHARM);

        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //A_c = 0.6680;

                double AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
                double delta_AFB_c = - 0.0147*myEW.S() + 0.0104*myEW.T();
                double delta_A_l = - 0.0284*myEW.S() + 0.0201*myEW.T();
                A_c *= 1.0 + delta_AFB_c/AFB_c - delta_A_l/myEW.A_l(SM.ELECTRON);
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                
                A_c -= 48.0*alpha*s2*(3.0-4.0*s2)/pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
                       *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
        
        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVf = SM.deltaGVq(SM.CHARM);
            double delGAf = SM.deltaGAq(SM.CHARM);
            if (delGVf!=0.0 || delGAf!=0.0) {
                double gVf = SM.StandardModel::gVq(SM.CHARM).real();
                double gAf = SM.StandardModel::gAq(SM.CHARM).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                A_c -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
            }
        }
    }
    
    return A_c;
}
        


