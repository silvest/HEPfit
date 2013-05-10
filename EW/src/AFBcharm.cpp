/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBcharm.h"


double AFBcharm::getThValue() 
{   
    double AFB_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        AFB_c = myEW.getMyEW_CHMN().AFB_q(SM.CHARM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        AFB_c = myEW.getMyEW_ABC().AFB_c(SM.epsilon1(),SM.epsilon3());
    else {
        AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    
        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //AFB_c = 0.0739;
                
                AFB_c += - 0.0147*myEW.S() + 0.0104*myEW.T();
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                double s6 = s4*s2;        
                double s8 = s6*s2;
                
                AFB_c -= 9.0*alpha*s2*(39.0-310.0*s2+992.0*s4-1600.0*s6+1024.0*s8)
                         /pow(1.0-4.0*s2+8.0*s4, 2.0)
                         /pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
                         *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVe = SM.deltaGVl(SM.ELECTRON);
            double delGAe = SM.deltaGAl(SM.ELECTRON);
            double delGVf = SM.deltaGVq(SM.BOTTOM);
            double delGAf = SM.deltaGAq(SM.BOTTOM);
            if (delGVe!=0.0 || delGAe!=0.0 || delGVf!=0.0 || delGAf!=0.0) {
                double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
                double Ge = gVe*gVe + gAe*gAe;
                double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;
                //
                double gVf = SM.StandardModel::gVq(SM.CHARM).real();
                double gAf = SM.StandardModel::gAq(SM.CHARM).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                AFB_c -= 3.0*gVf*gAf*(gVe*gVe - gAe*gAe)*gAe*gAe/Gf/Ge/Ge*delGVeOverGAe
                         + 3.0*gVe*gAe*(gVf*gVf - gAf*gAf)*gAf*gAf/Ge/Gf/Gf*delGVfOverGAf;
            }
        }
    }
    
    return AFB_c;
}
        
