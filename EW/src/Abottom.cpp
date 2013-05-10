/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Abottom.h"


double Abottom::getThValue() 
{
    double A_b;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        A_b = myEW.getMyEW_CHMN().A_q(SM.BOTTOM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        A_b = myEW.getMyEW_ABC().A_b(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {
        A_b = myEW.A_q(SM.BOTTOM);

        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //A_b = 0.93464;
                
                double AFB_b = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
                double delta_AFB_b = - 0.0188*myEW.S() + 0.0131*myEW.T();
                double delta_A_l = - 0.0284*myEW.S() + 0.0201*myEW.T();
                A_b *= 1.0 + delta_AFB_b/AFB_b - delta_A_l/myEW.A_l(SM.ELECTRON);
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                
                A_b -= 12.0*alpha*s2*(3.0-2.0*s2)/pow(9.0-12.0*s2+8.0*s4, 2.0)/(c2-s2)
                       *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVf = SM.deltaGVq(SM.BOTTOM);
            double delGAf = SM.deltaGAq(SM.BOTTOM);
            if (delGVf!=0.0 || delGAf!=0.0) {
                double gVf = SM.StandardModel::gVq(SM.BOTTOM).real();
                double gAf = SM.StandardModel::gAq(SM.BOTTOM).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                A_b -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
            }
        }
        
        /* TEST */
        //A_b -= myEW.A_q(SM.BOTTOM);
    }
    
    return A_b;
}
        

