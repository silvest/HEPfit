/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBlepton.h"


double AFBlepton::getThValue() 
{   
    double AFB_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        AFB_l = myEW.getMyEW_CHMN().AFB_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC) 
        AFB_l = myEW.getMyEW_ABC().AFB_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double AFB_l_Born = 0.01696*(1.0 - 34.0*delta_als);
        AFB_l = AFB_l_Born*(1.0 + 34.72*SM.epsilon1() - 45.15*SM.epsilon3());
    } else {    
        AFB_l = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //AFB_l = 0.01627;
            
                AFB_l += - 0.00677*myEW.S() + 0.00479*myEW.T();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                AFB_l -= 6.0*alpha*s2*(1.0-4.0*s2)/pow(1.0-4.0*s2+8.0*s4, 3.0)
                         *( myEW.S() - 4.0*c2*s2*myEW.T() );        
            } 
        }
    }
     
     return AFB_l;
}
        

