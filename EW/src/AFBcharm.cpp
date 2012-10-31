/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBcharm.h"


double AFBcharm::getThValue() {   
    double AFB_c;
    if (bCHMN)  
        AFB_c = myEW_CHMN.AFB_q(SM.CHARM);
    else {
        AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    
        if ( myEW.checkModelForSTU() ) {
            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //AFB_c = 0.0739;
                
                AFB_c += - 0.0147*myEW.S() + 0.0104*myEW.T();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                double s6 = s4*s2;        
                double s8 = s6*s2;
                
                AFB_c -= 9.0*alpha*s2*(39.0-310.0*s2+992.0*s4-1600.0*s6+1024.0*s8)
                         /pow(1.0-4.0*s2+8.0*s4, 2.0)
                         /pow(9.0-24.0*s2+32.0*s4, 2.0)/(c2-s2)
                         *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
    }
    
    return AFB_c;
}
        
