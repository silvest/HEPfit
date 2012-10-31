/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"


double Alepton::getThValue() { 
    double A_l;
    if (bCHMN)  
        A_l = myEW_CHMN.A_l(SM.ELECTRON);
    else {
        A_l = myEW.A_l(SM.ELECTRON);
        
        if ( myEW.checkModelForSTU() ) {
            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //A_l = 0.1473;
            
                A_l += - 0.0284*myEW.S() + 0.0201*myEW.T();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                A_l -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
                       *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
    }
    
    return A_l;
}
        

