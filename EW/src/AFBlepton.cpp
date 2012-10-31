/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBlepton.h"


double AFBlepton::getThValue() {   
    double AFB_l;
    if (bCHMN)  
        AFB_l = myEW_CHMN.AFB_l(SM.ELECTRON);
    else {    
        AFB_l = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

        if ( myEW.checkModelForSTU() ) {
            double alpha = SM.alphaMz();
            double c2 = myEW.c02();
            double s2 = myEW.s02();
            double s4 = s2*s2;
            
            AFB_l -= 6.0*alpha*s2*(1.0-4.0*s2)/pow(1.0-4.0*s2+8.0*s4, 3.0)
                     *( myEW.S() - 4.0*c2*s2*myEW.T() );        
        }
    }
     
     return AFB_l;
}
        

