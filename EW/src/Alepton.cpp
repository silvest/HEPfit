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
            double alpha = SM.alphaMz();
            double c2 = myEW.c02();
            double s2 = myEW.s02();
            double s4 = s2*s2;
            
            A_l -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
                   *( myEW.S() - 4.0*c2*s2*myEW.T() );
        }
    }
    
    return A_l;
}
        

