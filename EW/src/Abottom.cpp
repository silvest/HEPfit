/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Abottom.h"


double Abottom::getThValue() {
    double A_b;
    if (bCHMN)  
        A_b = myEW_CHMN.A_q(SM.BOTTOM);
    else {
        A_b = myEW.A_q(SM.BOTTOM);

        if ( myEW.checkModelForSTU() ) {
            double alpha = SM.alphaMz();
            double c2 = myEW.c02();
            double s2 = myEW.s02();
            double s4 = s2*s2;
            
            A_b -= 12.0*alpha*s2*(3.0-2.0*s2)/pow(9.0-12.0*s2+8.0*s4, 2.0)/(c2-s2)
                    *( myEW.S() - 4.0*c2*s2*myEW.T() );
        }
    }
   
    return A_b;
}
        

