/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaZ.h"


double GammaZ::getThValue() { 
    double Gamma_Z;
    if (bCHMN)  
        Gamma_Z = myEW_CHMN.GammaZ();
    else {
        Gamma_Z = myEW.Gamma_Z();

        if ( myEW.checkModelForSTU() ) {
            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //Gamma_Z = 2.4954; 

                Gamma_Z += - 0.00961*myEW.S() + 0.0263*myEW.T();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                Gamma_Z += alpha*alpha*SM.getMz()/72.0/c2/s2/(c2-s2)
                           *( -10.0*(3.0-8.0*s2)*myEW.S() 
                              + (63.0-126.0*s2-40.0*s4)*myEW.T() );
            }
        }
    }
  
    return Gamma_Z;
}
        
