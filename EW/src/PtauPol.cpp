/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PtauPol.h"


double PtauPol::getThValue() {  
    double P_tau_pol;
    if (myEWTYPE==EW::EWCHMN)  
        P_tau_pol = myEW.getMyEW_CHMN().A_l(SM.TAU);    
    else if (myEWTYPE==EW::EWABC) 
        P_tau_pol = myEW.getMyEW_ABC().A_l(SM.TAU,SM.epsilon1(),SM.epsilon3());
    else {
        P_tau_pol = myEW.A_l(SM.TAU);

        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //P_tau_pol = 0.1473;

                P_tau_pol += - 0.0284*myEW.S() + 0.0201*myEW.T();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                P_tau_pol -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
                             *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
    }
 
    return P_tau_pol;
}
        
