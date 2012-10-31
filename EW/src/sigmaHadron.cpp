/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"


double sigmaHadron::getThValue() { 
    double sigma_had;
    if (myEWTYPE==EW::EWCHMN)  
        sigma_had = myEW.getMyEW_CHMN().sigma0_had();
    else if (myEWTYPE==EW::EWABC) 
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {   
        sigma_had = myEW.sigma0_had();
        
        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //sigma_had = 41.479/GeVminus2_to_nb;                
                
                double delta_l = - 0.000192*myEW.S() + 0.000790*myEW.T();
                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                double delta_Z = - 0.00961*myEW.S() + 0.0263*myEW.T();
                sigma_had *= 1.0 + delta_l/myEW.Gamma_l(SM.ELECTRON)
                             + delta_had/myEW.Gamma_had()
                             - 2.0*delta_Z/myEW.Gamma_Z();
            } else {
                double alpha = myEW.alpha0();  
                double Mz = SM.getMz();
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                double s6 = s4*s2;        
                double s8 = s6*s2;
                
                sigma_had -= 72.0*M_PI*alpha
                             *(729.0-4788.0*s2+8352.0*s4-6176.0*s6+640.0*s8)
                             /Mz/Mz/pow(63.0-120.0*s2+160.0*s4, 3.0)/(c2-s2)
                             *( myEW.S() - 4.0*c2*s2*myEW.T() );         
            }
        }
    }
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


