/* 
 * File:   Rcharm.cpp
 * Author: mishima
 */

#include "Rcharm.h"


double Rcharm::getThValue() {   
    double R0_c;
    if (myEWTYPE==EW::EWCHMN)  
        R0_c = myEW.getMyEW_CHMN().R_c();
    else if (myEWTYPE==EW::EWABC) 
        R0_c = myEW.getMyEW_ABC().R_c(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {    
        R0_c = myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();
        
        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //R0_c = 0.17223;
                    
                double delta_c_over_Gamma_c = - 0.00649*myEW.S() + 0.0124*myEW.T();
                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                R0_c *= 1.0 + delta_c_over_Gamma_c - delta_had/myEW.Gamma_had();
            } else {            
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                R0_c -= 9.0*alpha*(9.0-36.0*s2+16.0*s4)
                        /pow(45.0-84.0*s2+88.0*s4, 2.0)/(c2-s2)
                        *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
    }

    return R0_c;
}
        

