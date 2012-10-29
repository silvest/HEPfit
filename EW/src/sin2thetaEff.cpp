/* 
 * File:   sin2thetaEff.cpp
 * Author: mishima
 */

#include "sin2thetaEff.h"


double sin2thetaEff::getThValue() { 
    double sin2_theta_eff;
    if (bCHMN)  
        sin2_theta_eff = myEW_CHMN.sin2thetaEff();
    else { 
        sin2_theta_eff = myEW.sin2thetaEff(SM.ELECTRON);
    
        if ( myEW.checkModelForSTU() ) {
            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //sin2_theta_eff = 0.23148;

                sin2_theta_eff += 0.00362*myEW.S() - 0.00256*myEW.T();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                
                sin2_theta_eff += alpha/4.0/(c2-s2)
                                  *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
    }

    return sin2_theta_eff;
}

