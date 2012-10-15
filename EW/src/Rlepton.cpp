/* 
 * File:   Rlepton.cpp
 * Author: mishima
 */

#include "Rlepton.h"


double Rlepton::getThValue() {
    double R0_l;
    if (bCHMN)  
        R0_l = myEW_CHMN.R_l(SM.ELECTRON);
    else {       
        R0_l = myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
        
        if ( myEW.checkModelForSTU() ) {
            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                R0_l = 20.740;

                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                double delta_l = - 0.000192*myEW.S() + 0.000790*myEW.T();
                R0_l *= 1.0 + delta_had/myEW.Gamma_had()
                        - delta_l/myEW.Gamma_l(SM.ELECTRON);
            } else {
                double alpha = SM.alphaMz();
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                
                R0_l += 8.0*alpha*(3.0-2.0*s2)*(1.0-5.0*s2)
                        /3.0/pow(1.0-4.0*s2+8.0*s4, 2.0)/(c2-s2)
                        *( myEW.S() - 4.0*c2*s2*myEW.T() );                
            }
        }
    }
 
    return R0_l;
}
        

