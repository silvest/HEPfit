/* 
 * File:   GammaW.cpp
 * Author: mishima
 */

#include "GammaW.h"


double GammaW::getThValue() {  
    double Gamma_W;
    if (bCHMN)  
        Gamma_W = myEW_CHMN.GammaW();
    else {
        Gamma_W = SM.GammaW();
        
        if ( myEW.checkModelForSTU() ) {
            double Wbar = 0.0;        
            if (SM.ModelName()=="NewPhysicsSTUVWXY") {
                Wbar = (myEW.V() - myEW.W())/SM.alphaMz();
            }

            if(bBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                Gamma_W = 2.091; 
                
                Gamma_W *= 1.0 - 0.00723*myEW.S() + 0.0111*myEW.T() + 0.00849*myEW.U() + 0.00781*Wbar;
            } else {
                double alpha = SM.alphaMz();
                double c = sqrt(myEW.c02());
                double c2 = myEW.c02();
                double s2 = myEW.s02();
            
                Gamma_W -= 3.0*alpha*alpha*c*SM.getMz()/8.0/s2/(c2-s2)
                           *( myEW.S() - 2.0*c2*myEW.T() - (c2-s2)*myEW.U()/2.0/s2 
                              - 2.0*(c2 - s2)*Wbar );
            }
        }
    }
 
    return Gamma_W;
}
