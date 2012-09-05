/* 
 * File:   sin2thetaEff.cpp
 * Author: mishima
 */

#include "sin2thetaEff.h"


double sin2thetaEff::getThValue() {  
    double sin2_theta_eff = myEW.sin2thetaEff(SM.ELECTRON);
    
    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        
        sin2_theta_eff += alpha/4.0/(c2-s2)
                          *( myEW.S() - 4.0*c2*s2*myEW.T() );
    }

    return sin2_theta_eff;
}

