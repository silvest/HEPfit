/* 
 * File:   sin2thetaEff.cpp
 * Author: mishima
 */

#include "sin2thetaEff.h"


sin2thetaEff::sin2thetaEff(const EW& EW_i) : ThObservable(EW_i) {
    sin2_theta_eff = EW_i.sin2thetaEff(SM.ELECTRON);
    
    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().alphaMz();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        
        sin2_theta_eff += alpha/4.0/(c2-s2)
                          *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double sin2thetaEff::getThValue() {   
    return sin2_theta_eff;
}

