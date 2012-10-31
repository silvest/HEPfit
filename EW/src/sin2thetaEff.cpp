/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sin2thetaEff.h"


double sin2thetaEff::getThValue() { 
    double sin2_theta_eff;
    if (bCHMN)  
        sin2_theta_eff = myEW_CHMN.sin2thetaEff();
    else { 
        sin2_theta_eff = myEW.sin2thetaEff(SM.ELECTRON);
    
        if ( myEW.checkModelForSTU() ) {
            double alpha = SM.alphaMz();
            double c2 = myEW.c02();
            double s2 = myEW.s02();
            
            sin2_theta_eff += alpha/4.0/(c2-s2)
                              *( myEW.S() - 4.0*c2*s2*myEW.T() );
        }
    }

    return sin2_theta_eff;
}

