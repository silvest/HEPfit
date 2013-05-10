/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sin2thetaEff.h"


double sin2thetaEff::getThValue() 
{ 
    double sin2_theta_eff;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        sin2_theta_eff = myEW.getMyEW_CHMN().sin2thetaEff();
    else if (myEWTYPE==EW::EWABC) 
        sin2_theta_eff = myEW.getMyEW_ABC().sin2thetaEff(SM.epsilon1(),SM.epsilon3());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double x0 = 0.075619 - 1.32*delta_als;
        double x = x0*(1.0 + 17.6*SM.epsilon1() - 22.9*SM.epsilon3());
        sin2_theta_eff = (1.0 - x)/4.0;
    } else { 
        sin2_theta_eff = myEW.sin2thetaEff(SM.ELECTRON);
    
        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //sin2_theta_eff = 0.23148;

                sin2_theta_eff += 0.00362*myEW.S() - 0.00256*myEW.T();
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                
                sin2_theta_eff += alpha/4.0/(c2-s2)
                                  *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVf = SM.deltaGVl(SM.ELECTRON);
            double delGAf = SM.deltaGAl(SM.ELECTRON);
            if (delGVf!=0.0 || delGAf!=0.0) {
                double gVf = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAf = SM.StandardModel::gAl(SM.ELECTRON).real();
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                sin2_theta_eff -= delGVfOverGAf/4.0;
            }
        }
    }

    return sin2_theta_eff;
}

