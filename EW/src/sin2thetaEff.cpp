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
    
        if(myEWTYPE==EW::EWBURGESS) {
            sin2_theta_eff += 0.00362*SM.obliqueS() - 0.00256*SM.obliqueT();
            return sin2_theta_eff;
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
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("sin2thetaEff::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //sin2_theta_eff -= myEW.sin2thetaEff(SM.ELECTRON);
    }

    return sin2_theta_eff;
}

