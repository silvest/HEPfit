/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PtauPol.h"


double PtauPol::getThValue() 
{  
    double P_tau_pol;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        P_tau_pol = myEW.getMyEW_CHMN().A_l(SM.TAU);    
    else if (myEWTYPE==EW::EWABC) 
        P_tau_pol = myEW.getMyEW_ABC().A_l(SM.TAU,SM.epsilon1(),SM.epsilon3());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double x0 = 0.075619 - 1.32*delta_alpha;
        double x = x0*(1.0 + 17.6*SM.epsilon1() - 22.9*SM.epsilon3());
        P_tau_pol = 2.0*x/(1.0 + x*x);
    } else {
        P_tau_pol = myEW.A_l(SM.TAU);

        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //P_tau_pol = 0.1473;

                P_tau_pol += - 0.0284*myEW.S() + 0.0201*myEW.T();
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                
                P_tau_pol -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
                             *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVf = SM.deltaGVl(SM.TAU);
            double delGAf = SM.deltaGAl(SM.TAU);
            if (delGVf!=0.0 || delGAf!=0.0) {
                double gVf = SM.StandardModel::gVl(SM.TAU).real();
                double gAf = SM.StandardModel::gAl(SM.TAU).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                P_tau_pol -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
            }
        }
    }
 
    return P_tau_pol;
}
        
