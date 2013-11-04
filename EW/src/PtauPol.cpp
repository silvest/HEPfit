/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PtauPol.h"
#include "EW_BURGESS.h"


double PtauPol::computeThValue() 
{  
    double P_tau_pol;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        P_tau_pol = myEW.getMyEW_CHMN().A_l(SM.TAU);    
    else if (myEWTYPE==EW::EWABC) 
        P_tau_pol = myEW.getMyEW_ABC().A_l(SM.TAU,SM.epsilon1(),SM.epsilon3(),false);
    else if (myEWTYPE==EW::EWABC2)
        P_tau_pol = myEW.getMyEW_ABC().A_l(SM.TAU,SM.epsilon1(),SM.epsilon3(),true);
    else {
        P_tau_pol = myEW.A_l(SM.TAU);

        if(myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().PtauPol(P_tau_pol);

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
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("PtauPol::computeThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //P_tau_pol -= myEW.A_l(SM.TAU);
    }
 
    return P_tau_pol;
}
        
