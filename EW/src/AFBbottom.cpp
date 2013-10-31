/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBbottom.h"


double AFBbottom::computeThValue() 
{   
    double AFB_b;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN) 
        AFB_b = myEW.getMyEW_CHMN().AFB_q(SM.BOTTOM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        AFB_b = myEW.getMyEW_ABC().AFB_b(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {
        AFB_b = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
        
        if(myEWTYPE==EW::EWBURGESS) {
            AFB_b += - 0.0188*SM.obliqueS() + 0.0131*SM.obliqueT();
            return AFB_b;
        }
              
        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVe = SM.deltaGVl(SM.ELECTRON);
            double delGAe = SM.deltaGAl(SM.ELECTRON);
            double delGVf = SM.deltaGVq(SM.BOTTOM);
            double delGAf = SM.deltaGAq(SM.BOTTOM);
            if (delGVe!=0.0 || delGAe!=0.0 || delGVf!=0.0 || delGAf!=0.0) {
                double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
                double Ge = gVe*gVe + gAe*gAe;
                double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;
                //
                double gVf = SM.StandardModel::gVq(SM.BOTTOM).real();
                double gAf = SM.StandardModel::gAq(SM.BOTTOM).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                AFB_b -= 3.0*gVf*gAf*(gVe*gVe - gAe*gAe)*gAe*gAe/Gf/Ge/Ge*delGVeOverGAe
                         + 3.0*gVe*gAe*(gVf*gVf - gAf*gAf)*gAf*gAf/Ge/Gf/Gf*delGVfOverGAf;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("AFBbottom::computeThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");
        
        /* Debug: extract pure NP contribution */
        //AFB_b -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
    }

    return AFB_b;
}
