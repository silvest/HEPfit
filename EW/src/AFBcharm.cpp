/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBcharm.h"


double AFBcharm::getThValue() 
{   
    double AFB_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        AFB_c = myEW.getMyEW_CHMN().AFB_q(SM.CHARM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        AFB_c = myEW.getMyEW_ABC().AFB_c(SM.epsilon1(),SM.epsilon3());
    else {
        AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    
        if(myEWTYPE==EW::EWBURGESS) {
            AFB_c += - 0.0147*SM.obliqueS() + 0.0104*SM.obliqueT();
            return AFB_c;
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVe = SM.deltaGVl(SM.ELECTRON);
            double delGAe = SM.deltaGAl(SM.ELECTRON);
            double delGVf = SM.deltaGVq(SM.CHARM);
            double delGAf = SM.deltaGAq(SM.CHARM);

            /* Oblique corrections */
            delGVe += myEW.delGVl_oblique(SM.ELECTRON);
            delGAe += myEW.delGAl_oblique(SM.ELECTRON);
            delGVf += myEW.delGVq_oblique(SM.CHARM);
            delGAf += myEW.delGAq_oblique(SM.CHARM);

            if (delGVe!=0.0 || delGAe!=0.0 || delGVf!=0.0 || delGAf!=0.0) {
                double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
                double Ge = gVe*gVe + gAe*gAe;
                double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;
                //
                double gVf = SM.StandardModel::gVq(SM.CHARM).real();
                double gAf = SM.StandardModel::gAq(SM.CHARM).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                AFB_c -= 3.0*gVf*gAf*(gVe*gVe - gAe*gAe)*gAe*gAe/Gf/Ge/Ge*delGVeOverGAe
                         + 3.0*gVe*gAe*(gVf*gVf - gAf*gAf)*gAf*gAf/Ge/Gf/Gf*delGVfOverGAf;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("AFBcharm::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //AFB_c -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    }
    
    return AFB_c;
}
        
