/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Acharm.h"


double Acharm::computeThValue() 
{ 
    double A_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        A_c = myEW.getMyEW_CHMN().A_q(SM.CHARM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        A_c = myEW.getMyEW_ABC().A_q(SM.CHARM,SM.epsilon1(),SM.epsilon3());
    else {
        A_c = myEW.A_q(SM.CHARM);

        if(myEWTYPE==EW::EWBURGESS) {
            double AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
            double delta_AFB_c = - 0.0147*SM.obliqueS() + 0.0104*SM.obliqueT();
            double delta_A_l = - 0.0284*SM.obliqueS() + 0.0201*SM.obliqueT();
            A_c *= 1.0 + delta_AFB_c/AFB_c - delta_A_l/myEW.A_l(SM.ELECTRON);
            return A_c;
        }
        
        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVf = SM.deltaGVq(SM.CHARM);
            double delGAf = SM.deltaGAq(SM.CHARM);
            if (delGVf!=0.0 || delGAf!=0.0) {
                double gVf = SM.StandardModel::gVq(SM.CHARM).real();
                double gAf = SM.StandardModel::gAq(SM.CHARM).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                A_c -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("Acharm::computeThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //A_c -= myEW.A_q(SM.CHARM);
    }
    
    return A_c;
}
        


