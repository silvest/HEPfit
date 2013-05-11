/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"


double Alepton::getThValue()
{
    double A_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)
        A_l = myEW.getMyEW_CHMN().A_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC)
        A_l = myEW.getMyEW_ABC().A_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double x0 = 0.075619 - 1.32*delta_alpha;
        double x = x0*(1.0 + 17.6*SM.epsilon1() - 22.9*SM.epsilon3());
        A_l = 2.0*x/(1.0 + x*x);
    } else {
        A_l = myEW.A_l(SM.ELECTRON);

        if(myEWTYPE==EW::EWBURGESS) {
            A_l += - 0.0284*SM.obliqueS() + 0.0201*SM.obliqueT();
            return A_l;
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVf = SM.deltaGVl(SM.ELECTRON);
            double delGAf = SM.deltaGAl(SM.ELECTRON);

            /* Oblique corrections */
            delGVf += myEW.delGVl_oblique(SM.ELECTRON);
            delGAf += myEW.delGAl_oblique(SM.ELECTRON);

            if (delGVf!=0.0 || delGAf!=0.0) {
                double gVf = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAf = SM.StandardModel::gAl(SM.ELECTRON).real();
                double Gf = gVf*gVf + gAf*gAf;
                double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

                A_l -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("Alepton::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //A_l -= myEW.A_l(SM.ELECTRON);
    }

    return A_l;
}


