/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBlepton.h"


double AFBlepton::getThValue() 
{   
    double AFB_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        AFB_l = myEW.getMyEW_CHMN().AFB_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC) 
        AFB_l = myEW.getMyEW_ABC().AFB_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double AFB_l_Born = 0.01696*(1.0 - 34.0*delta_als);
        AFB_l = AFB_l_Born*(1.0 + 34.72*SM.epsilon1() - 45.15*SM.epsilon3());
    } else {    
        AFB_l = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

        if(myEWTYPE==EW::EWBURGESS) {
            AFB_l += - 0.00677*SM.obliqueS() + 0.00479*SM.obliqueT();
            return AFB_l;
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            double delGVe = SM.deltaGVl(SM.ELECTRON);
            double delGAe = SM.deltaGAl(SM.ELECTRON);
            if (delGVe!=0.0 || delGAe!=0.0) {
                double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
                double Ge = gVe*gVe + gAe*gAe;
                double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;

                AFB_l -= 6.0*gVe*gAe*(gVe*gVe - gAe*gAe)*gAe*gAe/Ge/Ge/Ge*delGVeOverGAe;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("AFBlepton::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //AFB_l -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);
    }

    return AFB_l;
}
        

