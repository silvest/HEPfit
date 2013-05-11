/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Mw.h"


double Mw::getThValue() 
{
    double myMw;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();
    
    if (myEWTYPE==EW::EWCHMN)  
        myMw = myEW.getMyEW_CHMN().Mw();
    else if (myEWTYPE==EW::EWABC)
        myMw = myEW.getMyEW_ABC().Mw(SM.epsilon1(),SM.epsilon2(),SM.epsilon3());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double cW2_Born = 0.768905*(1.0 - 0.40*delta_alpha);
        double cW2 = cW2_Born*(1.0 + 1.43*SM.epsilon1() - 1.00*SM.epsilon2() - 0.86*SM.epsilon3());
        myMw = sqrt(cW2)*SM.getMz();
    } else {
        myMw = SM.Mw();    

        if(myEWTYPE==EW::EWBURGESS) {
            myMw *= 1.0 - 0.00723/2.0*SM.obliqueS() + 0.0111/2.0*SM.obliqueT() + 0.00849/2.0*SM.obliqueU();
            return myMw;
        }

        if (!SM.IsFlagNotLinearizedNP() ) {
            double alpha = myEW.alpha();
            double c2 = myEW.cW2_SM();
            double s2 = myEW.sW2_SM();

            myMw *= 1.0 - alpha/4.0/(c2-s2)
                    *( SM.obliqueS() - 2.0*c2*SM.obliqueT() - (c2-s2)*SM.obliqueU()/2.0/s2 );
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("Mw::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //myMw -= SM.Mw();
    }
    
    return myMw;
}

