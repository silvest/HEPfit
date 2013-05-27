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
    } else
        myMw = SM.Mw();

    return myMw;
}

