/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
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
        
        if (myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().AFBbottom(AFB_b);
              
        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            AFB_b = myEW.getMyEW_NPZff().AFBbottom(AFB_b);
        
        /* Debug: extract pure NP contribution */
        //AFB_b -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
    }

    return AFB_b;
}
