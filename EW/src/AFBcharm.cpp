/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "AFBcharm.h"


double AFBcharm::computeThValue() 
{   
    double AFB_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        AFB_c = myEW.getMyEW_CHMN().AFB_q(SM.CHARM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        AFB_c = myEW.getMyEW_ABC().AFB_c();
    else {
        AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    
        if (myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().AFBcharm(AFB_c);

        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            AFB_c = myEW.getMyEW_NPZff().AFBcharm(AFB_c);

        /* Debug: extract pure NP contribution */
        //AFB_c -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    }
    
    return AFB_c;
}
        
