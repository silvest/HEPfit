/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "AFBlepton.h"


double AFBlepton::computeThValue() 
{   
    double AFB_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        AFB_l = myEW.getMyEW_CHMN().AFB_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC) 
        AFB_l = myEW.getMyEW_ABC().AFB_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3(),false);
    else if (myEWTYPE==EW::EWABC2)
        AFB_l = myEW.getMyEW_ABC().AFB_l(SM.ELECTRON,SM.epsilon1(),SM.epsilon3(),true);
    else {    
        AFB_l = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

        if (myEWTYPE==EW::EWBURGESS)
            return myEW.getMyEW_BURGESS().AFBlepton(AFB_l);

        /* NP contribution to the Zff vertex */
        if (myEW.checkLEP1NP())
            AFB_l = myEW.getMyEW_NPZff().AFBlepton(AFB_l);

        /* Debug: extract pure NP contribution */
        //AFB_l -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);
    }

    return AFB_l;
}
        

