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
    double AFB_l = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

    /* NP contribution to the Zff vertex */
    if (myEW.checkLEP1NP())
        AFB_l = myEW.getMyEW_NPZff().AFBlepton(AFB_l);

    /* Debug: extract pure NP contribution */
    //AFB_l -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

    return AFB_l;
}
        

