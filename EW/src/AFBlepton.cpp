/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBlepton.h"


double AFBlepton::computeThValue() 
{   
    double AFB_l = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        AFB_l = myEW.getMyEW_NPZff().AFBlepton(AFB_l);

    return AFB_l;
}
        

