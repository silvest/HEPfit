/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBlepton.h"


double AFBlepton::computeThValue() 
{
    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        AFB_l = SM.getMyEW_NPZff().AFBlepton(AFB_l);

    return SM.AFB(SM.getLeptons(SM.ELECTRON));
}
        

