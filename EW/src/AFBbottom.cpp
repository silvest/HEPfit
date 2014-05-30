/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBbottom.h"


double AFBbottom::computeThValue() 
{   
    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        AFB_b = SM.getMyEW_NPZff().AFBbottom(AFB_b);

    return SM.AFB(SM.getQuarks(SM.BOTTOM));
}
