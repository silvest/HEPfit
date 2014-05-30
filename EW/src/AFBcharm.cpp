/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBcharm.h"


double AFBcharm::computeThValue() 
{   
//    /* NP contribution to the Zff vertex */
//    if (SM.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
//        AFB_c = SM.getMyEW_NPZff().AFBcharm(AFB_c);
    
    return SM.AFB(SM.getQuarks(SM.CHARM));
}
        
