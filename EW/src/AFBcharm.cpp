/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBcharm.h"


double AFBcharm::computeThValue() 
{   
    double AFB_c = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
    
    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        AFB_c = myEW.getMyEW_NPZff().AFBcharm(AFB_c);
    
    return AFB_c;
}
        
