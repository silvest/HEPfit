/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBbottom.h"


double AFBbottom::computeThValue() 
{   
    double AFB_b = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
              
    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        AFB_b = myEW.getMyEW_NPZff().AFBbottom(AFB_b);

    return AFB_b;
}
