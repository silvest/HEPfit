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
    double AFB_b = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
              
    /* NP contribution to the Zff vertex */
    if (myEW.checkLEP1NP())
        AFB_b = myEW.getMyEW_NPZff().AFBbottom(AFB_b);

    /* Debug: extract pure NP contribution */
    //AFB_b -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);

    return AFB_b;
}
