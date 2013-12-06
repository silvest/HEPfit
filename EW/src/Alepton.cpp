/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPZbbbar.h>
#include "Alepton.h"


double Alepton::computeThValue()
{
    double A_l = myEW.A_l(SM.ELECTRON);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff())
        A_l = myEW.getMyEW_NPZff().Alepton(A_l);

    /* Debug: extract pure NP contribution */
    //A_l -= myEW.A_l(SM.ELECTRON);

    return A_l;
}


