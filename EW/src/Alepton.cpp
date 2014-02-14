/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"


double Alepton::computeThValue()
{
    double A_l = myEW.A_l(SM.ELECTRON);

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        A_l = myEW.getMyEW_NPZff().Alepton(A_l);

    return A_l;
}


