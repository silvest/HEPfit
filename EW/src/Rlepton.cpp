/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rlepton.h"


double Rlepton::computeThValue() 
{
    double R0_l = myEW.R0_l();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        R0_l = myEW.getMyEW_NPZff().Rlepton(R0_l);
 
    return R0_l;
}
        

