/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rcharm.h"


double Rcharm::computeThValue() 
{   
    double R0_c = myEW.R0_c();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        R0_c = myEW.getMyEW_NPZff().Rcharm(R0_c);

    return R0_c;
}
        

