/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rbottom.h"


double Rbottom::computeThValue() 
{ 
    double R0_b = myEW.R0_b();

    /* NP contribution to the Zff vertex */
    if (myEW.checkNPZff_linearized() && SM.ModelName().compare("StandardModel") != 0)
        R0_b = myEW.getMyEW_NPZff().Rbottom(R0_b);

    return R0_b;
}
        

