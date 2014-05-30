/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rcharm.h"


double Rcharm::computeThValue() 
{   
    double R0_c = SM.R0_f(SM.getQuarks(SM.CHARM));

    return R0_c;
}
        

