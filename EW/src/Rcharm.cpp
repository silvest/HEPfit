/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rcharm.h"
#include "StandardModel.h"

double Rcharm::computeThValue()
{
    return SM.R0_f(SM.getQuarks(SM.CHARM));
}


