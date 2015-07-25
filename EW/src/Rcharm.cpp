/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rcharm.h"

double Rcharm::computeThValue()
{
    return SM.R0_f(SM.getQuarks(SM.CHARM));
}


