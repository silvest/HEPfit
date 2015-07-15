/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBcharm.h"

double AFBcharm::computeThValue()
{
    return SM.AFB(SM.getQuarks(SM.CHARM));
}

