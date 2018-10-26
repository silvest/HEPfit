/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBstrange.h"
#include "StandardModel.h"

double AFBstrange::computeThValue()
{
    return SM.AFB(SM.getQuarks(SM.STRANGE));
}

