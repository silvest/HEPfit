/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBbottom.h"
#include "StandardModel.h"

double AFBbottom::computeThValue()
{
    return SM.AFB(SM.getQuarks(SM.BOTTOM));
}
