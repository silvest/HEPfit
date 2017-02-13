/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rbottom.h"
#include "StandardModel.h"

double Rbottom::computeThValue()
{
    return SM.R0_f(SM.getQuarks(SM.BOTTOM));
}


