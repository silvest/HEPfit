/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rstrange.h"
#include "StandardModel.h"

double Rstrange::computeThValue()
{
    return SM.R0_f(SM.getQuarks(SM.STRANGE));
}


