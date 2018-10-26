/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Astrange.h"
#include "StandardModel.h"

double Astrange::computeThValue()
{
    return SM.A_f(SM.getQuarks(SM.STRANGE));
}



