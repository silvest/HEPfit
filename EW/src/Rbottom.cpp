/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rbottom.h"

double Rbottom::computeThValue()
{
    return SM.R0_f(SM.getQuarks(SM.BOTTOM));
}


