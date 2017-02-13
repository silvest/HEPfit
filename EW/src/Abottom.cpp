/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Abottom.h"
#include "StandardModel.h"

double Abottom::computeThValue()
{
    return SM.A_f(SM.getQuarks(SM.BOTTOM));
}


