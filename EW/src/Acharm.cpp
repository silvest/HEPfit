/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Acharm.h"
#include "StandardModel.h"

double Acharm::computeThValue()
{
    return SM.A_f(SM.getQuarks(SM.CHARM));
}



