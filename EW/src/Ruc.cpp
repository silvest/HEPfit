/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Ruc.h"
#include "StandardModel.h"

double Ruc::computeThValue()
{
    return 0.5 * (SM.R0_f(SM.getQuarks(SM.UP)) + SM.R0_f(SM.getQuarks(SM.CHARM)) );
}


