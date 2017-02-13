/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"
#include "StandardModel.h"

double Alepton::computeThValue()
{
    return SM.A_f(SM.getLeptons(SM.ELECTRON));
}


