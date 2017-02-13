/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaW.h"
#include "StandardModel.h"

double GammaW::computeThValue()
{
    return SM.GammaW();
}
