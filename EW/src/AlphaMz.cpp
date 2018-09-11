/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AlphaMz.h"
#include "StandardModel.h"

double AlphaEmMz::computeThValue()
{
    return SM.alphaMz();
}


