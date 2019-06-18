/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaZhad.h"
#include "StandardModel.h"

double GammaZhad::computeThValue()
{
    return SM.Gamma_had();
}

