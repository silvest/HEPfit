/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rinv.h"
#include "StandardModel.h"

double Rinv::computeThValue()
{
    return (SM.Gamma_inv())/(SM.GammaZ(SM.getLeptons(SM.ELECTRON)));
}

