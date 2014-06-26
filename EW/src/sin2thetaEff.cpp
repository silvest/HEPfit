/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sin2thetaEff.h"

double sin2thetaEff::computeThValue()
{
    return SM.sin2thetaEff(SM.getLeptons(SM.ELECTRON));
}

