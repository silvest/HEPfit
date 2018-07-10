/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sin2thetaEff.h"
#include "StandardModel.h"

double sin2thetaEff::computeThValue()
{
    return SM.sin2thetaEff(SM.getLeptons(SM.ELECTRON));
}



double sin2thetaEffel::computeThValue()
{
    return SM.sin2thetaEff(SM.getLeptons(SM.ELECTRON));
}


double sin2thetaEffmu::computeThValue()
{
    return SM.sin2thetaEff(SM.getLeptons(SM.MU));
}
