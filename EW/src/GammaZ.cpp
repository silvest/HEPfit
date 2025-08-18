/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaZ.h"
#include "StandardModel.h"

double GammaZ::computeThValue()
{
    return SM.Gamma_Z();
}



double GammaZee::computeThValue()
{
    return SM.GammaZ(SM.getLeptons(SM.ELECTRON));
}

double GammaZmumu::computeThValue()
{
    return SM.GammaZ(SM.getLeptons(SM.MU));
}

double GammaZtautau::computeThValue()
{
    return SM.GammaZ(SM.getLeptons(SM.TAU));
}

double GammaZuu::computeThValue()
{
    return SM.GammaZ(SM.getQuarks(SM.UP));
}

double GammaZcc::computeThValue()
{
    return SM.GammaZ(SM.getQuarks(SM.CHARM));
}

double GammaZdd::computeThValue()
{
    return SM.GammaZ(SM.getQuarks(SM.DOWN));
}

double GammaZss::computeThValue()
{
    return SM.GammaZ(SM.getQuarks(SM.STRANGE));
}

double GammaZbb::computeThValue()
{
    return SM.GammaZ(SM.getQuarks(SM.BOTTOM));
}

double GammaZinv::computeThValue()
{
    return SM.Gamma_inv();
}

