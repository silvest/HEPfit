/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBlepton.h"
#include "StandardModel.h"

double AFBlepton::computeThValue()
{
    double AFBe = SM.AFB(SM.getLeptons(SM.ELECTRON));
    double AFBmu = SM.AFB(SM.getLeptons(SM.MU));
    double AFBtau = SM.AFB(SM.getLeptons(SM.TAU));  

    return (AFBe + AFBmu + AFBtau)/3.;
}

double AFBelectron::computeThValue()
{
    return SM.AFB(SM.getLeptons(SM.ELECTRON));
}

double AFBmuon::computeThValue()
{
    return SM.AFB(SM.getLeptons(SM.MU));
}

double AFBtau::computeThValue()
{
    return SM.AFB(SM.getLeptons(SM.TAU));
}


