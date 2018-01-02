/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMEWPO.h"
#include "StandardModel.h"

Rb0GTHDM::Rb0GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Rb0GTHDM::computeThValue()
{
    double DeltaRb0=0.0;
    double Rb0SM=myGTHDM.R0_f(SM.getQuarks(SM.BOTTOM));
    return Rb0SM+DeltaRb0;
}
