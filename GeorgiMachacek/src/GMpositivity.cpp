/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GMpositivity.h"
#include "StandardModel.h"

GMpositivity1::GMpositivity1(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity1::computeThValue()
{
    double lambda2=0.;
    double lambda3=0.;

    return lambda2+lambda3;
}
