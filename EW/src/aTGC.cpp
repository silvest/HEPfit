/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "aTGC.h"
#include "StandardModel.h"

/* -------------------------------------*/

deltag1Z::deltag1Z(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltag1Z::~deltag1Z()
{}

double deltag1Z::computeThValue()
{
    return myNPbase->deltag1ZNP();
}

/* -------------------------------------*/

deltaKgamma::deltaKgamma(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}

deltaKgamma::~deltaKgamma()
{}

double deltaKgamma::computeThValue()
{
    return myNPbase->deltaKgammaNP();
}

/* -------------------------------------*/

lambdaZ::lambdaZ(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}

lambdaZ::~lambdaZ()
{}

double lambdaZ::computeThValue()
{
    return myNPbase->lambdaZNP();
}


