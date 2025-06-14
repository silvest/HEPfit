/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "aTGC.h"
#include "StandardModel.h"

/* -------------------------------------*/

deltag1Z::deltag1Z(const StandardModel& SM_i, const double mu_i)
:ThObservable(SM_i), mu(mu_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("deltag1Z called with a class whose parent is not NPbase");
}


deltag1Z::~deltag1Z()
{}

double deltag1Z::computeThValue()
{
    return myNPbase->deltag1ZNP(mu);
}

/* -------------------------------------*/

deltag1gamma::deltag1gamma(const StandardModel& SM_i, const double mu_i)
:ThObservable(SM_i), mu(mu_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("deltag1gamma called with a class whose parent is not NPbase");
}


deltag1gamma::~deltag1gamma()
{}

double deltag1gamma::computeThValue()
{
    return myNPbase->deltag1gaNP(mu);
}

/* -------------------------------------*/

deltaKgamma::deltaKgamma(const StandardModel& SM_i, const double mu_i)
:ThObservable(SM_i), mu(mu_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("deltaKgamma called with a class whose parent is not NPbase");
}

deltaKgamma::~deltaKgamma()
{}

double deltaKgamma::computeThValue()
{
    return myNPbase->deltaKgammaNP(mu);
}

/* -------------------------------------*/

lambdaZ::lambdaZ(const StandardModel& SM_i, const double mu_i)
:ThObservable(SM_i), mu(mu_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("lambdaZ called with a class whose parent is not NPbase");
}

lambdaZ::~lambdaZ()
{}

double lambdaZ::computeThValue()
{
    return myNPbase->lambdaZNP(mu);
}


