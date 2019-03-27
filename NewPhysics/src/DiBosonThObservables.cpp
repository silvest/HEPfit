/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DiBosonThObservables.h"
#include "NPbase.h"

mueeWW::mueeWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWW called with a class whose parent is not NPbase");
}

double mueeWW::computeThValue()
{
    return myNPbase->mueeWW(sqrt_s);
}

mueeWWPol::mueeWWPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWWPol called with a class whose parent is not NPbase");
}

double mueeWWPol::computeThValue()
{
    return myNPbase->mueeWWPol(sqrt_s, Pol_em, Pol_ep);
}