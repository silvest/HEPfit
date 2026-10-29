/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "phid.h"
#include "StandardModel.h"
#include "AmpDB2.h"

phid::phid(const StandardModel& SM_i) : ThObservable(SM_i)
{
    SM.getFlavour().getDB2(0);
}

double phid::computeThValue() 
{
    // Compute the phase of the B_d mixing amplitude relative to the b->c\bar{c}s decay amplitude
    // The phase is given in radians
    return remainder(-SM.getFlavour().getDB2(0).getM21(FULLNLO).arg() + 2.*(SM.getCKM().computelamc_s()*SM.getCKM().computelamc()).arg() + 2.*SM.getPhiBd(), 2.*M_PI);
}
