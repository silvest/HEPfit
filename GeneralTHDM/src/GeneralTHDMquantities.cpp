/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMquantities.h"

Re_sigma_u::Re_sigma_u(const StandardModel& SM_i)
: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double Re_sigma_u::computeThValue()
{
    return myGeneralTHDM->getabssigma_u()*cos(myGeneralTHDM->getphi_u());
}
