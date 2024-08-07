/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "C2beta.h"
#include "StandardModel.h"
#include "AmpDB2.h"

C2beta::C2beta(const StandardModel& SM_i) : ThObservable(SM_i)
{
    SM.getFlavour().getDB2(0);
}

double C2beta::computeThValue() 
{
    return cos(-SM.getFlavour().getDB2(0).getM21(FULLNLO).arg() + 2.*SM.getPhiBd());
}
