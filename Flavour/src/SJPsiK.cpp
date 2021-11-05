/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SJPsiK.h"
#include "StandardModel.h"

double SJPsiK::computeThValue() 
{
    return sin(-AmpBd(FULLNLO).arg() - 2.*(SM.getCKM().computelamc_s()*SM.getCKM().computelamc()).arg() + 2.*SM.getPhiBd());
}
