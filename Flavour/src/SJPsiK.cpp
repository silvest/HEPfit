/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SJPsiK.h"
#include "StandardModel.h"
#include "AmpDB2.h"

double SJPsiK::computeThValue() 
{
    return sin(-SM.getFlavour().getDB2(0).getM21(FULLNLO).arg() - 2.*(SM.getCKM().computelamc_s()*SM.getCKM().computelamc()).arg() + 2.*SM.getPhiBd());
}
