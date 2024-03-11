/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Phis_JPsiPhi.h"
#include "StandardModel.h"

double Phis_JPsiPhi::computeThValue() 
{
    return (-remainder((M21_Bs(FULLNLO).arg() - 2. * SM.getCKM().computelamc_s().arg() + 2.*SM.getPhiBs() ),2.*M_PI));
}
