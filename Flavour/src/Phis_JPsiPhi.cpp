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
    return (remainder((-AmpBs(FULLNLO).arg() + 2.*SM.getPhiBs() ),2.*M_PI));
}
