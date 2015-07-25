/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Betas_JPsiPhi.h"

double Betas_JPsiPhi::computeThValue() 
{
    return AmpBs(FULLNLO).arg()/2.*180./M_PI;
}
