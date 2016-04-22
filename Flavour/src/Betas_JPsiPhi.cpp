/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Betas_JPsiPhi.h"

double Betas_JPsiPhi::computeThValue() 
{
    return (remainder((AmpBs(FULLNLO).arg() + SM.getPhiBs() )/2.*180./M_PI,180.));
}
