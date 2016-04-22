/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"

double Alpha::computeThValue() 
{
    double alpha = (SM.computeAlpha() - SM.getPhiBd())/M_PI*180.;
    return(remainder(alpha,360.));
}
