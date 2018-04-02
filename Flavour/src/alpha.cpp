/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"
#include "StandardModel.h"

Alpha::Alpha(const StandardModel& SM_i) : ThObservable(SM_i) 
{}

double Alpha::computeThValue() 
{
    double alpha = (SM.getCKM().computeAlpha() - SM.getPhiBd())/M_PI*180.;
    return(remainder(alpha,360.));
}
