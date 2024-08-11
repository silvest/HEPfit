/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"
#include "StandardModel.h"
#include "AmpDB2.h"

Alpha::Alpha(const StandardModel& SM_i) : ThObservable(SM_i)
{        SM.getFlavour().getDB2(0);
}

double Alpha::computeThValue() 
{
    // alpha is really extracted as pi + 1/2 arg AmpDB2 - gamma 
    double alpha = (M_PI + SM.getFlavour().getDB2(0).getM21(FULLNLO).arg()/2. - SM.getCKM().computeGamma() - SM.getPhiBd())/M_PI*180.;
    return(remainder(alpha,360.));
}
