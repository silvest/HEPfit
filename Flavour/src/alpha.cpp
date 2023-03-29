/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"
#include "StandardModel.h"

Alpha::Alpha(const StandardModel& SM_i) : ThObservable(SM_i), AmpDB2(SM_i) 
{}

double Alpha::computeThValue() 
{
    // alpha is really extracted as pi + 1/2 arg AmpDB2 - gamma 
    double alpha = (M_PI + M12_Bd(FULLNLO).arg()/2. - SM.getCKM().computeGamma() - SM.getPhiBd())/M_PI*180.;
    return(remainder(alpha,360.));
}
