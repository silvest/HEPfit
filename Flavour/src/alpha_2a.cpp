/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha_2a.h"
#include "StandardModel.h"

double Alpha_2a::computeThValue() 
{ 
    // alpha is really extracted as pi + 1/2 arg AmpDB2 - gamma 
    double a_2a = (M_PI + AmpBd(FULLNLO).arg()/2. - SM.getCKM().computeGamma() - SM.getPhiBd())/M_PI*180.;
    if (a_2a < 0.)
        a_2a += 180.;
    else if (a_2a > 180.)
        a_2a -= 180.;
    return(a_2a);
}
