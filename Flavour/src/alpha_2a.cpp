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
    double a_2a = (SM.getCKM().computeAlpha()-SM.getPhiBd())/M_PI*180.;
    if (a_2a < 0.)
        a_2a += 180.;
    else if (a_2a > 180.)
        a_2a -= 180.;
    return(a_2a);
}
