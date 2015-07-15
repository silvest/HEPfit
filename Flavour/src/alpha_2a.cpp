/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha_2a.h"

double Alpha_2a::computeThValue() 
{ 
    double twoa = SM.computeAlpha()/M_PI*180.;
    return(twoa < 0. ? twoa + 180. : twoa);
}
