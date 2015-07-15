/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"

double Alpha::computeThValue() 
{
    return(SM.computeAlpha()/M_PI*180.);
}
