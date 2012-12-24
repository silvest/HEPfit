/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha_2a.h"

double Alpha_2a::getThValue() { 
    double twoa = SM.getCKM().getAlpha()/M_PI*180.;
    return(twoa < 0. ? twoa + 180. : twoa);
}
