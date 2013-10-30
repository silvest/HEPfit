/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"

double Alpha::computeThValue() {
    return(SM.getCKM().getAlpha()/M_PI*180.);
}
