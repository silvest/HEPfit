/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha.h"

double Alpha::getThValue() {
    return(SM.getCKM().getAlpha()/M_PI*180.);
}
