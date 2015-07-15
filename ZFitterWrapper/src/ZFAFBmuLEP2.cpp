/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFAFBmuLEP2.h"


double ZFAFBmuLEP2::computeThValue() { 
    double myZFAFBmuLEP2, dummy;
    myZF.calcXS_AFB(2, sqrt_s, &dummy, &myZFAFBmuLEP2);
    return myZFAFBmuLEP2;
}


