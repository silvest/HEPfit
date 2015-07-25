/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFAFBcharmLEP2.h"


double ZFAFBcharmLEP2::computeThValue() { 
    double myZFAFBcharmLEP2, dummy;
    myZF.calcXS_AFB(6, sqrt_s, &dummy, &myZFAFBcharmLEP2);
    return myZFAFBcharmLEP2;
}


